#include "simulation/barnes_hut_simulation_with_collisions.h"
#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include <omp.h>

void BarnesHutSimulationWithCollisions::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulationWithCollisions::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    // create QuadTree representing the current state of the universe
    BoundingBox this_bb = universe.parallel_cpu_get_bounding_box();
    Quadtree quadtree(universe, this_bb, int8_t(2));    

    // calculate cumulative masses and centers of mass
    quadtree.calculate_cumulative_masses();
    quadtree.calculate_center_of_mass();

    // calculate forces
    calculate_forces(universe, quadtree);

    // calculate velocities & positions
    NaiveParallelSimulation::calculate_velocities(universe);
    NaiveParallelSimulation::calculate_positions(universe);

    // find collisions
    find_collisions(universe);

    universe.current_simulation_epoch += 1;
    if (create_intermediate_plots &&
        universe.current_simulation_epoch % plot_intermediate_epochs == 0)
    {
        plotter.add_bodies_to_image(universe);
        plotter.write_and_clear();
    }
    
    return;
}


struct PairComparator {
    const Universe* uni;
    PairComparator(const Universe* u) : uni(u) {}
    bool operator()(const std::pair<int,int>& lhs, const std::pair<int,int>& rhs) const {
        if (lhs.first == rhs.first) return uni->weights[lhs.second] > uni->weights[rhs.second];
        return uni->weights[lhs.first] > uni->weights[rhs.first];
    }
};

void BarnesHutSimulationWithCollisions::find_collisions(Universe& universe){
    // set of pairs with custom ordering: pairs which contain the heaviest bodies appear first
    std::set<std::pair<int, int>, PairComparator> pair_set{PairComparator(&universe)};

    // insert indizes of colliding pairs into map. after this, themaps
    for (int i = 0; i < universe.num_bodies; i++) 
        for (int k = 0; k < universe.num_bodies; k++) {
            // skip loop if on the same body or if value already tried
            if (i == k) continue;
                
            // calculate distance
            Vector2d<double> direction_vector = universe.positions[i] - universe.positions[k];
            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));
            if (distance < 100000000000) {
                // insert into pair set (which automatically sorts inserted elements by key)
                std::pair idx_pair = std::make_pair(i, k);
                pair_set.insert(idx_pair);
            }
        }

    // ordered set of indizes of bodies that have collided
    std::set<int> collided_indices;

    // access first (heaviest) element of the set of pairs
    while (!pair_set.empty()) {
        // Retrieve the heaviest pair (based on the custom comparator)
        auto heaviest_pair_it = pair_set.begin();
        std::pair<int, int> heaviest_pair = *heaviest_pair_it;

        collided_indices.insert(heaviest_pair.second);

        double old_heavy_mass = universe.weights[heaviest_pair.first];
        double old_light_mass = universe.weights[heaviest_pair.second];

        // update weight & velocity of heavier body
        universe.weights[heaviest_pair.first] = old_heavy_mass + universe.weights[heaviest_pair.second];
        universe.velocities[heaviest_pair.first] = ((universe.velocities[heaviest_pair.second]*old_light_mass) + (universe.velocities[heaviest_pair.first]*old_heavy_mass)) / universe.weights[heaviest_pair.first];

        // Remove the heaviest pair from the set    
        pair_set.erase(heaviest_pair_it);

        // Remove all pairs where the first element is equal to second_body
        std::erase_if(pair_set, [&](const std::pair<int, int>& p) -> bool {
            return p.first == heaviest_pair.second;
        });
    }
    
    // delete collided indizes
    int i = 0;
    for (int idx : collided_indices) {
        universe.positions.erase(universe.positions.begin() + idx - i);
        universe.forces.erase(universe.forces.begin() + idx - i);
        universe.weights.erase(universe.weights.begin() + idx - i);
        universe.velocities.erase(universe.velocities.begin() + idx - i);
        
        // Decrement no. of bodies in universe
        universe.num_bodies--;
        i++;
    }
}

void BarnesHutSimulationWithCollisions::find_collisions_parallel(Universe& universe){
      // set of pairs with custom ordering: pairs which contain the heaviest bodies appear first
    std::set<std::pair<int, int>, PairComparator> pair_set{PairComparator(&universe)};

    // Parallelize insertion of indizes of colliding pairs into pair set. After this step, all possible colliding pairs are in the pair set
    #pragma omp parallel
    {
        std::set<std::pair<int, int>, PairComparator> private_pair_set{PairComparator(&universe)};

        #pragma omp for nowait
        for (int i = 0; i < universe.num_bodies; i++) {
            for (int k = 0; k < universe.num_bodies; k++) {
                if (i == k) continue;

                Vector2d<double> direction_vector = universe.positions[i] - universe.positions[k];
                double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));
                if (distance < 100000000000) {
                    private_pair_set.emplace(std::make_pair(i, k));
                }
            }
        }

        #pragma omp critical
        {
            pair_set.insert(private_pair_set.begin(), private_pair_set.end());
        }
    }

    // Ordered set of indizes of bodies that have collided
    std::set<int> collided_indices;

    // access first (heaviest) element of the set of pairs
    while (!pair_set.empty()) {
        // Retrieve the heaviest pair (based on the custom comparator)
        auto heaviest_pair_it = pair_set.begin();
        std::pair<int, int> heaviest_pair = *heaviest_pair_it;

        collided_indices.insert(heaviest_pair.second);

        double old_heavy_mass = universe.weights[heaviest_pair.first];
        double old_light_mass = universe.weights[heaviest_pair.second];

        // update weight & velocity of heavier body (those of lighter body are irrelevant as the lighter body will be removed in the following step)
        universe.weights[heaviest_pair.first] = old_heavy_mass + universe.weights[heaviest_pair.second];
        universe.velocities[heaviest_pair.first] = ((universe.velocities[heaviest_pair.second]*old_light_mass) + (universe.velocities[heaviest_pair.first]*old_heavy_mass)) / universe.weights[heaviest_pair.first];

        // Remove the heaviest pair from the set    
        pair_set.erase(heaviest_pair_it);

        // Remove all pairs where the first element is equal to second_body (since the second body will have already collided)
        std::erase_if(pair_set, [&](const std::pair<int, int>& p) -> bool {
            return p.first == heaviest_pair.second;
        });
    }
    
    // Convert collided set to vector (to facilitate for-loop)
    std::vector<int> sorted_collided_indices(collided_indices.begin(), collided_indices.end());

    // Parallelize the deletion of collided indices
    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(sorted_collided_indices.size()); i++) {
        int idx = sorted_collided_indices[i];
        // update values
        universe.forces.erase(universe.forces.begin() + idx - i);
        universe.positions.erase(universe.positions.begin() + idx - i);
        universe.weights.erase(universe.weights.begin() + idx - i);
        universe.velocities.erase(universe.velocities.begin() + idx - i);
        // Decrement no. of bodies in universe
        universe.num_bodies--;
    }
}