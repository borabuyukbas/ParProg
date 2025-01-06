#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"

#include <cmath>

void BarnesHutSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
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

    universe.current_simulation_epoch += 1;
    if (create_intermediate_plots &&
        universe.current_simulation_epoch % plot_intermediate_epochs == 0)
    {
        plotter.add_bodies_to_image(universe);
        plotter.write_and_clear();
    }
    
    return;
}

void BarnesHutSimulation::get_relevant_nodes(Universe& universe, Quadtree& quadtree, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double threshold_theta){
    // masses will have been set before this 
   quadtree.root->add_relevant_child(universe, relevant_nodes, body_position, body_index, threshold_theta);
}

void BarnesHutSimulation::calculate_forces(Universe& universe, Quadtree& quadtree){
    return;
}