#include "simulation/naive_sequential_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"
#include "simulation/constants.h"

#include <cmath>

// calculates and sets a 2D force vector ([F] = N) for every body contained in
// the referenced Universe object.
void NaiveSequentialSimulation::calculate_forces(Universe &universe) {
  for (int i = 0; i < universe.num_bodies; i++) {
    Vector2d<double> sum_forces(0.0, 0.0);
    for (int j = 0; j < universe.num_bodies; j++) {
      // don't add forces between the same mass
      if (j == i)
        continue;

      // calculate norm of vector between two masses (distance)
      // vector points from m_i to m_j ("j pulls on i")
      Vector2d<double> p_ij = universe.positions[j] - universe.positions[i];
      double distance = sqrt(pow(p_ij[0], 2.0) + pow(p_ij[1], 2.0));

      // calculate angle in radian (atan2 to guarantee full range of angles)
      double phi = atan2(p_ij[1], p_ij[0]);

      // absolute force on i
      double force = gravitational_force(universe.weights[i],
                                         universe.weights[j], distance);

      // get force components in x and y directions
      double force_x = force * cos(phi);
      double force_y = force * sin(phi);

      sum_forces = sum_forces + Vector2d<double>(force_x, force_y);
    }
    universe.forces[i] = sum_forces;
  }
}

void NaiveSequentialSimulation::calculate_velocities(Universe &universe) {
  for (int i = 0; i < universe.num_bodies; i++) {
    Vector2d<double> acceleration =
        calculate_acceleration(universe.forces[i], universe.weights[i]);

    universe.velocities[i] = calculate_velocity(universe.velocities[i],
                                                acceleration, epoch_in_seconds);
  }
}

void NaiveSequentialSimulation::calculate_positions(Universe &universe) {
  for (int i = 0; i < universe.num_bodies; i++) {
    Vector2d<double> shift = universe.velocities[i] * epoch_in_seconds;

    universe.positions[i] = universe.positions[i] + shift;
  }
}

void NaiveSequentialSimulation::simulate_epoch(
    Plotter &plotter, Universe &universe, bool create_intermediate_plots,
    std::uint32_t plot_intermediate_epochs) {
  calculate_forces(universe);
  calculate_velocities(universe);
  calculate_positions(universe);
  universe.current_simulation_epoch += 1;

  if (create_intermediate_plots &&
      universe.current_simulation_epoch % plot_intermediate_epochs == 0) {
    plotter.add_bodies_to_image(universe);
    plotter.write_and_clear();
  }
}

void NaiveSequentialSimulation::simulate_epochs(
    Plotter &plotter, Universe &universe, std::uint32_t num_epochs,
    bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs) {
  for (int i = 0; i < num_epochs; i++) {
    simulate_epoch(plotter, universe, create_intermediate_plots,
                   plot_intermediate_epochs);
  }
}