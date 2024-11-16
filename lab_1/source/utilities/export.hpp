#pragma once

#include "structures/universe.h"
#include <filesystem>
#include <format>
#include <fstream>

static void save_universe(std::filesystem::path save_universe_path,
                          Universe &universe) {

  std::ofstream universe_file(save_universe_path);
  std::string line;

  if (!universe_file.is_open()) {
    throw std::invalid_argument("Could not save universe to given path!");
  }

  universe_file << "### Bodies" << std::endl;
  universe_file << std::to_string(universe.num_bodies) << std::endl;

  universe_file << "### Positions" << std::endl;
  for (Vector2d<double> position : universe.positions) {
    universe_file << std::format("{:.6f}", position[0]) << " "
                  << std::format("{:.6f}", position[1]) << std::endl;
  }

  universe_file << "### Weights" << std::endl;
  for (double weight : universe.weights) {
    universe_file << std::format("{:.6f}", weight) << std::endl;
  }

  universe_file << "### Velocities" << std::endl;
  for (Vector2d<double> velocity : universe.velocities) {
    universe_file << std::format("{:.6f}", velocity[0]) << " "
                  << std::format("{:.6f}", velocity[1]) << std::endl;
  }

  universe_file << "### Forces" << std::endl;
  for (Vector2d<double> force : universe.forces) {
    universe_file << std::format("{:.6f}", force[0]) << " "
                  << std::format("{:.6f}", force[1]) << std::endl;
  }

  universe_file.close();
}
