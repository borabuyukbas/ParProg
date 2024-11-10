#pragma once
#include <cstdint>
#include <vector>

#include "structures/bounding_box.h"
#include "structures/vector2d.h"

class Universe {
public:
  std::uint32_t num_bodies; // Gesamtanzahl der Himmelskörper
  std::uint32_t
      current_simulation_epoch; // Anzahl der bereits simulierten Epochen
  std::vector<double> weights;  // Massen der Himmelskörper in kg
  std::vector<Vector2d<double>>
      forces; // Momentan auf die Körper wirkende Kräfte in N
  std::vector<Vector2d<double>>
      velocities; // Momentane Geschwindigkeiten der Körper in m/s
  std::vector<Vector2d<double>>
      positions; // Momentane Positionen der Körper in m

  Universe() : num_bodies(0), current_simulation_epoch(0) {}

  BoundingBox get_bounding_box();
};