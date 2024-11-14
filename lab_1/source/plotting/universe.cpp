#include "plotting/plotter.h"

void Plotter::add_bodies_to_image(Universe& universe) {
    // mark position of each orb
    for (Vector2d<double> position : universe.positions) {
        mark_position(position, 255, 255, 255);
    }
}