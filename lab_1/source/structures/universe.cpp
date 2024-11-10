#include "structures/universe.h"
#include <limits>

BoundingBox Universe::get_bounding_box() {
  // Gebe kleinste BoundingBox, falls keinen Himmelskörper vorhanden
  if (num_bodies == 0)
    return BoundingBox(0, 0, 0, 0);

  // Setze zuerst niedrigste bzw. höchste Werte für x- und y-{min, max}
  double x_max = std::numeric_limits<double>::lowest();
  double y_max = std::numeric_limits<double>::lowest();

  double x_min = std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max();

  // Gehe über alle Positionen durch, aktualisiere x- und y-{min, max}
  for (Vector2d<double> position : positions) {
    if (position[0] > x_max)
      x_max = position[0];
    if (position[1] > y_max)
      y_max = position[1];

    if (position[0] < x_min)
      x_min = position[0];
    if (position[1] < y_min)
      y_min = position[1];
  }

  return BoundingBox(x_min, x_max, y_min, y_max);
}
