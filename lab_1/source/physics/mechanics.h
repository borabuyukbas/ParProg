#pragma once

#include "structures/vector2d.h"
#include <stdexcept>

template <typename T>
static Vector2d<T> calculate_acceleration(Vector2d<T> F, double m) {
  if (m == 0) {
    throw std::invalid_argument("Die Masse m darf nicht null sein.");
  }

  return F / m;
}

template <typename T>
static Vector2d<T> calculate_velocity(Vector2d<T> v0, Vector2d<T> a, double t) {
  return v0 + (a * t);
}