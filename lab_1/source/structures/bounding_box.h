#pragma once

#include <string>
#include "vector2d.h"

class BoundingBox
{
public:

    double x_min, x_max, y_min, y_max;
    BoundingBox() :
        x_min(0),
        x_max(0),
        y_min(0),
        y_max(0)
    {}

    BoundingBox(double p_x_min, double p_x_max, double p_y_min, double p_y_max) :
        x_min(p_x_min), x_max(p_x_max), y_min(p_y_min), y_max(p_y_max)
    {}

    // member function declarations
    std::string get_string();
    double get_diagonal();
    void plotting_sanity_check();
    BoundingBox get_scaled(std::uint32_t scaling_factor);
    bool contains(Vector2d<double> point);
    BoundingBox get_quadrant(std::uint8_t index);
};