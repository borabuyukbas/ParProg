#include "quadtreeNode.h"
#include "structures/vector2d.h"

#include <iostream>


double QuadtreeNode::calculate_node_cumulative_mass(){
    // return if cumulative mass is ready
    if (this->cumulative_mass_ready)
    {
        return this->cumulative_mass;
    }

    // calculate cumulative mass sum of each child
    double sum = 0;
    #pragma omp parallel for reduction(+:sum)
    for (QuadtreeNode* child : children)
    {
        sum += child->calculate_node_cumulative_mass();
    }
    
    // set the cumulative mass ready for next uses and return
    this->cumulative_mass = sum;
    this->cumulative_mass_ready = true;
    return this->cumulative_mass;
}

QuadtreeNode::QuadtreeNode(BoundingBox arg_bounding_box){
    bounding_box = arg_bounding_box;
}

QuadtreeNode::~QuadtreeNode(){
    for (auto c : children) {
        delete c;
    }
}

Vector2d<double> QuadtreeNode::calculate_node_center_of_mass(){
    // return if center of mass is ready
    if (this->center_of_mass_ready)
    {
        return this->center_of_mass;
    }

    // calculate center of mass sum of each child
    // divided the vector2d into two doubles for easy reduction
    double numerator_x = 0;
    double numerator_y = 0;
    #pragma omp parallel for reduction(+:numerator_x) reduction(+:numerator_y)
    for (QuadtreeNode* child : children)
    {
        Vector2d<double> calculated_numerator = child->calculate_node_center_of_mass() * child->calculate_node_cumulative_mass();
        numerator_x += calculated_numerator[0];
        numerator_y += calculated_numerator[1];
    }
    
    // set the center of mass ready for next uses and return
    this->center_of_mass = Vector2d<double>(numerator_x, numerator_y) / this->calculate_node_cumulative_mass();
    this->center_of_mass_ready = true;
    return this->center_of_mass;
}