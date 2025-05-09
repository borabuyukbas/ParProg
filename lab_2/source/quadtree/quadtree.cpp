#include "quadtree.h"

#include "quadtreeNode.h"
#include "structures/bounding_box.h"
#include "structures/universe.h"
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <vector>

#define CUTOFF 20000

std::vector<std::int32_t> get_indices(Universe& universe, BoundingBox bounding_box, std::vector<std::int32_t> body_indices) {
    // create a vector that holds incides for current quadrant
    std::vector<std::int32_t> indices;

    for (std::int32_t indice : body_indices) {
        if (bounding_box.contains(universe.positions[indice])) {
            indices.push_back(indice);
        }
    }
    
    return indices;
}

Quadtree::Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode){
    if (!(construct_mode == 0 || construct_mode == 1 || construct_mode == 2)) {
        throw std::runtime_error("Invalid 'construct_mode' has been passed to the Quadtree constructor!");
    }

    root = new QuadtreeNode(bounding_box);

    // construct body_indices vector with all bodies inside the bounding_box
    std::vector<std::int32_t> body_indices;
    for (int i = 0; i < universe.num_bodies; ++i) {
        if (bounding_box.contains(universe.positions[i])) {
            body_indices.push_back(i);
        }
    }

    switch (construct_mode) {
        case 0:
            root->children = construct(universe, bounding_box, body_indices);
            break;
        case 1:
            root->children = construct_task(universe, bounding_box, body_indices);
            break;
        case 2:
            root->children = construct_task_with_cutoff(universe, bounding_box, body_indices);
            break;
    }
}

Quadtree::~Quadtree(){
    delete root;
}

void Quadtree::calculate_cumulative_masses(){
    root->calculate_node_cumulative_mass();
}

void Quadtree::calculate_center_of_mass(){
    root->calculate_node_center_of_mass();
}

std::vector<QuadtreeNode*> Quadtree::construct(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    // initialize a new vector that contains subquadrants
    std::vector<QuadtreeNode*> subquadrants;

    // loop over all quadrants
    for (int i = 0; i <= 3; ++i) {
        BoundingBox quadrant = BB.get_quadrant(i); // quadrant bounding box
        
        // create a vector that holds incides for current quadrant
        std::vector<std::int32_t> indices = get_indices(universe, quadrant, body_indices);

        int size = indices.size();
        // if there is only one body, set the body identifier with the mass values and stop the recursion
        if (size == 1) {
            QuadtreeNode* node = new QuadtreeNode(quadrant);
            node->body_identifier = indices[0];
            node->cumulative_mass = universe.weights[indices[0]];
            node->cumulative_mass_ready = true;
            node->center_of_mass = universe.positions[indices[0]];
            node->center_of_mass_ready = true;
            subquadrants.push_back(node); // append the node in the subquadrants list
        }
        // if there are more then one bodies, search for the new children and copy the values back to the node
        else if (size > 1) {
            QuadtreeNode* node = new QuadtreeNode(quadrant);
            node->children = construct(universe, quadrant, indices);
            subquadrants.push_back(node); // append the node in the subquadrants list
        }
    }

    return subquadrants;
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    // initialize a new vector that contains subquadrants
    std::vector<QuadtreeNode*> subquadrants;

    #pragma omp parallel
    {
        #pragma omp single
        {
            // loop over all quadrants
            for (int i = 0; i <= 3; ++i) {
                BoundingBox quadrant = BB.get_quadrant(i); // quadrant bounding box
                
                // create a vector that holds incides for current quadrant
                std::vector<std::int32_t> indices = get_indices(universe, quadrant, body_indices);

                int size = indices.size();

                // if there is only one body, set the body identifier with the mass values and stop the recursion
                if (size == 1) {
                    QuadtreeNode* node = new QuadtreeNode(quadrant);
                    node->body_identifier = indices[0];
                    node->cumulative_mass = universe.weights[indices[0]];
                    node->cumulative_mass_ready = true;
                    node->center_of_mass = universe.positions[indices[0]];
                    node->center_of_mass_ready = true;

                    #pragma omp critical
                    subquadrants.push_back(node); // append the node in the subquadrants list
                }
                // if there are more then one bodies, search for the new children and copy the values back to the node
                else if (size > 1) {
                    #pragma omp task
                    {
                        QuadtreeNode* node = new QuadtreeNode(quadrant);
                        node->children = construct_task(universe, quadrant, indices);

                        #pragma omp critical
                        subquadrants.push_back(node); // append the node in the subquadrants list
                    }
                }
            }
        }
    }

    return subquadrants;
}

std::vector<QuadtreeNode*> Quadtree::construct_task_with_cutoff(Universe& universe, BoundingBox& BB, std::vector<std::int32_t>& body_indices){
    // initialize a new vector that contains subquadrants
    std::vector<QuadtreeNode*> subquadrants;

    #pragma omp parallel
    {
        #pragma omp single
        {
            // loop over all quadrants
            for (int i = 0; i <= 3; ++i) {
                BoundingBox quadrant = BB.get_quadrant(i); // quadrant bounding box
                
                // create a vector that holds incides for current quadrant
                std::vector<std::int32_t> indices = get_indices(universe, quadrant, body_indices);

                int size = indices.size();

                // if there is only one body, set the body identifier with the mass values and stop the recursion
                if (size == 1) {
                    QuadtreeNode* node = new QuadtreeNode(quadrant);
                    node->body_identifier = indices[0];
                    node->cumulative_mass = universe.weights[indices[0]];
                    node->cumulative_mass_ready = true;
                    node->center_of_mass = universe.positions[indices[0]];
                    node->center_of_mass_ready = true;

                    #pragma omp critical
                    subquadrants.push_back(node); // append the node in the subquadrants list
                }
                // if there are more then one bodies, search for the new children and copy the values back to the node
                else if (size > 1) {
                    #pragma omp task final (size <= CUTOFF)
                    {
                        QuadtreeNode* node = new QuadtreeNode(quadrant);

                        if (size > CUTOFF) {
                            node->children = construct_task_with_cutoff(universe, quadrant, indices);
                        }
                        else {
                            node->children = construct(universe, quadrant, indices);
                        }

                        #pragma omp critical
                        subquadrants.push_back(node); // append the node in the subquadrants list
                    }
                }
            }
        }
    }

    return subquadrants;
}

std::vector<BoundingBox> Quadtree::get_bounding_boxes(QuadtreeNode* qtn){
    // traverse quadtree and collect bounding boxes
    std::vector<BoundingBox> result;
    // collect bounding boxes from children
    for(auto child: qtn->children){
        for(auto bb: get_bounding_boxes(child)){
            result.push_back(bb);
        }
    }
    result.push_back(qtn->bounding_box);
    return result;
}








