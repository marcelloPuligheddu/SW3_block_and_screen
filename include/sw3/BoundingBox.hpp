#pragma once

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>

struct BoundingBox {
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    void grow(const BoundingBox& other);
    double surface_area() const;
    void show() const;
};
