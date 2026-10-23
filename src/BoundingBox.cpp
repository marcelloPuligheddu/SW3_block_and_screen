#include "sw3/BoundingBox.hpp"

#include <algorithm>
#include <iostream>

void BoundingBox::grow(const BoundingBox& other) {
    min_x = std::min(min_x, other.min_x);
    max_x = std::max(max_x, other.max_x);
    min_y = std::min(min_y, other.min_y);
    max_y = std::max(max_y, other.max_y);
    min_z = std::min(min_z, other.min_z);
    max_z = std::max(max_z, other.max_z);
}

double BoundingBox::surface_area() const {
    const double dx = max_x - min_x;
    const double dy = max_y - min_y;
    const double dz = max_z - min_z;
    return 2.0 * (dx * dy + dy * dz + dz * dx);
}

void BoundingBox::show() const {
    std::cout << " --- X " << min_x << " to " << max_x << std::endl;
    std::cout << " --- Y " << min_y << " to " << max_y << std::endl;
    std::cout << " --- Z " << min_z << " to " << max_z << std::endl;
}
