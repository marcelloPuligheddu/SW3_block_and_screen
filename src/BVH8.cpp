#include "sw3/BVH8.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <stdexcept>

void BVH8Node::show() {
    std::cout << " BVH8 node: " << std::endl;
    for (int c = 0; c < 8; ++c) {
        std::cout << child_index[c] << " & " << child_count[c] << std::endl;
        std::cout << " --- " << min_x[c] << " " << max_x[c] << std::endl;
        std::cout << " --- " << min_y[c] << " " << max_y[c] << std::endl;
        std::cout << " --- " << min_z[c] << " " << max_z[c] << std::endl;
        std::cout << " --- --> " << (max_x[c] - min_x[c]) * (max_y[c] - min_y[c]) * (max_z[c] - min_z[c]) << std::endl;
    }
}

BVH8::BVH8(const std::vector<ParticleBlock>& input_blocks, uint32_t max_leaf_size_)
    : max_leaf_size(max_leaf_size_) {
    nodes.clear();
    sorted_ = input_blocks;
    if (sorted_.empty()) {
        return;
    }

    std::vector<BinaryNode> binary_nodes;
    binary_nodes.reserve(sorted_.size() * 2);
    BoundingBox root_box;
    for (const auto& s : sorted_) {
        root_box.grow(s.box);
    }
    build_binary_recursive(binary_nodes, 0, static_cast<uint32_t>(sorted_.size()), root_box);
    nodes.reserve(binary_nodes.size() / 4);
    collapse_binary_to_bvh8(binary_nodes, 0);
}

void BVH8::show_bvh() {
    for (const auto& node : nodes) {
        node.show();
    }
}

void BVH8::FindCollidingBlocksBVH(const ParticleBlock& query, std::vector<uint32_t>& out_collisions) {
    if (nodes.empty()) {
        return;
    }
    uint32_t traversal_stack[128];
    int stack_ptr = 0;
    traversal_stack[stack_ptr++] = 0;
    uint32_t hit_lanes[8];

    while (stack_ptr > 0) {
        const uint32_t current_node_idx = traversal_stack[--stack_ptr];
        const BVH8Node& node = nodes[current_node_idx];
        const uint32_t hit_count = Intersect_Block_BVH8Node(query, node, hit_lanes);
        for (uint32_t h = 0; h < hit_count; ++h) {
            const uint32_t i = hit_lanes[h];
            const uint32_t child_id = node.child_index[i];
            if (child_id & (1u << 31)) {
                const uint32_t element_offset = child_id & ~(1u << 31);
                const uint32_t count = node.child_count[i];
                for (uint32_t s = 0; s < count; ++s) {
                    const ParticleBlock& target = sorted_[element_offset + s];
                    if (Intersect_Block_Block(target, query)) {
                        out_collisions.push_back(target.id);
                    }
                }
            } else {
                if (stack_ptr < 128) {
                    traversal_stack[stack_ptr++] = child_id;
                } else {
                    throw std::runtime_error("Broken Tree, stack limit exceeded, somehow");
                }
            }
        }
    }
}

uint32_t BVH8::build_binary_recursive(std::vector<BinaryNode>& b_nodes,
                                     uint32_t start,
                                     uint32_t count,
                                     const BoundingBox& bbox) {
    const uint32_t node_idx = static_cast<uint32_t>(b_nodes.size());
    b_nodes.emplace_back();
    b_nodes[node_idx].bbox = bbox;

    if (count <= max_leaf_size) {
        b_nodes[node_idx].start_index = start;
        b_nodes[node_idx].count = count;
        return node_idx;
    }

    constexpr int BINS = 12;
    double best_cost = static_cast<double>(count);
    int best_axis = -1;
    int best_split_bin = -1;

    for (int axis = 0; axis < 3; ++axis) {
        const double bounds_min = (axis == 0) ? bbox.min_x : ((axis == 1) ? bbox.min_y : bbox.min_z);
        const double bounds_max = (axis == 0) ? bbox.max_x : ((axis == 1) ? bbox.max_y : bbox.max_z);
        const double scale = static_cast<double>(BINS) / (bounds_max - bounds_min + 1e-5);

        struct Bin {
            BoundingBox box;
            int count = 0;
        } bins[BINS];

        for (uint32_t i = 0; i < count; ++i) {
            const auto& s = sorted_[start + i];
            const double val = (axis == 0) ? s.box.min_x + s.box.max_x : ((axis == 1) ? s.box.min_y + s.box.max_y : s.box.min_z + s.box.max_z);
            const int bin_idx = std::min(BINS - 1, std::max(0, static_cast<int>((val - bounds_min) * scale)));
            bins[bin_idx].count++;
            bins[bin_idx].box.grow(s.box);
        }

        for (int i = 0; i < BINS - 1; ++i) {
            BoundingBox box_left, box_right;
            int count_left = 0, count_right = 0;
            for (int j = 0; j <= i; ++j) {
                box_left.grow(bins[j].box);
                count_left += bins[j].count;
            }
            for (int j = i + 1; j < BINS; ++j) {
                box_right.grow(bins[j].box);
                count_right += bins[j].count;
            }

            if (count_left == 0 || count_right == 0) {
                continue;
            }

            const double cost = (count_left * box_left.surface_area() + count_right * box_right.surface_area()) / bbox.surface_area();
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split_bin = i;
            }
        }
    }

    if (best_axis == -1) {
        b_nodes[node_idx].start_index = start;
        b_nodes[node_idx].count = count;
        return node_idx;
    }

    const double bounds_min = (best_axis == 0) ? bbox.min_x : ((best_axis == 1) ? bbox.min_y : bbox.min_z);
    const double bounds_max = (best_axis == 0) ? bbox.max_x : ((best_axis == 1) ? bbox.max_y : bbox.max_z);
    const double scale = static_cast<double>(BINS) / (bounds_max - bounds_min + 1e-5);

    auto split_iterator = std::partition(sorted_.begin() + start,
                                        sorted_.begin() + start + static_cast<std::ptrdiff_t>(count),
                                        [=](const ParticleBlock& s) {
                                            const double val = (best_axis == 0) ? s.box.min_x + s.box.max_x
                                                : ((best_axis == 1) ? s.box.min_y + s.box.max_y : s.box.min_z + s.box.max_z);
                                            const int bin_idx = std::min(BINS - 1, std::max(0, static_cast<int>((val - bounds_min) * scale)));
                                            return bin_idx <= best_split_bin;
                                        });

    uint32_t left_count = static_cast<uint32_t>(std::distance(sorted_.begin() + start, split_iterator));
    if (left_count == 0 || left_count == count) {
        left_count = count / 2;
    }

    BoundingBox left_box, right_box;
    for (uint32_t i = 0; i < left_count; ++i) {
        left_box.grow(sorted_[start + i].box);
    }
    for (uint32_t i = left_count; i < count; ++i) {
        right_box.grow(sorted_[start + i].box);
    }

    b_nodes[node_idx].left_child = build_binary_recursive(b_nodes, start, left_count, left_box);
    b_nodes[node_idx].right_child = build_binary_recursive(b_nodes, start + left_count, count - left_count, right_box);
    return node_idx;
}

uint32_t BVH8::collapse_binary_to_bvh8(const std::vector<BinaryNode>& b_nodes, uint32_t b_root_idx) {
    const uint32_t bvh8_idx = static_cast<uint32_t>(nodes.size());
    nodes.emplace_back();

    for (int i = 0; i < 8; ++i) {
        nodes[bvh8_idx].child_index[i] = 0xFFFFFFFF;
        nodes[bvh8_idx].child_count[i] = 0;
        nodes[bvh8_idx].min_x[i] = std::numeric_limits<double>::infinity();
        nodes[bvh8_idx].max_x[i] = -std::numeric_limits<double>::infinity();
        nodes[bvh8_idx].min_y[i] = std::numeric_limits<double>::infinity();
        nodes[bvh8_idx].max_y[i] = -std::numeric_limits<double>::infinity();
        nodes[bvh8_idx].min_z[i] = std::numeric_limits<double>::infinity();
        nodes[bvh8_idx].max_z[i] = -std::numeric_limits<double>::infinity();
    }

    auto comp = [&](uint32_t idx_a, uint32_t idx_b) {
        return b_nodes[idx_a].bbox.surface_area() < b_nodes[idx_b].bbox.surface_area();
    };
    std::priority_queue<uint32_t, std::vector<uint32_t>, decltype(comp)> cv_queue(comp);
    cv_queue.push(b_root_idx);

    while (cv_queue.size() < 8) {
        const uint32_t top_idx = cv_queue.top();
        if (b_nodes[top_idx].count > 0) {
            break;
        }
        cv_queue.pop();
        cv_queue.push(b_nodes[top_idx].left_child);
        cv_queue.push(b_nodes[top_idx].right_child);
    }

    int slot = 0;
    while (!cv_queue.empty()) {
        const uint32_t src_b_idx = cv_queue.top();
        cv_queue.pop();

        const auto& src_node = b_nodes[src_b_idx];
        nodes[bvh8_idx].min_x[slot] = src_node.bbox.min_x;
        nodes[bvh8_idx].max_x[slot] = src_node.bbox.max_x;
        nodes[bvh8_idx].min_y[slot] = src_node.bbox.min_y;
        nodes[bvh8_idx].max_y[slot] = src_node.bbox.max_y;
        nodes[bvh8_idx].min_z[slot] = src_node.bbox.min_z;
        nodes[bvh8_idx].max_z[slot] = src_node.bbox.max_z;

        if (src_node.count > 0) {
            nodes[bvh8_idx].child_index[slot] = src_node.start_index | (1u << 31);
            nodes[bvh8_idx].child_count[slot] = src_node.count;
        } else {
            const uint32_t child8_id = collapse_binary_to_bvh8(b_nodes, src_b_idx);
            nodes[bvh8_idx].child_index[slot] = child8_id;
            nodes[bvh8_idx].child_count[slot] = 0;
        }
        ++slot;
    }
    return bvh8_idx;
}

BVHCollisionQuery::Iterator::Iterator() : query(*(ParticleBlock*)nullptr), bvh_nodes(*(std::vector<BVH8Node>*)nullptr), leaf(*(std::vector<ParticleBlock>*)nullptr) {}

BVHCollisionQuery::Iterator::Iterator(const ParticleBlock& q,
                                     const std::vector<BVH8Node>& nodes,
                                     const std::vector<ParticleBlock>& leaves)
    : query(q), bvh_nodes(nodes), leaf(leaves) {
    if (!bvh_nodes.empty()) {
        traversal_stack[stack_ptr++] = 0;
        advance();
    }
}

void BVHCollisionQuery::Iterator::advance() {
    while (true) {
        if (processing_leaf) {
            while (current_leaf_offset < current_leaf_count) {
                const ParticleBlock& target = leaf[current_leaf_start + current_leaf_offset++];
                if (Intersect_Block_Block(target, query)) {
                    current_value = target.id;
                    return;
                }
            }
            processing_leaf = false;
        }

        if (processing_node_hits) {
            while (current_hit_lane_idx < hit_count) {
                const uint32_t i = hit_lanes[current_hit_lane_idx++];
                const BVH8Node& node = bvh_nodes[current_node_idx];
                const uint32_t child_id = node.child_index[i];

                if (child_id & (1u << 31)) {
                    current_leaf_start = child_id & ~(1u << 31);
                    current_leaf_count = node.child_count[i];
                    current_leaf_offset = 0;
                    processing_leaf = true;
                    break;
                } else if (stack_ptr < 128) {
                    traversal_stack[stack_ptr++] = child_id;
                }
            }
            if (processing_leaf) {
                continue;
            }
            processing_node_hits = false;
        }

        if (stack_ptr > 0) {
            current_node_idx = traversal_stack[--stack_ptr];
            const BVH8Node& node = bvh_nodes[current_node_idx];
            hit_count = Intersect_Block_BVH8Node(query, node, hit_lanes);
            current_hit_lane_idx = 0;
            processing_node_hits = true;
        } else {
            break;
        }
    }
    current_value = SafeSentinel;
}

BVHCollisionQuery::Iterator::reference BVHCollisionQuery::Iterator::operator*() const {
    return current_value;
}

BVHCollisionQuery::Iterator::pointer BVHCollisionQuery::Iterator::operator->() const {
    return &current_value;
}

BVHCollisionQuery::Iterator& BVHCollisionQuery::Iterator::operator++() {
    advance();
    return *this;
}

BVHCollisionQuery::Iterator BVHCollisionQuery::Iterator::operator++(int) {
    Iterator tmp = *this;
    advance();
    return tmp;
}

bool BVHCollisionQuery::Iterator::operator==(const Iterator& other) const {
    return current_value == other.current_value;
}

bool BVHCollisionQuery::Iterator::operator!=(const Iterator& other) const {
    return !(*this == other);
}

BVHCollisionQuery::BVHCollisionQuery(const ParticleBlock& q,
                                     const std::vector<BVH8Node>& nodes,
                                     const std::vector<ParticleBlock>& leaves)
    : query(q), bvh_nodes(nodes), leaf(leaves) {}

BVHCollisionQuery::Iterator BVHCollisionQuery::begin() const {
    return Iterator(query, bvh_nodes, leaf);
}

BVHCollisionQuery::Iterator BVHCollisionQuery::end() const {
    return Iterator();
}
