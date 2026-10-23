#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "sw3/ParticleBlock.hpp"

struct BinaryNode {
    BoundingBox bbox;
    uint32_t left_child = 0xFFFFFFFF;
    uint32_t right_child = 0xFFFFFFFF;
    uint32_t start_index = 0;
    uint32_t count = 0;
};

struct alignas(64) BVH8Node {
    double min_x[8];
    double max_x[8];
    double min_y[8];
    double max_y[8];
    double min_z[8];
    double max_z[8];

    uint32_t child_index[8];
    uint32_t child_count[8];

    void show() const ;
};

inline uint32_t Intersect_Block_BVH8Node(const ParticleBlock& block, const BVH8Node& node, uint32_t* out_lanes) {
    uint32_t count = 0;
    for (int i = 0; i < 8; ++i) {
        if (node.child_index[i] == 0xFFFFFFFF) {
            continue;
        }

        const bool intersect =
            node.max_x[i] > block.box.min_x &&
            node.min_x[i] < block.box.max_x &&
            node.max_y[i] > block.box.min_y &&
            node.min_y[i] < block.box.max_y &&
            node.max_z[i] > block.box.min_z &&
            node.min_z[i] < block.box.max_z;

        if (intersect) {
            out_lanes[count++] = i;
        }
    }
    return count;
}

inline bool Intersect_Block_Block(const ParticleBlock& A, const ParticleBlock& B) {
    return A.box.max_x > B.box.min_x && A.box.min_x < B.box.max_x &&
           A.box.max_y > B.box.min_y && A.box.min_y < B.box.max_y &&
           A.box.max_z > B.box.min_z && A.box.min_z < B.box.max_z;
}

class BVH8 {
public:
    std::vector<BVH8Node> nodes;
    std::vector<ParticleBlock> sorted_;

    BVH8(const std::vector<ParticleBlock>& input_blocks, uint32_t max_leaf_size_);
    void show_bvh() const ;
    void FindCollidingBlocksBVH(const ParticleBlock& query, std::vector<uint32_t>& out_collisions);

private:
    uint32_t max_leaf_size;
    uint32_t build_binary_recursive(std::vector<BinaryNode>& b_nodes, uint32_t start, uint32_t count, const BoundingBox& bbox);
    uint32_t collapse_binary_to_bvh8(const std::vector<BinaryNode>& b_nodes, uint32_t b_root_idx);
};

class BVHCollisionQuery {
public:
    static constexpr uint32_t SafeSentinel = 0xFFFFFFFF;

    class Iterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = uint32_t;
        using difference_type = std::ptrdiff_t;
        using pointer = const uint32_t*;
        using reference = const uint32_t&;

        const ParticleBlock& query;
        const std::vector<BVH8Node>& bvh_nodes;
        const std::vector<ParticleBlock>& leaf;

        uint32_t traversal_stack[128];
        int stack_ptr = 0;

        uint32_t current_node_idx = 0;
        uint32_t hit_lanes[8];
        uint32_t hit_count = 0;
        uint32_t current_hit_lane_idx = 0;

        uint32_t current_leaf_start = 0;
        uint32_t current_leaf_count = 0;
        uint32_t current_leaf_offset = 0;

        bool processing_leaf = false;
        bool processing_node_hits = false;

        uint32_t current_value = SafeSentinel;

        Iterator();
        Iterator(const ParticleBlock& q, const std::vector<BVH8Node>& nodes, const std::vector<ParticleBlock>& leaves);

        void advance();
        reference operator*() const;
        pointer operator->() const;
        Iterator& operator++();
        Iterator operator++(int);
        bool operator==(const Iterator& other) const;
        bool operator!=(const Iterator& other) const;
    };

    BVHCollisionQuery(const ParticleBlock& q, const std::vector<BVH8Node>& nodes, const std::vector<ParticleBlock>& leaves);
    Iterator begin() const;
    Iterator end() const;

private:
    const ParticleBlock& query;
    const std::vector<BVH8Node>& bvh_nodes;
    const std::vector<ParticleBlock>& leaf;
};
