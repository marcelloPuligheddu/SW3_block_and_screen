#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <cstddef>
#include <cstdint>
#include <queue>


constexpr int VEC_WIDTH = 1;

struct BoundingBox {
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    double min_z = std::numeric_limits<double>::infinity();
    double max_z = -std::numeric_limits<double>::infinity();

    void grow(const BoundingBox& other) {
        min_x = std::min(min_x, other.min_x);
        max_x = std::max(max_x, other.max_x);
        min_y = std::min(min_y, other.min_y);
        max_y = std::max(max_y, other.max_y);
        min_z = std::min(min_z, other.min_z);
        max_z = std::max(max_z, other.max_z);
    }

    double surface_area() const {
        double dx = max_x - min_x;
        double dy = max_y - min_y;
        double dz = max_z - min_z;
        return 2.0f * (dx * dy + dy * dz + dz * dx);
    }
    void show() const {
            std::cout << " --- X " << min_x << " to " << max_x << std::endl;
            std::cout << " --- Y " << min_y << " to " << max_y << std::endl;
            std::cout << " --- Z " << min_z << " to " << max_z << std::endl;
    }
};

struct ParticleBlock {
    alignas(64) double pos[3][VEC_WIDTH];
    alignas(64) double vel[3][VEC_WIDTH];
    alignas(64) double acc[3][VEC_WIDTH];
    alignas(64) int64_t type[VEC_WIDTH]; 
    bool is_ghost[VEC_WIDTH];
    BoundingBox box;
    uint32_t id;
    uint32_t idx[8];
    void show() const {
        std::cout << " Block " << id << std::endl;
        box.show();
        for ( int l = 0 ; l < VEC_WIDTH; l++ ){
            std::cout << " --- " << pos[0][l] << " " << pos[1][l] << " " << pos[2][l] << std::endl;
        }
        std::cout << std::endl;
    }
};


// Intermediate binary node structure
struct BinaryNode {
    BoundingBox bbox;
    uint32_t left_child = 0xFFFFFFFF;
    uint32_t right_child = 0xFFFFFFFF;
    uint32_t start_index = 0;
    uint32_t count = 0; // If count > 0, it's a leaf node
};

struct alignas(64) BVH8Node {
    // Spatial dimensions for all 8 children lanes
    double min_x[8]; double max_x[8];
    double min_y[8]; double max_y[8];
    double min_z[8]; double max_z[8];

    // Node metadata for all 8 children lanes
    // High bit set on child_index[i] means lane i points to a leaf node.
    uint32_t child_index[8];
    uint32_t child_count[8];
    void show(){
        std::cout << " BVH8 node: " << std::endl;
        for ( int c = 0 ; c < 8 ; c++ ){
            std::cout << child_index[c] << " & " << child_count[c] << std::endl;
            std::cout << " --- " << min_x[c] << " " << max_x[c] << std::endl;
            std::cout << " --- " << min_y[c] << " " << max_y[c] << std::endl;
            std::cout << " --- " << min_z[c] << " " << max_z[c] << std::endl;
            std::cout << " --- --> " << (max_x[c]-min_x[c])*(max_y[c]-min_y[c])*(max_z[c]-min_z[c]) << std::endl;
        }
        
    }
};

// Returns the number of childrens of node that overlap block. 
// Sets out_lanes[ 0 : count-1 ] = idx [0:8] of overlapping child sequentially
inline uint32_t Intersect_Block_BVH8Node(const ParticleBlock& block, const BVH8Node& node, uint32_t* out_lanes) {
    uint32_t count = 0;
    for (int i = 0; i < 8; ++i) {
        if (node.child_index[i] == 0xFFFFFFFF) continue;

        const bool intersect =
            node.max_x[i] > block.box.min_x &&
            node.min_x[i] < block.box.max_x &&
            node.max_y[i] > block.box.min_y &&
            node.min_y[i] < block.box.max_y &&
            node.max_z[i] > block.box.min_z &&
            node.min_z[i] < block.box.max_z;

        if ( intersect ) { out_lanes[count++] = i; }
    }
    return count;
}

// Return true of two boxes overlap
inline bool Intersect_Block_Block(const ParticleBlock& A, const ParticleBlock& B ) {
    return  A.box.max_x > B.box.min_x && A.box.min_x < B.box.max_x &&
            A.box.max_y > B.box.min_y && A.box.min_y < B.box.max_y &&
            A.box.max_z > B.box.min_z && A.box.min_z < B.box.max_z;
}


class BVH8 {
public:
    std::vector<BVH8Node> nodes;
    std::vector<ParticleBlock> sorted_;
    uint32_t max_leaf_size;

    BVH8 ( const std::vector<ParticleBlock>& input_blocks, const uint32_t max_leaf_size_ ) : max_leaf_size { max_leaf_size_ } {
        nodes.clear();

        // ! !! !!! HARD COPY !!! !! ! //
        sorted_ = input_blocks;
        if (sorted_.empty()) return;

        // Step 1: Build a high-quality binary bvh using SAH
        std::vector<BinaryNode> binary_nodes;
        binary_nodes.reserve(sorted_.size() * 2);
        BoundingBox root_box;
        for (const auto& s : sorted_) root_box.grow(s.box);
        build_binary_recursive(binary_nodes, 0, sorted_.size(), root_box);

        // Step 2: Collapse the binary tree into an 8-wide tree layout
        nodes.reserve(binary_nodes.size() / 4);
        collapse_binary_to_bvh8(binary_nodes, 0);
    }

    void show_bvh(){
        for ( auto node : nodes ){
            node.show();
        }
    }

    void FindCollidingBlocksBVH( const ParticleBlock& query, std::vector<uint32_t>& out_collisions) {
        if (nodes.empty()) return;
        // The max stack size is ( MaxDepth * ( BranchingFactor - 1 ) ) + 1
        // 128 means a MaxDepth of (128-1) / (8-1) = 18
        // A bvh8 of depth 18 can accomodate 8**18 leaves, or > 10**17 > UINT32_MAX
        uint32_t traversal_stack[128];
        int stack_ptr = 0;
        traversal_stack[stack_ptr++] = 0;
        uint32_t hit_lanes[8];
        while (stack_ptr > 0) {
            uint32_t current_node_idx = traversal_stack[--stack_ptr];
            const BVH8Node& node = nodes[current_node_idx];
            uint32_t hit_count = Intersect_Block_BVH8Node(query, node, hit_lanes);
            for (uint32_t h = 0; h < hit_count; ++h) {
                uint32_t i = hit_lanes[h];
                uint32_t child_id = node.child_index[i];
                if (child_id & (1 << 31)) {
                    uint32_t element_offset = child_id & ~(1 << 31);
                    uint32_t count = node.child_count[i];
                    for (uint32_t s = 0; s < count; ++s) {
                        const ParticleBlock& target = sorted_[element_offset + s];
                        if ( Intersect_Block_Block ( target, query ) ){
                            out_collisions.push_back(target.id);
                        }
                    }
                } else {
                    if ( stack_ptr < 128) {
                        traversal_stack[stack_ptr++] = child_id;
                    } else {
                       throw std::runtime_error( " Broken Tree, stack limit exceeded, somehow " );
                    }
                }
            }
        }
    }
private:
    // SAH Object-Binning Binary Splitter
    uint32_t build_binary_recursive(std::vector<BinaryNode>& b_nodes, uint32_t start, uint32_t count, const BoundingBox& bbox) {
        uint32_t node_idx = b_nodes.size();
        b_nodes.emplace_back();
        b_nodes[node_idx].bbox = bbox;

        // Base case Leaf node
        if (count <= max_leaf_size) {
            b_nodes[node_idx].start_index = start;
            b_nodes[node_idx].count = count;
            return node_idx;
        }
        // SAH Binning configuration
        const int BINS = 12;
        double best_cost = static_cast<double>(count); // Cost of not splitting
        int best_axis = -1;
        int best_split_bin = -1;

        for (int axis = 0; axis < 3; ++axis) {
            double bounds_min = (axis == 0) ? bbox.min_x : ((axis == 1) ? bbox.min_y : bbox.min_z);
            double bounds_max = (axis == 0) ? bbox.max_x : ((axis == 1) ? bbox.max_y : bbox.max_z);
            double scale = static_cast<double>(BINS) / (bounds_max - bounds_min + 1e-5f);

            struct Bin { BoundingBox box; int count = 0; } bins[BINS];

            for (uint32_t i = 0; i < count; ++i) {
                const auto& s = sorted_[start + i];
                double val = (axis == 0) ? s.box.min_x+s.box.max_x : ((axis == 1) ? s.box.min_y+s.box.max_y : s.box.min_z+s.box.max_z);
                int bin_idx = std::min(BINS - 1, std::max(0, static_cast<int>((val - bounds_min) * scale)));
                bins[bin_idx].count++;
                bins[bin_idx].box.grow(s.box);
            }

            // Sweep over split positions to track lowest Surface Area Heuristic Cost
            for (int i = 0; i < BINS - 1; ++i) {
                BoundingBox box_left, box_right;
                int count_left = 0, count_right = 0;
                for (int j = 0; j <= i; ++j) { box_left.grow(bins[j].box); count_left += bins[j].count; }
                for (int j = i + 1; j < BINS; ++j) { box_right.grow(bins[j].box); count_right += bins[j].count; }

                if (count_left == 0 || count_right == 0) continue;
                double cost = (count_left * box_left.surface_area() + count_right * box_right.surface_area()) / bbox.surface_area();
                if (cost < best_cost) {
                    best_cost = cost;
                    best_axis = axis;
                    best_split_bin = i;
                }
            }
        }
        // Fallback to spatial split if SAH sorting fails to yield optimization gains
        if (best_axis == -1) {
            b_nodes[node_idx].start_index = start;
            b_nodes[node_idx].count = count;
            return node_idx;
        }

        // Re-partition elements based on the winner configurations found above
        double bounds_min = (best_axis == 0) ? bbox.min_x : ((best_axis == 1) ? bbox.min_y : bbox.min_z);
        double bounds_max = (best_axis == 0) ? bbox.max_x : ((best_axis == 1) ? bbox.max_y : bbox.max_z);
        double scale = static_cast<double>(BINS) / (bounds_max - bounds_min + 1e-5f);

        auto split_iterator = std::partition(sorted_.begin() + start, sorted_.begin() + start + count,
            [=](const ParticleBlock& s) {
                double val = (best_axis == 0) ? s.box.min_x+s.box.max_x : ((best_axis == 1) ? s.box.min_y+s.box.max_y : s.box.min_z+s.box.max_z);
                int bin_idx = std::min(BINS - 1, std::max(0, static_cast<int>((val - bounds_min) * scale)));
                return bin_idx <= best_split_bin;
            });

        uint32_t left_count = std::distance(sorted_.begin() + start, split_iterator);
        // Safety fallback partition
        if (left_count == 0 || left_count == count) {
            left_count = count / 2;
        }
        BoundingBox left_box, right_box;
        for (uint32_t i = 0; i < left_count; ++i) left_box.grow(sorted_[start + i].box);
        for (uint32_t i = left_count; i < count; ++i) right_box.grow(sorted_[start + i].box);

        b_nodes[node_idx].left_child = build_binary_recursive(b_nodes, start, left_count, left_box);
        b_nodes[node_idx].right_child = build_binary_recursive(b_nodes, start + left_count, count - left_count, right_box);
        return node_idx;
    }
    // Collapses binary nodes into 8 wide slots based on Surface Area priorities
    uint32_t collapse_binary_to_bvh8(const std::vector<BinaryNode>& b_nodes, uint32_t b_root_idx) {
        uint32_t bvh8_idx = nodes.size();
        nodes.emplace_back();

        // Initialize empty slot states
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

        // Use a priority queue to expand the binary nodes with the largest surface area first
        auto comp = [&](uint32_t idx_a, uint32_t idx_b) {
            return b_nodes[idx_a].bbox.surface_area() < b_nodes[idx_b].bbox.surface_area();
        };
        std::priority_queue<uint32_t, std::vector<uint32_t>, decltype(comp)> cv_queue(comp);

        cv_queue.push(b_root_idx);

        // Keep expanding internal nodes until we gather up to 8 distinct branches
        while (cv_queue.size() < 8) {
            uint32_t top_idx = cv_queue.top();
            // Largest remaining node is a leaf, stop expanding
            if (b_nodes[top_idx].count > 0) break;
            cv_queue.pop();
            cv_queue.push(b_nodes[top_idx].left_child);
            cv_queue.push(b_nodes[top_idx].right_child);
        }
        // Transfer the gathered binary nodes into the 8 layout slots of the BVH8 node
        int slot = 0;
        while (!cv_queue.empty()) {
            uint32_t src_b_idx = cv_queue.top();
            cv_queue.pop();

            const auto& src_node = b_nodes[src_b_idx];
            nodes[bvh8_idx].min_x[slot] = src_node.bbox.min_x; nodes[bvh8_idx].max_x[slot] = src_node.bbox.max_x;
            nodes[bvh8_idx].min_y[slot] = src_node.bbox.min_y; nodes[bvh8_idx].max_y[slot] = src_node.bbox.max_y;
            nodes[bvh8_idx].min_z[slot] = src_node.bbox.min_z; nodes[bvh8_idx].max_z[slot] = src_node.bbox.max_z;

            if (src_node.count > 0) {
                // Point directly to primitive elements list
                nodes[bvh8_idx].child_index[slot] = src_node.start_index | (1 << 31);
                nodes[bvh8_idx].child_count[slot] = src_node.count;
            } else {
                // Inner child layout: dive deeper inside next recursion hierarchy pass step
                uint32_t child8_id = collapse_binary_to_bvh8(b_nodes, src_b_idx);
                nodes[bvh8_idx].child_index[slot] = child8_id;
                nodes[bvh8_idx].child_count[slot] = 0;
            }
            slot++;
        }
        return bvh8_idx;
    }
};
class BVHCollisionQuery {
public:
    // The Sentinel value representing "no more collisions"
    static constexpr uint32_t SafeSentinel = 0xFFFFFFFF;

    struct Iterator {
        using iterator_category = std::input_iterator_tag;
        using value_type        = uint32_t;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const uint32_t*;
        using reference         = const uint32_t&;

        // References to the scene data
        const ParticleBlock& query;
        const std::vector<BVH8Node>& bvh_nodes;
        const std::vector<ParticleBlock>& leaf;

        // Internal State Machine (No heap allocation)
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

        // The current cached result yielded to the user
        uint32_t current_value = SafeSentinel;

        // Constructor for the 'begin' iterator
        Iterator(const ParticleBlock& q, const std::vector<BVH8Node>& nodes, const std::vector<ParticleBlock>& leaves)
            : query(q), bvh_nodes(nodes), leaf(leaves)
        {
            if (!bvh_nodes.empty()) {
                traversal_stack[stack_ptr++] = 0; // Push root
                advance(); // Find the very first collision immediately
            }
        }
        // Constructor for the 'end' (sentinel) iterator
        Iterator() : query(*(ParticleBlock*)nullptr), bvh_nodes(*(std::vector<BVH8Node>*)nullptr), leaf(*(std::vector<ParticleBlock>*)nullptr) {}

        // Core state machine advancement logic
        void advance() {
            while (true) {
                // Step 1: Evaluate primitives inside an active leaf node
                if (processing_leaf) {
                    while (current_leaf_offset < current_leaf_count) {
                        const ParticleBlock& target = leaf[current_leaf_start + current_leaf_offset++];
                        if ( Intersect_Block_Block( target, query ) ){
                            current_value = target.id;
                            return; // Found a collision! Yield control back to user.
                        }
                    }
                    processing_leaf = false;
                }

                // Step 2: Unpack hit lanes from the current BVH8 node
                if (processing_node_hits) {
                    while (current_hit_lane_idx < hit_count) {
                        uint32_t i = hit_lanes[current_hit_lane_idx++];
                        const BVH8Node& node = bvh_nodes[current_node_idx];
                        uint32_t child_id = node.child_index[i];

                        if (child_id & (1 << 31)) { // Leaf node bit flag
                            current_leaf_start = child_id & ~(1 << 31);
                            current_leaf_count = node.child_count[i];
                            current_leaf_offset = 0;
                            processing_leaf = true;
                            break;
                        } else {
                            if (stack_ptr < 128) {
                                traversal_stack[stack_ptr++] = child_id;
                            }
                        }
                    }
                    if (processing_leaf) continue;
                    processing_node_hits = false;
                }
                // Step 3: Pop the next candidate node from the traversal stack
                if (stack_ptr > 0) {
                    current_node_idx = traversal_stack[--stack_ptr];
                    const BVH8Node& node = bvh_nodes[current_node_idx];

                    hit_count = Intersect_Block_BVH8Node(query, node, hit_lanes);
                    current_hit_lane_idx = 0;
                    processing_node_hits = true;
                } else {
                    break; // Tree completely exhausted
                }
            }
            current_value = SafeSentinel; // Signal end of iteration
        }

        // Iterator operator contracts
        reference operator*() const { return current_value; }
        pointer operator->() const { return &current_value; }

        Iterator& operator++() {
            advance();
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = *this;
            advance();
            return tmp;
        }

        // Comparison against the end/sentinel iterator
        bool operator==(const Iterator& other) const {
            return current_value == other.current_value;
        }
        bool operator!=(const Iterator& other) const {
            return !(*this == other);
        }
    };
    // Iterable entry point wrapper
    BVHCollisionQuery(const ParticleBlock& q, const std::vector<BVH8Node>& nodes, const std::vector<ParticleBlock>& leaves)
        : query(q), bvh_nodes(nodes), leaf(leaves) {}

    Iterator begin() const { return Iterator(query, bvh_nodes, leaf); }
    Iterator end() const   { return Iterator(); }

private:
    const ParticleBlock& query;
    const std::vector<BVH8Node>& bvh_nodes;
    const std::vector<ParticleBlock>& leaf;
};

struct ParamSW {
   int64_t n_types;
   const double * as;             
   const double * sg;             
   const double * costheta0;       
   const double * lambda_epsilon;   
   const double * lambda_epsilon_2; 
};

void force_SW_triplet_8x8x8_tile(
      const ParticleBlock * __restrict__ block_i,
      const ParticleBlock * __restrict__ block_j,
      const ParticleBlock * __restrict__ block_k,
      const ParamSW * __restrict__ param,
      double * const __restrict__ force_i_block ) {
         

   alignas(64) double local_accum[3 * VEC_WIDTH] = {0.0};
   const int64_t n = param->n_types;
   const int64_t n2 = n * n;

   // =========================================================================
   // STAGE 1: HOIST & PRE-CACHE ALL IJ & IK TWO-BODY TERMS
   // =========================================================================
   
   // Temporary stack buffers for all 8 interactions with Block J
   alignas(64) double dx_ij[3][VEC_WIDTH][VEC_WIDTH]; // [dim][j][i]
   alignas(64) double r_ij[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double r_ij_inv[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double r_ij_inv_sq[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double ep_ij[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double scaled_rma_ij_sq_rinv[VEC_WIDTH][VEC_WIDTH];
   // Temporary stack buffers for all 8 interactions with Block K
   alignas(64) double dx_ik[3][VEC_WIDTH][VEC_WIDTH]; // [dim][k][i]
   alignas(64) double r_ik[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double r_ik_inv[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double r_ik_inv_sq[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double ep_ik[VEC_WIDTH][VEC_WIDTH];
   alignas(64) double scaled_rma_ik_sq_rinv[VEC_WIDTH][VEC_WIDTH];

   // Evaluates all IJ pairs
   for (int j = 0; j < VEC_WIDTH; j++) {
      const int64_t type_j = block_j->type[j];
      for (int i = 0; i < VEC_WIDTH; i++) {
         double dx0 = block_i->pos[0][i] - block_j->pos[0][j];
         double dx1 = block_i->pos[1][i] - block_j->pos[1][j];
         double dx2 = block_i->pos[2][i] - block_j->pos[2][j];
         dx_ij[0][j][i] = dx0; dx_ij[1][j][i] = dx1; dx_ij[2][j][i] = dx2;
         double r_sq = dx0*dx0 + dx1*dx1 + dx2*dx2 ;
         r_ij[j][i] = std::sqrt(r_sq);
         r_ij_inv[j][i] = 1.0 / r_ij[j][i];
         r_ij_inv_sq[j][i] = r_ij_inv[j][i] * r_ij_inv[j][i];
         double local_as = param->as[block_i->type[i] * n + type_j];
         double local_sg = param->sg[block_i->type[i] * n + type_j];
         double r_minus_as_inv = 1.0 / (r_ij[j][i] - local_as);
         double scaled_inv = local_sg * r_minus_as_inv;
         if ( r_ij[j][i] < local_as ) {
             ep_ij[j][i] = std::exp(scaled_inv);
         } else {
             ep_ij[j][i] = 0.0;
         }
         scaled_rma_ij_sq_rinv[j][i] = scaled_inv * r_minus_as_inv * r_ij_inv[j][i];
      }
   }
   // Removes i == j case
   for ( unsigned int j = 0; j < VEC_WIDTH; j++) {
      for ( unsigned int i = 0; i < VEC_WIDTH; i++) {
         if ( block_i->is_ghost[i] or block_j->is_ghost[j] or (block_i == block_j and i == j) ) {
            r_ij_inv[j][i] = 0.0 ;
            r_ij_inv_sq[j][i] = 0.0;
            ep_ij[j][i] = 0.0;
            scaled_rma_ij_sq_rinv[j][i] = 0.0;
         }
      }
   }


   // Evaluates all IK pairs
   for (int k = 0; k < VEC_WIDTH; k++) {
      const int64_t type_k = block_k->type[k];
      for (int i = 0; i < VEC_WIDTH; i++) {
         double dx0 = block_i->pos[0][i] - block_k->pos[0][k];
         double dx1 = block_i->pos[1][i] - block_k->pos[1][k];
         double dx2 = block_i->pos[2][i] - block_k->pos[2][k];
         dx_ik[0][k][i] = dx0; dx_ik[1][k][i] = dx1; dx_ik[2][k][i] = dx2;
         double r_sq = dx0*dx0 + dx1*dx1 + dx2*dx2 ;
         r_ik[k][i] = std::sqrt(r_sq);
         r_ik_inv[k][i] = 1.0 / r_ik[k][i];
         r_ik_inv_sq[k][i] = r_ik_inv[k][i] * r_ik_inv[k][i];
         double local_as = param->as[block_i->type[i] * n + type_k];
         double local_sg = param->sg[block_i->type[i] * n + type_k];
         double r_minus_as_inv = 1.0 / (r_ik[k][i] - local_as);
         double scaled_inv = local_sg * r_minus_as_inv;
         if ( r_ik[k][i] < local_as ) {
             ep_ik[k][i] = std::exp(scaled_inv);
         } else {
             ep_ik[k][i] = 0.0;
         }
         scaled_rma_ik_sq_rinv[k][i] = scaled_inv * r_minus_as_inv * r_ik_inv[k][i];
      }
   }
   // Removes i == k case
   for ( unsigned int k = 0; k < VEC_WIDTH; k++) {
      for ( unsigned int i = 0; i < VEC_WIDTH; i++) {
         if ( block_i->is_ghost[i] or block_k->is_ghost[k] or (block_i == block_k and i == k) ) {
            r_ik_inv[k][i] = 0.0 ;
            r_ik_inv_sq[k][i] = 0.0;
            ep_ik[k][i] = 0.0;
            scaled_rma_ik_sq_rinv[k][i] = 0.0;
         }
      }
   }

   // =========================================================================
   // STAGE 2: COMPUTES BLOCK ** 3 I J K
   // =========================================================================
   for (int j = 0; j < VEC_WIDTH; j++) {
      const int64_t type_j = block_j->type[j];
      if ( block_j->is_ghost[j] ) continue;
      for (int k = 0; k < VEC_WIDTH; k++) {
         if ( block_k->is_ghost[k] or (block_j == block_k and j == k) ) continue ; // remove j == k cases
         const int64_t offset_ijk_base = type_j * n + block_k->type[k] * n2;

         alignas(64) double fjx_lane[VEC_WIDTH], fjy_lane[VEC_WIDTH], fjz_lane[VEC_WIDTH];
         alignas(64) double fkx_lane[VEC_WIDTH], fky_lane[VEC_WIDTH], fkz_lane[VEC_WIDTH];

         for (int i = 0; i < VEC_WIDTH; i++) {
            double r_ij_ik_inv = r_ij_inv[j][i] * r_ik_inv[k][i];
            
            double dot = dx_ij[0][j][i]*dx_ik[0][k][i] + 
                         dx_ij[1][j][i]*dx_ik[1][k][i] + 
                         dx_ij[2][j][i]*dx_ik[2][k][i];
            // Gathers final pure 3 body coeffiecients
            double cos_theta_0 = param->costheta0[block_i->type[i] + offset_ijk_base];
            double lambda_epsilon = param->lambda_epsilon[block_i->type[i] + offset_ijk_base];
            double lambda_epsilon_2 = param->lambda_epsilon_2[block_i->type[i] + offset_ijk_base];
            // Combines them and 2 body info into 3 body force
            double cos_ijk = dot * r_ij_ik_inv;
            double cos_ijk_minus_cos0 = cos_ijk - cos_theta_0;
            double cos_ijk_minus_cos0_sq = cos_ijk_minus_cos0 * cos_ijk_minus_cos0;
            double exp_ij_ik = ep_ij[j][i] * ep_ik[k][i];
            double scaled_exp_cos = lambda_epsilon * exp_ij_ik * cos_ijk_minus_cos0_sq;
            double exp_term_ij = scaled_exp_cos * scaled_rma_ij_sq_rinv[j][i];
            double exp_term_ik = scaled_exp_cos * scaled_rma_ik_sq_rinv[k][i];
            double scaled_exp_cos_ij_ik = lambda_epsilon_2 * exp_ij_ik * cos_ijk_minus_cos0;
            double ang_term = r_ij_ik_inv * scaled_exp_cos_ij_ik;
            double cos_term_ij_ik = cos_ijk * scaled_exp_cos_ij_ik;
            double cos_term_ij = r_ij_inv_sq[j][i] * cos_term_ij_ik;
            double cos_term_ik = r_ik_inv_sq[k][i] * cos_term_ij_ik;

            fjx_lane[i] = dx_ij[0][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[0][k][i] * ang_term;
            fjy_lane[i] = dx_ij[1][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[1][k][i] * ang_term;
            fjz_lane[i] = dx_ij[2][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[2][k][i] * ang_term;
            fkx_lane[i] = dx_ik[0][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[0][j][i] * ang_term;
            fky_lane[i] = dx_ik[1][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[1][j][i] * ang_term;
            fkz_lane[i] = dx_ik[2][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[2][j][i] * ang_term;
            local_accum[0 * VEC_WIDTH + i] -= (fjx_lane[i] + fkx_lane[i]) ;
            local_accum[1 * VEC_WIDTH + i] -= (fjy_lane[i] + fky_lane[i]) ;
            local_accum[2 * VEC_WIDTH + i] -= (fjz_lane[i] + fkz_lane[i]) ;
         }
      }
   }
   // Writes to out
   for (int idx = 0; idx < 3 * VEC_WIDTH; idx++) { force_i_block[idx] += local_accum[idx]; }
}

struct ParticleDomain {
    int Ncore;
    int Nhalo;
    int Nblocks; // = Ncore + Nhalo
    std::vector<ParticleBlock> blocks; // At least Nblocks
    std::vector<size_t> neighbours;    // At most Ncore * Nblocks, likely C * Ncore 
    std::vector<size_t> offset;        // Ncore + 1
    ParamSW param;
};

void compute_force( ParticleDomain * P  ){
    for (int ib = 0; ib < P -> Ncore; ++ib) {
        const ParticleBlock &bi = P -> blocks[ib];

        for (int jj = P -> offset[ib]; jj < P -> offset[ib + 1]; ++jj) {
            const int jb = P -> neighbours[jj];
            const ParticleBlock &bj = P -> blocks[jb];

            for (int kk = jj ; kk < P -> offset[ib + 1]; ++kk) {
                const int kb = P -> neighbours[kk];
                const ParticleBlock &bk = P -> blocks[kb];

                force_SW_triplet_8x8x8_tile(&bi, &bj, &bk,
                                            &(P->param), &(P->blocks[ib].acc[0][0]) );
            }
        }
    }
}
/*
void find_neigh_from bvh( ParticleDomain * P, BVH * bvh ){

    for (int ib = 0; ib < P -> Ncore; ++ib) {
        P -> offset[ib] = P -> neighbours.size();
        bvh -> query_bvh( P -> bounds[ib], Rlist, candidates);

        for (int jb : candidates) {
            if ( box_distance_squared(bounds[ib],bounds[jb])< Rlist * Rlist ) {
               neighbours.push_back(jb);
            }
        }
    }
}
*/


// TODO, proper max( a * sigma ) in this block / 2
constexpr double cutoff = 2.2;
void build_block_box( ParticleBlock& block ){
    for (int i = 0; i < VEC_WIDTH; ++i) {
        BoundingBox atom_box ;
        atom_box.min_x = block.pos[0][i] - cutoff;
        atom_box.max_x = block.pos[0][i] + cutoff;
        atom_box.min_y = block.pos[1][i] - cutoff;
        atom_box.max_y = block.pos[1][i] + cutoff;
        atom_box.min_z = block.pos[2][i] - cutoff;
        atom_box.max_z = block.pos[2][i] + cutoff;
        block.box.grow( atom_box );
    }
}

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>

// Expand 21 bits so there are two zeroes between consecutive bits.
static inline uint64_t spread3(uint32_t x)
{
    uint64_t v = x & 0x1fffff;
    v = (v | (v << 32)) & 0x1f00000000ffff;
    v = (v | (v << 16)) & 0x1f0000ff0000ff;
    v = (v | (v <<  8)) & 0x100f00f00f00f00f;
    v = (v | (v <<  4)) & 0x10c30c30c30c30c3;
    v = (v | (v <<  2)) & 0x1249249249249249;
    return v;
}

// Map coordinates to a 21-bit integer per dimension and interleave bits.
static inline uint64_t morton3D(double x, double y, double z,
                                double xmin, double ymin, double zmin,
                                double scale)
{
    constexpr double MAX = 2097151.0; // 2^21 - 1

    auto quantize = [&](double v, double lo) -> uint32_t {
        double q = (v - lo) * scale;
        q = std::clamp(q, 0.0, MAX);
        return static_cast<uint32_t>(q);
    };

    const uint32_t ix = quantize(x, xmin);
    const uint32_t iy = quantize(y, ymin);
    const uint32_t iz = quantize(z, zmin);

    return spread3(ix) | (spread3(iy) << 1) | (spread3(iz) << 2);
}


int main() {

    const double N_third = 2*5 ;
    const double lattice = 3.1 ;
    const int n_types = 5;
    const uint32_t N = N_third*N_third*N_third;
    const uint32_t num_blocks = N / VEC_WIDTH ;

    double epsilon = 2.1683;
    double sigma =  2.0951;
    double a =  1.80;
    double lambda = 21.0;
    double gamma = 1.20;
    double cos_theta_0 = -0.333333333333;

    // Build parameter matrix allocations
    std::vector<double> as(n_types * n_types, a * sigma);
    std::vector<double> sg(n_types * n_types, sigma * gamma );
    std::vector<double> costheta0(n_types * n_types * n_types, cos_theta_0);
    std::vector<double> lambda_epsilon(n_types * n_types * n_types, lambda * epsilon);
    std::vector<double> lambda_epsilon_2(n_types * n_types * n_types, lambda * epsilon/2.0);

    ParamSW param{n_types, as.data(), sg.data(), costheta0.data(), lambda_epsilon.data(), lambda_epsilon_2.data()};

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist_pos( -0.1 , 0.1 );
    std::uniform_int_distribution<int> dist_type(0, n_types - 1 - 1);

    // Bounding box of the simulation domain
    double xmin = 0.0;
    double ymin = 0.0;
    double zmin = 0.0;
    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;
    // Unsorted positions
    uint32_t n = 0;
    std::vector<double> pos ( 3 * N );
    for (int i = 0; i < N_third; ++i) {
        for (int j = 0; j < N_third; ++j) {
            for (int k = 0; k < N_third; ++k) {
                pos[0+3*n] = i * lattice + dist_pos(gen);
                pos[1+3*n] = j * lattice + dist_pos(gen);
                pos[2+3*n] = k * lattice + dist_pos(gen);
                xmax = std::max( pos[0+3*n], xmax );
                xmin = std::min( pos[0+3*n], xmin );
                ymax = std::max( pos[1+3*n], ymax );
                ymin = std::min( pos[1+3*n], ymin );
                zmax = std::max( pos[2+3*n], zmax );
                zmin = std::min( pos[2+3*n], zmin );
                n ++;
             }
        }
    }

    // Z curve sort of positions
    const double scale = 2097151.0 / std::max({xmax - xmin, ymax - ymin, zmax - zmin});
    std::vector<uint32_t> order(N);
    std::vector<uint64_t> key(N);
    for (size_t i = 0; i < N; ++i) {
        key[i] = morton3D(pos[3*i], pos[3*i+1], pos[3*i+2], xmin, ymin, zmin, scale);
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return key[a] < key[b]; });
    //std::cout << " Zcurve: " << std::endl;
    //for (size_t i = 0; i < N; ++i) {
    //    std::cout << order[i] <<  " " ;
    //} std::cout << std::endl;

    // Packing atoms in blocks
    std::vector<ParticleBlock> blocks(num_blocks);
    for ( uint32_t i = 0; i < N ; i++ ){
        uint32_t b = i / VEC_WIDTH;
        uint32_t l = i % VEC_WIDTH;
        blocks[b].type[l] = 1 + dist_type(gen);
        blocks[b].pos[0][l] = pos[ 3*order[i]+0 ];
        blocks[b].pos[1][l] = pos[ 3*order[i]+1 ];
        blocks[b].pos[2][l] = pos[ 3*order[i]+2 ];
        blocks[b].is_ghost[l] = false;
        blocks[b].id = b ;
    }
    for ( auto& b : blocks ){ build_block_box( b ); }

    //std::cout << " BLOCKS: " << std::endl;
    //for ( const auto &b : blocks ){ b.show(); }

    std::vector<double> global_force(3 * VEC_WIDTH, 0.0);

    std::cout << "=========================================================\n";
    std::cout << "Starting Dense Triplet O(N^3) Stillinger-Weber Benchmark\n";
    std::cout << "Evaluating " << num_blocks << "^3 = " << num_blocks*num_blocks*num_blocks << " Block-Triplets\n";
    std::cout << "=========================================================\n\n";

    auto brute_ni = std::min( uint32_t(1+12*12*12/VEC_WIDTH),num_blocks);
    auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < brute_ni; ++i) {
        for (auto j = 0; j < brute_ni; ++j) {
            for (auto k = 0; k < brute_ni; ++k) {
                force_SW_triplet_8x8x8_tile(&blocks[i], &blocks[j], &blocks[k], &param, global_force.data());
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    double sw_ms = std::chrono::duration<double, std::milli>(end - start).count();
    size_t NiNN = size_t(VEC_WIDTH)*VEC_WIDTH*VEC_WIDTH*brute_ni*brute_ni*brute_ni;
    std::cout << " " << NiNN / 1.e6 / (sw_ms / 1.e3) << " MSW / s " << global_force[0] << std::endl;
    if ( brute_ni != num_blocks ) { std::cout << " Force unreliable, brute force calc. was cut for time " << std::endl;}
    std::cout << "=========================================================\n";

    BVH8 bvh( blocks, 4 );
    std::vector<double> global_force2(3 * VEC_WIDTH, 0.0);
    auto start2 = std::chrono::high_resolution_clock::now();
    uint64_t n3t=0, n2t=0 ;
    for ( const auto& i_block : blocks ){
        auto nn_j = BVHCollisionQuery( i_block, bvh.nodes, bvh.sorted_ );
        auto nn_k = BVHCollisionQuery( i_block, bvh.nodes, bvh.sorted_ );
        //std::cout << " ----------------------------------------------------- " << std::endl;
        //i_block.show();
        for ( auto j_block_id : nn_j ){
            n2t ++ ;
            //blocks[j_block_id].show();
            for ( auto k_block_id : nn_k ){
                force_SW_triplet_8x8x8_tile(&i_block, &blocks[j_block_id], &blocks[k_block_id], &param, global_force2.data());
                n3t ++ ;
            }
        }
        //std::cout << " ----------------------------------------------------- " << std::endl;
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    double sw2_ms = std::chrono::duration<double, std::milli>(end2 - start2).count();
    std::cout << " " << N / (sw2_ms / 1.e3) << " SW / s " << global_force2[0] << std::endl;
    std::cout << " " << double(VEC_WIDTH)*VEC_WIDTH*n2t / (VEC_WIDTH*num_blocks) << " doublets / atom " << std::endl;
    std::cout << " " << double(VEC_WIDTH)*VEC_WIDTH*VEC_WIDTH*n3t / (VEC_WIDTH*num_blocks) << " triplets / atom " << std::endl;
    std::cout << " " << 1.e-6 * VEC_WIDTH*VEC_WIDTH*VEC_WIDTH*n3t / (sw2_ms / 1.e3) << " MTs / s " << std::endl;
    std::cout << "=========================================================\n";

}
