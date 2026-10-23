#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "sw3/BVH8.hpp"
#include "sw3/StillingerWeber.hpp"

namespace {

// Expand 21 bits so there are two zeroes between consecutive bits.
uint64_t spread3(uint32_t x) {
    uint64_t v = x & 0x1fffff;
    v = (v | (v << 32)) & 0x1f00000000ffff;
    v = (v | (v << 16)) & 0x1f0000ff0000ff;
    v = (v | (v << 8)) & 0x100f00f00f00f00f;
    v = (v | (v << 4)) & 0x10c30c30c30c30c3;
    v = (v | (v << 2)) & 0x1249249249249249;
    return v;
}

uint64_t morton3D(double x, double y, double z, double xmin, double ymin, double zmin, double scale) {
    constexpr double MAX = 2097151.0;
    const auto quantize = [&](double v, double lo) -> uint32_t {
        double q = (v - lo) * scale;
        q = std::clamp(q, 0.0, MAX);
        return static_cast<uint32_t>(q);
    };

    const uint32_t ix = quantize(x, xmin);
    const uint32_t iy = quantize(y, ymin);
    const uint32_t iz = quantize(z, zmin);
    return spread3(ix) | (spread3(iy) << 1) | (spread3(iz) << 2);
}

}  // namespace

int main() {
    const double N_third = 2.0 * 5.0;
    const double lattice = 3.1;
    const int n_types = 5;
    const uint32_t N = static_cast<uint32_t>(N_third * N_third * N_third);
    const uint32_t num_blocks = N / VEC_WIDTH;

    double epsilon = 2.1683;
    double sigma = 2.0951;
    double a = 1.80;
    double lambda = 21.0;
    double gamma = 1.20;
    double cos_theta_0 = -0.333333333333;

    std::vector<double> as(n_types * n_types, a * sigma);
    std::vector<double> sg(n_types * n_types, sigma * gamma);
    std::vector<double> costheta0(n_types * n_types * n_types, cos_theta_0);
    std::vector<double> lambda_epsilon(n_types * n_types * n_types, lambda * epsilon);
    std::vector<double> lambda_epsilon_2(n_types * n_types * n_types, lambda * epsilon / 2.0);

    ParamSW param{n_types, as.data(), sg.data(), costheta0.data(), lambda_epsilon.data(), lambda_epsilon_2.data()};

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist_pos(-0.1, 0.1);
    std::uniform_int_distribution<int> dist_type(0, n_types - 1 - 1);

    double xmin = 0.0;
    double ymin = 0.0;
    double zmin = 0.0;
    double xmax = 0.0;
    double ymax = 0.0;
    double zmax = 0.0;

    uint32_t n = 0;
    std::vector<double> pos(3 * N);
    for (int i = 0; i < N_third; ++i) {
        for (int j = 0; j < N_third; ++j) {
            for (int k = 0; k < N_third; ++k) {
                pos[0 + 3 * n] = i * lattice + dist_pos(gen);
                pos[1 + 3 * n] = j * lattice + dist_pos(gen);
                pos[2 + 3 * n] = k * lattice + dist_pos(gen);
                xmax = std::max(pos[0 + 3 * n], xmax);
                xmin = std::min(pos[0 + 3 * n], xmin);
                ymax = std::max(pos[1 + 3 * n], ymax);
                ymin = std::min(pos[1 + 3 * n], ymin);
                zmax = std::max(pos[2 + 3 * n], zmax);
                zmin = std::min(pos[2 + 3 * n], zmin);
                ++n;
            }
        }
    }

    const double scale = 2097151.0 / std::max({xmax - xmin, ymax - ymin, zmax - zmin});
    std::vector<uint32_t> order(N);
    std::vector<uint64_t> key(N);
    for (size_t i = 0; i < N; ++i) {
        key[i] = morton3D(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2], xmin, ymin, zmin, scale);
        order[i] = static_cast<uint32_t>(i);
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return key[a] < key[b]; });

    std::vector<ParticleBlock> blocks(num_blocks);
    for (uint32_t i = 0; i < N; ++i) {
        const uint32_t b = i / VEC_WIDTH;
        const uint32_t l = i % VEC_WIDTH;
        blocks[b].type[l] = 1 + dist_type(gen);
        blocks[b].pos[0][l] = pos[3 * order[i] + 0];
        blocks[b].pos[1][l] = pos[3 * order[i] + 1];
        blocks[b].pos[2][l] = pos[3 * order[i] + 2];
        blocks[b].is_ghost[l] = false;
        blocks[b].id = b;
    }
    for (auto& b : blocks) {
        build_block_box(b);
    }

    std::vector<double> global_force(3 * VEC_WIDTH, 0.0);

    std::cout << "=========================================================" << std::endl;
    std::cout << "Starting Dense Triplet O(N^3) Stillinger-Weber Benchmark" << std::endl;
    std::cout << "Evaluating " << num_blocks << "^3 = " << num_blocks * num_blocks * num_blocks << " Block-Triplets" << std::endl;
    std::cout << "=========================================================" << std::endl << std::endl;

    const auto brute_ni = std::min<uint32_t>(1 + 12 * 12 * 12 / VEC_WIDTH, num_blocks);
    const auto start = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < brute_ni; ++i) {
        for (auto j = 0; j < brute_ni; ++j) {
            for (auto k = 0; k < brute_ni; ++k) {
                force_SW_triplet_8x8x8_tile(&blocks[i], &blocks[j], &blocks[k], &param, global_force.data());
            }
        }
    }
    const auto end = std::chrono::high_resolution_clock::now();
    const double sw_ms = std::chrono::duration<double, std::milli>(end - start).count();
    const size_t NiNN = size_t(VEC_WIDTH) * VEC_WIDTH * VEC_WIDTH * brute_ni * brute_ni * brute_ni;
    std::cout << " " << NiNN / 1.e6 / (sw_ms / 1.e3) << " MSW / s " << global_force[0] << std::endl;
    if (brute_ni != num_blocks) {
        std::cout << " Force unreliable, brute force calc. was cut for time " << std::endl;
    }
    std::cout << "=========================================================" << std::endl;

    BVH8 bvh(blocks, 4);
    std::vector<double> global_force2(3 * VEC_WIDTH, 0.0);
    const auto start2 = std::chrono::high_resolution_clock::now();
    uint64_t n3t = 0, n2t = 0;
    for (const auto& i_block : blocks) {
        auto nn_j = BVHCollisionQuery(i_block, bvh.nodes, bvh.sorted_);
        auto nn_k = BVHCollisionQuery(i_block, bvh.nodes, bvh.sorted_);
        for (auto j_block_id : nn_j) {
            ++n2t;
            for (auto k_block_id : nn_k) {
                force_SW_triplet_8x8x8_tile(&i_block, &blocks[j_block_id], &blocks[k_block_id], &param, global_force2.data());
                ++n3t;
            }
        }
    }
    const auto end2 = std::chrono::high_resolution_clock::now();
    const double sw2_ms = std::chrono::duration<double, std::milli>(end2 - start2).count();
    std::cout << " " << N / (sw2_ms / 1.e3) << " SW / s " << global_force2[0] << std::endl;
    std::cout << " " << double(VEC_WIDTH) * VEC_WIDTH * n2t / (VEC_WIDTH * num_blocks) << " doublets / atom " << std::endl;
    std::cout << " " << double(VEC_WIDTH) * VEC_WIDTH * VEC_WIDTH * n3t / (VEC_WIDTH * num_blocks) << " triplets / atom " << std::endl;
    std::cout << " " << 1.e-6 * VEC_WIDTH * VEC_WIDTH * VEC_WIDTH * n3t / (sw2_ms / 1.e3) << " MTs / s " << std::endl;
    std::cout << "=========================================================" << std::endl;
    return 0;
}
