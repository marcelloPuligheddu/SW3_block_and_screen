#pragma once

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "sw3/BoundingBox.hpp"

constexpr int VEC_WIDTH = 1;

struct ParticleBlock {
    alignas(64) double pos[3][VEC_WIDTH];
    alignas(64) double vel[3][VEC_WIDTH];
    alignas(64) double acc[3][VEC_WIDTH];
    alignas(64) int64_t type[VEC_WIDTH];
    bool is_ghost[VEC_WIDTH];
    BoundingBox box;
    uint32_t id;
    uint32_t idx[8];

    void show() const;
};

struct ParamSW {
    int64_t n_types;
    const double* as;
    const double* sg;
    const double* costheta0;
    const double* lambda_epsilon;
    const double* lambda_epsilon_2;
};

struct ParticleDomain {
    int Ncore;
    int Nhalo;
    int Nblocks;
    std::vector<ParticleBlock> blocks;
    std::vector<size_t> neighbours;
    std::vector<size_t> offset;
    ParamSW param;
};

void build_block_box(ParticleBlock& block);
