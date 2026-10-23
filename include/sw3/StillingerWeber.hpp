#pragma once

#include "sw3/BVH8.hpp"
#include "sw3/ParticleBlock.hpp"

void force_SW_triplet_8x8x8_tile(
    const ParticleBlock* __restrict__ block_i,
    const ParticleBlock* __restrict__ block_j,
    const ParticleBlock* __restrict__ block_k,
    const ParamSW* __restrict__ param,
    double* const __restrict__ force_i_block);

void compute_force(ParticleDomain* P);
