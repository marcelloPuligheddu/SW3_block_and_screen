#include "sw3/ParticleBlock.hpp"

#include <iostream>

namespace {
constexpr double kBlockCutoff = 2.2;
}

void ParticleBlock::show() const {
    std::cout << " Block " << id << std::endl;
    box.show();
    for (int l = 0; l < VEC_WIDTH; ++l) {
        std::cout << " --- " << pos[0][l] << " " << pos[1][l] << " " << pos[2][l] << std::endl;
    }
    std::cout << std::endl;
}

void build_block_box(ParticleBlock& block) {
    for (int i = 0; i < VEC_WIDTH; ++i) {
        BoundingBox atom_box;
        atom_box.min_x = block.pos[0][i] - kBlockCutoff;
        atom_box.max_x = block.pos[0][i] + kBlockCutoff;
        atom_box.min_y = block.pos[1][i] - kBlockCutoff;
        atom_box.max_y = block.pos[1][i] + kBlockCutoff;
        atom_box.min_z = block.pos[2][i] - kBlockCutoff;
        atom_box.max_z = block.pos[2][i] + kBlockCutoff;
        block.box.grow(atom_box);
    }
}
