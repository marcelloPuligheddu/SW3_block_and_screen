#include "sw3/StillingerWeber.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

void force_SW_triplet_8x8x8_tile(
    const ParticleBlock* __restrict__ block_i,
    const ParticleBlock* __restrict__ block_j,
    const ParticleBlock* __restrict__ block_k,
    const ParamSW* __restrict__ param,
    double* const __restrict__ force_i_block) {

    alignas(64) double local_accum[3 * VEC_WIDTH] = {0.0};
    const int64_t n = param->n_types;
    const int64_t n2 = n * n;

    alignas(64) double dx_ij[3][VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ij[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ij_inv[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ij_inv_sq[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double ep_ij[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double scaled_rma_ij_sq_rinv[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double dx_ik[3][VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ik[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ik_inv[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double r_ik_inv_sq[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double ep_ik[VEC_WIDTH][VEC_WIDTH];
    alignas(64) double scaled_rma_ik_sq_rinv[VEC_WIDTH][VEC_WIDTH];

    for (int j = 0; j < VEC_WIDTH; ++j) {
        const int64_t type_j = block_j->type[j];
        for (int i = 0; i < VEC_WIDTH; ++i) {
            double dx0 = block_i->pos[0][i] - block_j->pos[0][j];
            double dx1 = block_i->pos[1][i] - block_j->pos[1][j];
            double dx2 = block_i->pos[2][i] - block_j->pos[2][j];
            dx_ij[0][j][i] = dx0;
            dx_ij[1][j][i] = dx1;
            dx_ij[2][j][i] = dx2;
            const double r_sq = dx0 * dx0 + dx1 * dx1 + dx2 * dx2;
            r_ij[j][i] = std::sqrt(r_sq);
            r_ij_inv[j][i] = 1.0 / r_ij[j][i];
            r_ij_inv_sq[j][i] = r_ij_inv[j][i] * r_ij_inv[j][i];
            const double local_as = param->as[block_i->type[i] * n + type_j];
            const double local_sg = param->sg[block_i->type[i] * n + type_j];
            const double r_minus_as_inv = 1.0 / (r_ij[j][i] - local_as);
            const double scaled_inv = local_sg * r_minus_as_inv;
            if (r_ij[j][i] < local_as) {
                ep_ij[j][i] = std::exp(scaled_inv);
            } else {
                ep_ij[j][i] = 0.0;
            }
            scaled_rma_ij_sq_rinv[j][i] = scaled_inv * r_minus_as_inv * r_ij_inv[j][i];
        }
    }

    for (unsigned int j = 0; j < VEC_WIDTH; ++j) {
        for (unsigned int i = 0; i < VEC_WIDTH; ++i) {
            if (block_i->is_ghost[i] || block_j->is_ghost[j] || (block_i == block_j && i == j)) {
                r_ij_inv[j][i] = 0.0;
                r_ij_inv_sq[j][i] = 0.0;
                ep_ij[j][i] = 0.0;
                scaled_rma_ij_sq_rinv[j][i] = 0.0;
            }
        }
    }

    for (int k = 0; k < VEC_WIDTH; ++k) {
        const int64_t type_k = block_k->type[k];
        for (int i = 0; i < VEC_WIDTH; ++i) {
            double dx0 = block_i->pos[0][i] - block_k->pos[0][k];
            double dx1 = block_i->pos[1][i] - block_k->pos[1][k];
            double dx2 = block_i->pos[2][i] - block_k->pos[2][k];
            dx_ik[0][k][i] = dx0;
            dx_ik[1][k][i] = dx1;
            dx_ik[2][k][i] = dx2;
            const double r_sq = dx0 * dx0 + dx1 * dx1 + dx2 * dx2;
            r_ik[k][i] = std::sqrt(r_sq);
            r_ik_inv[k][i] = 1.0 / r_ik[k][i];
            r_ik_inv_sq[k][i] = r_ik_inv[k][i] * r_ik_inv[k][i];
            const double local_as = param->as[block_i->type[i] * n + type_k];
            const double local_sg = param->sg[block_i->type[i] * n + type_k];
            const double r_minus_as_inv = 1.0 / (r_ik[k][i] - local_as);
            const double scaled_inv = local_sg * r_minus_as_inv;
            if (r_ik[k][i] < local_as) {
                ep_ik[k][i] = std::exp(scaled_inv);
            } else {
                ep_ik[k][i] = 0.0;
            }
            scaled_rma_ik_sq_rinv[k][i] = scaled_inv * r_minus_as_inv * r_ik_inv[k][i];
        }
    }

    for (unsigned int k = 0; k < VEC_WIDTH; ++k) {
        for (unsigned int i = 0; i < VEC_WIDTH; ++i) {
            if (block_i->is_ghost[i] || block_k->is_ghost[k] || (block_i == block_k && i == k)) {
                r_ik_inv[k][i] = 0.0;
                r_ik_inv_sq[k][i] = 0.0;
                ep_ik[k][i] = 0.0;
                scaled_rma_ik_sq_rinv[k][i] = 0.0;
            }
        }
    }

    for (int j = 0; j < VEC_WIDTH; ++j) {
        const int64_t type_j = block_j->type[j];
        if (block_j->is_ghost[j]) {
            continue;
        }
        for (int k = 0; k < VEC_WIDTH; ++k) {
            if (block_k->is_ghost[k] || (block_j == block_k && j == k)) {
                continue;
            }
            const int64_t offset_ijk_base = type_j * n + block_k->type[k] * n2;

            alignas(64) double fjx_lane[VEC_WIDTH], fjy_lane[VEC_WIDTH], fjz_lane[VEC_WIDTH];
            alignas(64) double fkx_lane[VEC_WIDTH], fky_lane[VEC_WIDTH], fkz_lane[VEC_WIDTH];

            for (int i = 0; i < VEC_WIDTH; ++i) {
                const double r_ij_ik_inv = r_ij_inv[j][i] * r_ik_inv[k][i];
                const double dot = dx_ij[0][j][i] * dx_ik[0][k][i] +
                                   dx_ij[1][j][i] * dx_ik[1][k][i] +
                                   dx_ij[2][j][i] * dx_ik[2][k][i];
                const double cos_theta_0 = param->costheta0[block_i->type[i] + offset_ijk_base];
                const double lambda_epsilon = param->lambda_epsilon[block_i->type[i] + offset_ijk_base];
                const double lambda_epsilon_2 = param->lambda_epsilon_2[block_i->type[i] + offset_ijk_base];
                const double cos_ijk = dot * r_ij_ik_inv;
                const double cos_ijk_minus_cos0 = cos_ijk - cos_theta_0;
                const double cos_ijk_minus_cos0_sq = cos_ijk_minus_cos0 * cos_ijk_minus_cos0;
                const double exp_ij_ik = ep_ij[j][i] * ep_ik[k][i];
                const double scaled_exp_cos = lambda_epsilon * exp_ij_ik * cos_ijk_minus_cos0_sq;
                const double exp_term_ij = scaled_exp_cos * scaled_rma_ij_sq_rinv[j][i];
                const double exp_term_ik = scaled_exp_cos * scaled_rma_ik_sq_rinv[k][i];
                const double scaled_exp_cos_ij_ik = lambda_epsilon_2 * exp_ij_ik * cos_ijk_minus_cos0;
                const double ang_term = r_ij_ik_inv * scaled_exp_cos_ij_ik;
                const double cos_term_ij_ik = cos_ijk * scaled_exp_cos_ij_ik;
                const double cos_term_ij = r_ij_inv_sq[j][i] * cos_term_ij_ik;
                const double cos_term_ik = r_ik_inv_sq[k][i] * cos_term_ij_ik;

                fjx_lane[i] = dx_ij[0][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[0][k][i] * ang_term;
                fjy_lane[i] = dx_ij[1][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[1][k][i] * ang_term;
                fjz_lane[i] = dx_ij[2][j][i] * (exp_term_ij + cos_term_ij) - dx_ik[2][k][i] * ang_term;
                fkx_lane[i] = dx_ik[0][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[0][j][i] * ang_term;
                fky_lane[i] = dx_ik[1][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[1][j][i] * ang_term;
                fkz_lane[i] = dx_ik[2][k][i] * (exp_term_ik + cos_term_ik) - dx_ij[2][j][i] * ang_term;

                local_accum[0 * VEC_WIDTH + i] -= (fjx_lane[i] + fkx_lane[i]);
                local_accum[1 * VEC_WIDTH + i] -= (fjy_lane[i] + fky_lane[i]);
                local_accum[2 * VEC_WIDTH + i] -= (fjz_lane[i] + fkz_lane[i]);
            }
        }
    }

    for (int idx = 0; idx < 3 * VEC_WIDTH; ++idx) {
        force_i_block[idx] += local_accum[idx];
    }
}

void compute_force(ParticleDomain* P) {
    for (int ib = 0; ib < P->Ncore; ++ib) {
        const ParticleBlock& bi = P->blocks[ib];
        for (int jj = static_cast<int>(P->offset[ib]); jj < static_cast<int>(P->offset[ib + 1]); ++jj) {
            const int jb = static_cast<int>(P->neighbours[jj]);
            const ParticleBlock& bj = P->blocks[jb];
            for (int kk = jj; kk < static_cast<int>(P->offset[ib + 1]); ++kk) {
                const int kb = static_cast<int>(P->neighbours[kk]);
                const ParticleBlock& bk = P->blocks[kb];
                force_SW_triplet_8x8x8_tile(&bi, &bj, &bk, &(P->param), &(P->blocks[ib].acc[0][0]));
            }
        }
    }
}
