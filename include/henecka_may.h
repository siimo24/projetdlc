#ifndef HENECKA_MAY_H
#define HENECKA_MAY_H

#include "common.h"
#include "init_phase.h"

/* ========================================================================
 *  henecka_may.h — Block-Threshold RSA Key Reconstruction (HMM10)
 *
 *  Error correction: works when bits may be WRONG (δ < 0.237 error rate),
 *  unlike HS09 which requires knowing which bits are correct.
 *
 *  Algorithm overview (HMM10):
 *   - Same equations (8)-(11) as HS09 constrain each bit slice
 *   - But instead of pruning with known bits, process blocks of t bits
 *   - Within each block: expand all 2^t candidate paths
 *   - For each path: count how many bits match the degraded key
 *   - Threshold prune: keep only paths with match_count ≥ C
 *
 *  Key parameters:
 *   t (block size):  Θ(ln(n) / ε²), typically 10-20 for 1024-bit keys
 *   C (threshold):   ⌊5t · (1/2 + γ₀)⌋  where γ₀ ≈ 0.263
 *
 *  The correct candidate matches ~5t·(1−δ) bits (high),
 *  wrong candidates match ~5t/2 bits (random). C separates these.
 * ======================================================================== */

/* --- Configuration parameters --- */
typedef struct {
    int block_size_t;       /* Block size t (0 = auto-compute from delta) */
    int threshold_C;        /* Threshold C (0 = auto-compute from t) */
    int max_candidates;     /* Maximum candidates before forced pruning */
    double target_delta;    /* Target error rate for auto parameter selection */
} hmm_params_t;

/* Initialize with sensible defaults */
void hmm_params_default(hmm_params_t *params);

/* Auto-compute t and C from error rate and key size */
void hmm_params_auto(hmm_params_t *params, int key_bits, double error_rate);

/* Print parameters */
void hmm_params_print(const hmm_params_t *params);

/* --- Reconstruction result --- */
typedef struct {
    rsa_key_t recovered_key;
    bool success;

    /* Statistics */
    long long total_candidates_generated;
    long long total_candidates_pruned;
    int max_candidates_at_once;
    int blocks_processed;
    double elapsed_seconds;
} hmm_result_t;

void hmm_result_init(hmm_result_t *res);
void hmm_result_clear(hmm_result_t *res);
void hmm_result_print(const hmm_result_t *res);

/* --- Main entry point --- */

/*
 * Run the Henecka-May-Meurer block-threshold reconstruction.
 *
 * @param dk          Degraded key (uses raw degraded values, NOT known-bits)
 * @param init        Init phase results (k, kp, kq, τ values)
 * @param params      Algorithm parameters (block size, threshold)
 * @param result      Output: recovered key and statistics
 * @param try_swap    If true, retry with swapped kp/kq on failure
 * @param verbose     Print per-block progress
 * @return            0 on success, -1 on failure
 */
int run_henecka_may(const degraded_key_t *dk,
                    const init_result_t *init,
                    const hmm_params_t *params,
                    hmm_result_t *result,
                    bool try_swap,
                    bool verbose);

#endif /* HENECKA_MAY_H */
