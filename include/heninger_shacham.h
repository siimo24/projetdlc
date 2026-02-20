#ifndef HENINGER_SHACHAM_H
#define HENINGER_SHACHAM_H

#include "common.h"
#include "init_phase.h"

/* ========================================================================
 *  heninger_shacham.h — Branch & Prune RSA Key Reconstruction (HS09)
 *
 *  Erasure correction: works when a fraction δ ≥ 0.27 of bits are known.
 *
 *  Algorithm (HS09, Section 4):
 *   For each bit position i from 1 to n/2:
 *     - Slice(i) = (p[i], q[i], d[i+τ(k)], dp[i+τ(kp)], dq[i+τ(kq)])
 *     - 4 equations constrain the 5 unknowns → exactly 2 candidates
 *     - Prune candidates that disagree with known bits
 *     - Recurse on surviving candidates
 *   Termination:
 *     - Compute p_candidate = (e·dp − 1)/kp + 1
 *     - If p_candidate · q_candidate = N → success
 * ======================================================================== */

/* --- Reconstruction result --- */
typedef struct {
    rsa_key_t recovered_key;    /* The fully recovered RSA key */
    bool success;               /* true if key was recovered */

    /* Statistics */
    long long branches_explored;    /* Total recursive calls */
    long long branches_pruned;      /* Pruned by equations */
    long long branches_known_pruned;/* Pruned by known bits */
    int max_depth_reached;          /* Deepest recursion level */
    double elapsed_seconds;         /* Wall clock time */
} hs_result_t;

void hs_result_init(hs_result_t *res);
void hs_result_clear(hs_result_t *res);
void hs_result_print(const hs_result_t *res);

/* --- Main entry point --- */

/*
 * Run the Heninger-Shacham branch-and-prune reconstruction.
 *
 * Tries the given kp/kq assignment first. If it fails and try_swap is true,
 * automatically retries with kp and kq swapped.
 *
 * @param dk          Degraded key (with known-bits arrays)
 * @param init        Init phase results (k, kp, kq, τ values, partial values)
 * @param result      Output: recovered key and statistics
 * @param try_swap    If true, retry with swapped kp/kq on failure
 * @param verbose     Print per-slice progress
 * @return            0 on success, -1 on failure
 */
int run_heninger_shacham(const degraded_key_t *dk,
                         const init_result_t *init,
                         hs_result_t *result,
                         bool try_swap,
                         bool verbose);

#endif /* HENINGER_SHACHAM_H */
