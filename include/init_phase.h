#ifndef INIT_PHASE_H
#define INIT_PHASE_H

#include "common.h"

/* ========================================================================
 *  init_phase.h — Initialization Phase for RSA Key Reconstruction
 *
 *  Shared by both Heninger-Shacham and Henecka-May-Meurer algorithms.
 *
 *  Steps:
 *   1. Find k from degraded d by comparing MSBs of d̃(k') with degraded d
 *   2. Correct the MSB half of d using d̃(k)
 *   3. Solve quadratic equation mod e to find kp, kq (via Tonelli-Shanks)
 *   4. Compute τ(k), τ(kp), τ(kq)
 *   5. Correct LSBs of d, dp, dq using modular inverses
 *   6. Set initial known bits: p[0]=1, q[0]=1 (primes are odd)
 *
 *  Background (from HS09 paper):
 *   - ed = 1 + k(N - (p+q) + 1)  where k ∈ [1, e-1]
 *   - d̃(k') = ⌊(k'(N+1)+1)/e⌋ approximates d when k'=k
 *     (because p+q ≈ 0 relative to N, the MSBs of d̃(k) match d)
 *   - k² − [k(N−1)+1]kp − k ≡ 0 (mod e)  has exactly 2 solutions
 *     when e is prime: one is kp, the other is kq
 *   - τ(x) = v₂(x) = largest power of 2 dividing x
 *   - LSBs of d, dp, dq are determined: e⁻¹ mod 2^(τ+offset)
 * ======================================================================== */

/* --- Result structure for the init phase --- */
typedef struct {
    mpz_t k;            /* Secret multiplier: ed = 1 + k·φ(N) */
    mpz_t kp;           /* CRT multiplier: e·dp = 1 + kp·(p-1) */
    mpz_t kq;           /* CRT multiplier: e·dq = 1 + kq·(q-1) */

    int tau_k;          /* τ(k)  = v₂(k)  */
    int tau_kp;         /* τ(kp) = v₂(kp) */
    int tau_kq;         /* τ(kq) = v₂(kq) */

    mpz_t d_msb_corrected;  /* d with MSB half corrected via d̃(k) */

    /* Initial partial values for branch-and-prune (bits 0..τ corrected) */
    mpz_t p0;           /* Partial p:  bit 0 = 1 (odd prime) */
    mpz_t q0;           /* Partial q:  bit 0 = 1 (odd prime) */
    mpz_t d0;           /* Partial d:  LSBs corrected up to τ(k)+1 */
    mpz_t dp0;          /* Partial dp: LSBs corrected up to τ(kp) */
    mpz_t dq0;          /* Partial dq: LSBs corrected up to τ(kq) */

    bool valid;         /* true if initialization succeeded */
    bool kp_kq_swapped; /* true if we had to try swapped kp/kq assignment */
} init_result_t;

/* Initialize/clear the result structure */
void init_result_init(init_result_t *res);
void init_result_clear(init_result_t *res);

/* Print a summary of the init phase results */
void init_result_print(const init_result_t *res);

/* --- Main init phase entry point --- */

/*
 * Run the full initialization phase on a degraded key.
 *
 * @param dk       The degraded key (with known-bits arrays)
 * @param result   Output: populated init_result_t
 * @param verbose  Print detailed progress information
 * @return         0 on success, -1 on failure
 */
int run_init_phase(const degraded_key_t *dk, init_result_t *result, bool verbose);

/* --- Individual steps (exposed for testing) --- */

/*
 * Step 1: Find k by comparing d̃(k') MSBs with degraded d.
 *
 * Iterates k' from 1 to e-1, computes d̃(k') = ⌊(k'(N+1)+1)/e⌋,
 * counts matching MSBs with degraded d. Best match gives k.
 *
 * @param N           RSA modulus
 * @param e           Public exponent
 * @param degraded_d  Degraded private exponent
 * @param out_k       Output: found k value
 * @param out_d_tilde Output: d̃(k) for the winning k
 * @param verbose     Print progress
 * @return            Number of matching MSBs (higher = more confident)
 */
int find_k(const mpz_t N, const mpz_t e, const mpz_t degraded_d,
           mpz_t out_k, mpz_t out_d_tilde, bool verbose);

/*
 * Step 2: Correct MSB half of degraded d.
 *
 * Takes MSB n/2-2 bits from d̃(k), LSB bits from degraded d.
 *
 * @param d_tilde     d̃(k) computed from the found k
 * @param degraded_d  Original degraded d
 * @param out_d       Output: d with corrected MSBs
 * @param n_bits      Total bit length of d
 */
void correct_d_msb(const mpz_t d_tilde, const mpz_t degraded_d,
                   mpz_t out_d, int n_bits);

/*
 * Step 3: Solve for kp and kq.
 *
 * Solves: kp² − [k(N−1)+1]·kp − k ≡ 0 (mod e)
 * Rewritten as: (kp − A)² ≡ A² + k (mod e)
 * where A = [k(N−1)+1+e] / 2 mod e
 *
 * Uses Tonelli-Shanks for modular square root.
 *
 * @param N     RSA modulus
 * @param e     Public exponent (must be prime)
 * @param k     Secret multiplier found in step 1
 * @param kp    Output: first solution (candidate kp)
 * @param kq    Output: second solution (candidate kq)
 * @return      0 on success, -1 if no solution
 */
int solve_kp_kq(const mpz_t N, const mpz_t e, const mpz_t k,
                mpz_t kp, mpz_t kq);

/*
 * Tonelli-Shanks algorithm: compute square root of n mod p.
 *
 * @param n     Value to take square root of
 * @param p     Prime modulus
 * @param root1 Output: first square root
 * @param root2 Output: second square root (= p - root1)
 * @return      0 on success, -1 if n is not a QR mod p
 */
int tonelli_shanks(const mpz_t n, const mpz_t p, mpz_t root1, mpz_t root2);

/*
 * Step 4: Correct LSBs of d, dp, dq.
 *
 * The lowest τ(k)+1 bits of d (and τ(kp), τ(kq) bits of dp, dq)
 * are deterministic: they equal e⁻¹ mod 2^power.
 *
 * @param out     Output: value with corrected LSBs
 * @param e       Public exponent
 * @param k_val   k, kp, or kq respectively
 * @param tau_val τ(k), τ(kp), or τ(kq)
 * @param is_d    true for d (uses power=τ+2), false for dp/dq (uses power=τ+1)
 */
void correct_lsb(mpz_t out, const mpz_t e, const mpz_t k_val,
                 int tau_val, bool is_d);

/*
 * Verify that k, kp, kq are consistent with the original key.
 * (Only usable for testing — requires the original key.)
 *
 * @param original  Original (non-degraded) RSA key
 * @param result    Init phase results to verify
 * @return          true if all values match
 */
bool verify_init_result(const rsa_key_t *original, const init_result_t *result);

#endif /* INIT_PHASE_H */
