#include "init_phase.h"

/* ========================================================================
 *  init_phase.c — Initialization Phase for RSA Key Reconstruction
 *
 *  Mathematical background (HS09, Section 3):
 *
 *  From RSA: e·d = 1 + k·φ(N) = 1 + k·(N − p − q + 1)
 *  Since p+q is negligible compared to N (both are ~n/2 bits, N is ~n bits):
 *    d ≈ (k·(N+1) + 1) / e    for the correct k
 *  This approximation is exact in the top n/2 bits of d.
 *
 *  For CRT exponents:
 *    e·dp = 1 + kp·(p−1),  e·dq = 1 + kq·(q−1)
 *  Combining with the main equation gives:
 *    kp² − [k(N−1)+1]·kp − k ≡ 0  (mod e)
 *  which has exactly 2 solutions mod e (when e is prime): kp and kq.
 * ======================================================================== */

/* --- Init/Clear --- */

void init_result_init(init_result_t *res) {
    mpz_inits(res->k, res->kp, res->kq, NULL);
    mpz_inits(res->d_msb_corrected, NULL);
    mpz_inits(res->p0, res->q0, res->d0, res->dp0, res->dq0, NULL);
    res->tau_k = res->tau_kp = res->tau_kq = 0;
    res->valid = false;
    res->kp_kq_swapped = false;
}

void init_result_clear(init_result_t *res) {
    mpz_clears(res->k, res->kp, res->kq, NULL);
    mpz_clears(res->d_msb_corrected, NULL);
    mpz_clears(res->p0, res->q0, res->d0, res->dp0, res->dq0, NULL);
}

void init_result_print(const init_result_t *res) {
    printf("=== Init Phase Results ===\n");
    gmp_printf("  k    = %Zd\n", res->k);
    gmp_printf("  kp   = %Zd\n", res->kp);
    gmp_printf("  kq   = %Zd\n", res->kq);
    printf("  τ(k) = %d,  τ(kp) = %d,  τ(kq) = %d\n",
           res->tau_k, res->tau_kp, res->tau_kq);
    printf("  Valid: %s\n", res->valid ? "YES" : "NO");
    if (res->kp_kq_swapped)
        printf("  (kp/kq assignment swapped)\n");

    printf("\n  Initial partial values:\n");
    gmp_printf("    p0  = %Zd\n", res->p0);
    gmp_printf("    q0  = %Zd\n", res->q0);
    gmp_printf("    d0  = %Zd\n", res->d0);
    gmp_printf("    dp0 = %Zd\n", res->dp0);
    gmp_printf("    dq0 = %Zd\n", res->dq0);
    printf("  Slice(0): p[0]=%d, q[0]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n",
           get_bit(res->p0, 0), get_bit(res->q0, 0),
           res->tau_k, get_bit(res->d0, res->tau_k),
           res->tau_kp, get_bit(res->dp0, res->tau_kp),
           res->tau_kq, get_bit(res->dq0, res->tau_kq));
    printf("==========================\n");
}

/* ========================================================================
 *  Step 1: Find k
 *
 *  For each candidate k' ∈ [1, e-1]:
 *    d̃(k') = ⌊(k'·(N+1) + 1) / e⌋
 *  Compare the top (n/2 − 2) bits of d̃(k') with degraded d.
 *  The k' with the most matching MSBs is our k.
 *
 *  Why n/2 − 2 bits? Because p+q has at most n/2 + 1 bits, so
 *  d̃(k) and d agree on at least the top n/2 − 2 bits.
 * ======================================================================== */

int find_k(const mpz_t N, const mpz_t e, const mpz_t degraded_d,
           mpz_t out_k, mpz_t out_d_tilde, bool verbose)
{
    mpz_t k_candidate, d_tilde, temp, N_plus_1;
    mpz_inits(k_candidate, d_tilde, temp, N_plus_1, NULL);

    mpz_add_ui(N_plus_1, N, 1);

    size_t total_bits = mpz_sizeinbase(degraded_d, 2);
    size_t msb_count = (total_bits / 2) - 2;

    int best_match = 0;
    int best_match_k_count = 0;

    if (verbose) {
        printf("[INIT] Finding k: comparing top %zu bits of d̃(k') with degraded d\n",
               msb_count);
        printf("[INIT] Searching k' in [1, e-1] ...\n");
    }

    /* Iterate all candidates k' from 1 to e-1 */
    for (mpz_set_ui(k_candidate, 1);
         mpz_cmp(k_candidate, e) < 0;
         mpz_add_ui(k_candidate, k_candidate, 1))
    {
        /* d̃(k') = ⌊(k'·(N+1) + 1) / e⌋ */
        mpz_mul(temp, k_candidate, N_plus_1);
        mpz_add_ui(temp, temp, 1);
        mpz_fdiv_q(d_tilde, temp, e);

        /* Count matching MSBs */
        int matching = 0;
        for (size_t i = 0; i < msb_count; i++) {
            size_t bit_idx = total_bits - 1 - i;  /* From MSB downward */
            if (mpz_tstbit(d_tilde, bit_idx) == mpz_tstbit(degraded_d, bit_idx))
                matching++;
        }

        if (matching > best_match) {
            best_match = matching;
            mpz_set(out_k, k_candidate);
            mpz_set(out_d_tilde, d_tilde);
            best_match_k_count = 1;
        } else if (matching == best_match) {
            best_match_k_count++;
        }
    }

    if (verbose) {
        gmp_printf("[INIT] Best k = %Zd  (%d / %zu MSBs match, %.1f%%)\n",
                   out_k, best_match, msb_count,
                   100.0 * best_match / msb_count);
        if (best_match_k_count > 1) {
            printf("[INIT] WARNING: %d candidates tied for best match\n",
                   best_match_k_count);
        }
    }

    mpz_clears(k_candidate, d_tilde, temp, N_plus_1, NULL);
    return best_match;
}

/* ========================================================================
 *  Step 2: Correct MSB half of d
 *
 *  Take top (n/2 − 2) bits from d̃(k), bottom bits from degraded d.
 *  Result: d with corrected MSBs, degraded LSBs.
 * ======================================================================== */

void correct_d_msb(const mpz_t d_tilde, const mpz_t degraded_d,
                   mpz_t out_d, int n_bits)
{
    mpz_t msb_mask, lsb_mask, msb_part, lsb_part;
    mpz_inits(msb_mask, lsb_mask, msb_part, lsb_part, NULL);

    size_t total_bits = (size_t)n_bits;
    if (total_bits == 0)
        total_bits = mpz_sizeinbase(d_tilde, 2);

    size_t msb_bits = (total_bits / 2) - 2;
    size_t lsb_bits = total_bits - msb_bits;

    /* MSB mask: msb_bits ones shifted left by lsb_bits */
    mpz_set_ui(msb_mask, 1);
    mpz_mul_2exp(msb_mask, msb_mask, msb_bits);
    mpz_sub_ui(msb_mask, msb_mask, 1);
    mpz_mul_2exp(msb_mask, msb_mask, lsb_bits);

    /* LSB mask: lsb_bits ones */
    mpz_set_ui(lsb_mask, 1);
    mpz_mul_2exp(lsb_mask, lsb_mask, lsb_bits);
    mpz_sub_ui(lsb_mask, lsb_mask, 1);

    /* Combine: MSBs from d̃(k), LSBs from degraded d */
    mpz_and(msb_part, d_tilde, msb_mask);
    mpz_and(lsb_part, degraded_d, lsb_mask);
    mpz_ior(out_d, msb_part, lsb_part);

    mpz_clears(msb_mask, lsb_mask, msb_part, lsb_part, NULL);
}

/* ========================================================================
 *  Step 3: Tonelli-Shanks algorithm
 *
 *  Computes r such that r² ≡ n (mod p), where p is an odd prime.
 *  Returns both roots: r and p − r.
 *
 *  Algorithm:
 *   1. Write p−1 = Q · 2^S (Q odd)
 *   2. Find z: a quadratic non-residue mod p
 *   3. Set M=S, c=z^Q, t=n^Q, R=n^((Q+1)/2)
 *   4. Loop: find smallest i where t^(2^i) ≡ 1, update c, t, R, M
 *   5. When t=1: R is the square root
 * ======================================================================== */

int tonelli_shanks(const mpz_t n, const mpz_t p, mpz_t root1, mpz_t root2) {
    /* Check that n is a quadratic residue mod p */
    if (mpz_legendre(n, p) != 1) {
        return -1;  /* Not a QR, no solution */
    }

    mpz_t Q, z, M, c, t, R, temp, b, exponent;
    mpz_inits(Q, z, M, c, t, R, temp, b, exponent, NULL);

    /* Factor p−1 = Q · 2^S */
    mpz_sub_ui(Q, p, 1);
    unsigned long S = 0;
    while (mpz_even_p(Q)) {
        mpz_fdiv_q_2exp(Q, Q, 1);
        S++;
    }

    /* Find quadratic non-residue z */
    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) {
        mpz_add_ui(z, z, 1);
    }

    /* Initialize: c = z^Q mod p, t = n^Q mod p, R = n^((Q+1)/2) mod p */
    mpz_powm(c, z, Q, p);
    mpz_powm(t, n, Q, p);
    mpz_add_ui(temp, Q, 1);
    mpz_fdiv_q_2exp(temp, temp, 1);
    mpz_powm(R, n, temp, p);
    mpz_set_ui(M, S);

    /* Main loop */
    while (mpz_cmp_ui(t, 0) != 0 && mpz_cmp_ui(t, 1) != 0) {
        /* Find smallest i such that t^(2^i) ≡ 1 mod p */
        unsigned long i = 0;
        mpz_set(temp, t);
        while (mpz_cmp_ui(temp, 1) != 0) {
            mpz_powm_ui(temp, temp, 2, p);
            i++;
        }

        /* If i == M, no solution (shouldn't happen if n is QR) */
        if (i == mpz_get_ui(M)) {
            mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
            return -1;
        }

        /* b = c^(2^(M−i−1)) mod p */
        unsigned long exp_val = mpz_get_ui(M) - i - 1;
        mpz_ui_pow_ui(exponent, 2, exp_val);
        mpz_powm(b, c, exponent, p);

        /* Update: c = b², t = t·c, R = R·b, M = i */
        mpz_powm_ui(c, b, 2, p);
        mpz_mul(t, t, c);
        mpz_mod(t, t, p);
        mpz_mul(R, R, b);
        mpz_mod(R, R, p);
        mpz_set_ui(M, i);
    }

    /* Set solutions */
    if (mpz_cmp_ui(t, 1) == 0) {
        mpz_set(root1, R);
        mpz_sub(root2, p, R);  /* root2 = p − R */
    } else {
        mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
        return -1;
    }

    mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
    return 0;
}

/* ========================================================================
 *  Step 3b: Solve quadratic for kp, kq
 *
 *  Equation: kp² − [k(N−1)+1]·kp − k ≡ 0 (mod e)
 *
 *  Complete the square:
 *    Let A = [k(N−1) + 1 + e] / 2   (mod e)
 *    Then: (kp − A)² ≡ A² + k       (mod e)
 *
 *  So B = A² + k mod e, and we need √B mod e.
 *  Then kp = √B + A mod e,  kq = −√B + A mod e.
 * ======================================================================== */

/* ========================================================================
 * Step 3b: Solve quadratic for kp, kq
 *
 * Equation: kp² − [k(N−1)+1]·kp − k ≡ 0 (mod e)
 *
 * Complete the square:
 * Let A = [k(N−1) + 1 + e] / 2   (mod e)
 * Then: (kp − A)² ≡ A² + k       (mod e)
 *
 * So B = A² + k mod e, and we need √B mod e.
 * Then kp = √B + A mod e,  kq = −√B + A mod e.
 * ======================================================================== */

int solve_kp_kq(const mpz_t N, const mpz_t e, const mpz_t k,
                mpz_t kp, mpz_t kq)
{
    mpz_t A, B, sqrt1, sqrt2;
    mpz_inits(A, B, sqrt1, sqrt2, NULL);

    /* A = [k·(N−1) + 1 + e] / 2  mod e */
    mpz_sub_ui(A, N, 1);       /* A = N − 1          */
    mpz_mul(A, A, k);          /* A = k·(N−1)        */
    mpz_add_ui(A, A, 1);       /* A = k·(N−1) + 1    */

    /* FIX: Make A even before exact division by 2 */
    /* Since e is odd, adding e to an odd A makes it even. */
    if (mpz_odd_p(A)) {
        mpz_add(A, A, e);
    }

    mpz_divexact_ui(A, A, 2);  /* A = [...] / 2      */
    mpz_mod(A, A, e);          /* A = A mod e        */

    /* B = A² + k  mod e */
    mpz_mul(B, A, A);
    mpz_add(B, B, k);
    mpz_mod(B, B, e);

    /* Solve: x² ≡ B (mod e) using Tonelli-Shanks */
    if (tonelli_shanks(B, e, sqrt1, sqrt2) != 0) {
        fprintf(stderr, "[INIT] ERROR: No solution for kp/kq quadratic equation\n");
        mpz_clears(A, B, sqrt1, sqrt2, NULL);
        return -1;
    }

    /* kp = sqrt + A mod e,  kq = −sqrt + A mod e */
    mpz_add(kp, sqrt1, A);
    mpz_mod(kp, kp, e);

    mpz_add(kq, sqrt2, A);
    mpz_mod(kq, kq, e);

    mpz_clears(A, B, sqrt1, sqrt2, NULL);
    return 0;
}

/* ========================================================================
 *  Step 4: Correct LSBs of d, dp, dq
 *
 *  The lowest bits of d, dp, dq are deterministic:
 *    d  mod 2^(τ(k)+2)  = e⁻¹ mod 2^(τ(k)+2)
 *    dp mod 2^(τ(kp)+1) = e⁻¹ mod 2^(τ(kp)+1)
 *    dq mod 2^(τ(kq)+1) = e⁻¹ mod 2^(τ(kq)+1)
 *
 *  This is because:
 *    ed ≡ 1 (mod 2^(τ(k)+2))     [since k = 2^τ(k) · k' with k' odd]
 *    e·dp ≡ 1 (mod 2^(τ(kp)+1))  [similarly]
 * ======================================================================== */

void correct_lsb(mpz_t out, const mpz_t e, const mpz_t k_val,
                 int tau_val, bool is_d)
{
    mpz_t modulus;
    mpz_init(modulus);

    /* For d: power = τ(k) + 2,  for dp/dq: power = τ(kx) + 1 */
    int power = is_d ? (tau_val + 2) : (tau_val + 1);

    mpz_ui_pow_ui(modulus, 2, (unsigned long)power);
    mpz_invert(out, e, modulus);

    mpz_clear(modulus);
}

/* ========================================================================
 *  Verification (requires original key — for testing only)
 * ======================================================================== */

bool verify_init_result(const rsa_key_t *original, const init_result_t *result) {
    bool ok = true;
    mpz_t tmp, phi, p1, q1;
    mpz_inits(tmp, phi, p1, q1, NULL);

    /* Compute expected k: k = (ed - 1) / φ(N) */
    mpz_sub_ui(p1, original->p, 1);
    mpz_sub_ui(q1, original->q, 1);
    mpz_mul(phi, p1, q1);

    mpz_mul(tmp, original->e, original->d);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, phi);

    if (mpz_cmp(tmp, result->k) != 0) {
        gmp_printf("[VERIFY] FAIL: k mismatch. Expected %Zd, got %Zd\n",
                   tmp, result->k);
        ok = false;
    } else {
        printf("[VERIFY] k: OK\n");
    }

    /* Compute expected kp: kp = (e·dp - 1) / (p - 1) */
    mpz_mul(tmp, original->e, original->dp);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, p1);

    /* kp and kq might be swapped — check both assignments */
    if (mpz_cmp(tmp, result->kp) == 0) {
        printf("[VERIFY] kp: OK\n");
    } else if (mpz_cmp(tmp, result->kq) == 0) {
        printf("[VERIFY] kp: OK (swapped with kq)\n");
    } else {
        gmp_printf("[VERIFY] FAIL: kp mismatch. Expected %Zd, got kp=%Zd kq=%Zd\n",
                   tmp, result->kp, result->kq);
        ok = false;
    }

    /* Compute expected kq: kq = (e·dq - 1) / (q - 1) */
    mpz_mul(tmp, original->e, original->dq);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, q1);

    if (mpz_cmp(tmp, result->kp) == 0 || mpz_cmp(tmp, result->kq) == 0) {
        printf("[VERIFY] kq: OK\n");
    } else {
        gmp_printf("[VERIFY] FAIL: kq mismatch. Expected %Zd, got kp=%Zd kq=%Zd\n",
                   tmp, result->kp, result->kq);
        ok = false;
    }

    /* Verify τ values */
    int expected_tau_k = tau(result->k);
    if (result->tau_k != expected_tau_k) {
        printf("[VERIFY] FAIL: τ(k) = %d, expected %d\n",
               result->tau_k, expected_tau_k);
        ok = false;
    }

    mpz_clears(tmp, phi, p1, q1, NULL);
    return ok;
}

/* ========================================================================
 *  Main entry point: run full init phase
 * ======================================================================== */

int run_init_phase(const degraded_key_t *dk, init_result_t *result, bool verbose) {

    printf("╔══════════════════════════════════════╗\n");
    printf("║   Init Phase: Computing k, kp, kq    ║\n");
    printf("╚══════════════════════════════════════╝\n\n");

    /* ---- Step 1: Find k ---- */
    printf("[INIT] Step 1/4: Finding k from degraded d ...\n");

    mpz_t d_tilde;
    mpz_init(d_tilde);

    int match_count = find_k(dk->key.N, dk->key.e, dk->key.d,
                             result->k, d_tilde, verbose);

    if (match_count == 0) {
        fprintf(stderr, "[INIT] ERROR: Could not find k (no MSB matches)\n");
        mpz_clear(d_tilde);
        return -1;
    }

    gmp_printf("[INIT] Found k = %Zd  (%d MSBs matched)\n\n", result->k, match_count);

    /* ---- Step 2: Correct MSB half of d ---- */
    printf("[INIT] Step 2/4: Correcting MSB half of d ...\n");

    correct_d_msb(d_tilde, dk->key.d, result->d_msb_corrected, dk->full_bits);

    if (verbose) {
        gmp_printf("[INIT] d̃(k)           = %Zd\n", d_tilde);
        gmp_printf("[INIT] d_msb_corrected = %Zd\n", result->d_msb_corrected);
    }
    printf("[INIT] MSB correction done\n\n");

    mpz_clear(d_tilde);

    /* ---- Step 3: Solve for kp, kq ---- */
    printf("[INIT] Step 3/4: Solving quadratic for kp, kq ...\n");

    if (solve_kp_kq(dk->key.N, dk->key.e, result->k,
                    result->kp, result->kq) != 0) {
        fprintf(stderr, "[INIT] ERROR: Failed to solve for kp/kq\n");
        return -1;
    }

    gmp_printf("[INIT] kp = %Zd\n", result->kp);
    gmp_printf("[INIT] kq = %Zd\n\n", result->kq);

    /* ---- Step 4: Compute τ values and correct LSBs ---- */
    printf("[INIT] Step 4/4: Computing τ values, correcting LSBs ...\n");

    result->tau_k  = tau(result->k);
    result->tau_kp = tau(result->kp);
    result->tau_kq = tau(result->kq);

    printf("[INIT] τ(k)=%d, τ(kp)=%d, τ(kq)=%d\n",
           result->tau_k, result->tau_kp, result->tau_kq);

    /* Initialize p0, q0 with known bit 0 (primes are odd) */
    mpz_set_ui(result->p0, 0);
    mpz_set_ui(result->q0, 0);
    mpz_setbit(result->p0, 0);  /* p[0] = 1 */
    mpz_setbit(result->q0, 0);  /* q[0] = 1 */

    /* Correct LSBs of d, dp, dq */
    correct_lsb(result->d0,  dk->key.e, result->k,  result->tau_k,  true);
    correct_lsb(result->dp0, dk->key.e, result->kp, result->tau_kp, false);
    correct_lsb(result->dq0, dk->key.e, result->kq, result->tau_kq, false);

    if (verbose) {
        printf("[INIT] LSB corrections:\n");
        gmp_printf("  d0  = %Zd  (mod 2^%d)\n", result->d0, result->tau_k + 2);
        gmp_printf("  dp0 = %Zd  (mod 2^%d)\n", result->dp0, result->tau_kp + 1);
        gmp_printf("  dq0 = %Zd  (mod 2^%d)\n", result->dq0, result->tau_kq + 1);
    }

    result->valid = true;
    printf("\n[INIT] Initialization complete.\n");

    if (verbose) {
        init_result_print(result);
    }

    return 0;
}
