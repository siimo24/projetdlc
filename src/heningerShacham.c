#include "heninger_shacham.h"
#include <time.h>

/* ========================================================================
 *  heninger_shacham.c — Branch & Prune RSA Key Reconstruction
 *
 *  From HS09 paper, equations (8)–(11):
 *
 *  At bit position i, with current partial values p0, q0, d0, dp0, dq0,
 *  and candidate bits p_i, q_i, d_i, dp_i, dq_i:
 *
 *  (8)  bit_i(N − p0·q0)                        = p_i ⊕ q_i
 *  (9)  bit_{i+τ(k)}(k(N+1) − k(p0+q0) − e·d0 + 1) = d_i ⊕ p_i ⊕ q_i
 *  (10) bit_{i+τ(kp)}(kp(p0−1) + 1 − e·dp0)    = dp_i ⊕ p_i
 *  (11) bit_{i+τ(kq)}(kq(q0−1) + 1 − e·dq0)    = dq_i ⊕ q_i
 *
 *  These 4 equations on 5 unknowns yield exactly 2 valid candidates
 *  per slice (in the absence of known-bit pruning). Known bits further
 *  reduce this to ≤ 2, often to 1, keeping the tree narrow.
 *
 *  Termination check at each level:
 *    p_cand = (e·dp0 − 1) / kp + 1
 *    q_cand = (e·dq0 − 1) / kq + 1
 *    If p_cand · q_cand = N → found the factorization.
 * ======================================================================== */

/* --- Result init/clear --- */

void hs_result_init(hs_result_t *res) {
    rsa_key_init(&res->recovered_key);
    res->success = false;
    res->branches_explored = 0;
    res->branches_pruned = 0;
    res->branches_known_pruned = 0;
    res->max_depth_reached = 0;
    res->elapsed_seconds = 0.0;
}

void hs_result_clear(hs_result_t *res) {
    rsa_key_clear(&res->recovered_key);
}

void hs_result_print(const hs_result_t *res) {
    printf("=== Heninger-Shacham Results ===\n");
    printf("  Success:              %s\n", res->success ? "YES" : "NO");
    printf("  Branches explored:    %lld\n", res->branches_explored);
    printf("  Pruned (equations):   %lld\n", res->branches_pruned);
    printf("  Pruned (known bits):  %lld\n", res->branches_known_pruned);
    printf("  Max depth:            %d\n", res->max_depth_reached);
    printf("  Time:                 %.3f s\n", res->elapsed_seconds);
    if (res->success) {
        gmp_printf("  p = %Zd\n", res->recovered_key.p);
        gmp_printf("  q = %Zd\n", res->recovered_key.q);
    }
    printf("================================\n");
}

/* ---- All 32 possibilities for a 5-bit slice ---- */
static const char POSSIBILITIES[32][5] = {
    {0,0,0,0,0}, {0,0,0,0,1}, {0,0,0,1,0}, {0,0,0,1,1},
    {0,0,1,0,0}, {0,0,1,0,1}, {0,0,1,1,0}, {0,0,1,1,1},
    {0,1,0,0,0}, {0,1,0,0,1}, {0,1,0,1,0}, {0,1,0,1,1},
    {0,1,1,0,0}, {0,1,1,0,1}, {0,1,1,1,0}, {0,1,1,1,1},
    {1,0,0,0,0}, {1,0,0,0,1}, {1,0,0,1,0}, {1,0,0,1,1},
    {1,0,1,0,0}, {1,0,1,0,1}, {1,0,1,1,0}, {1,0,1,1,1},
    {1,1,0,0,0}, {1,1,0,0,1}, {1,1,0,1,0}, {1,1,0,1,1},
    {1,1,1,0,0}, {1,1,1,0,1}, {1,1,1,1,0}, {1,1,1,1,1}
};

/* ---- Equation checks ---- */

/* Equation (8): bit_i(N − p0·q0) = p_i ⊕ q_i */
static bool check_eq8(const mpz_t N, const mpz_t p0, const mpz_t q0,
                      char p_i, char q_i, int i)
{
    mpz_t term;
    mpz_init(term);
    mpz_mul(term, p0, q0);
    mpz_sub(term, N, term);
    bool ok = (get_bit(term, i) == (p_i ^ q_i));
    mpz_clear(term);
    return ok;
}

/* Equation (9): bit_{i+τ(k)}(k(N+1) − k(p0+q0) − e·d0 + 1) = d_i ⊕ p_i ⊕ q_i */
static bool check_eq9(const mpz_t N, const mpz_t e, const mpz_t k, int tau_k,
                      const mpz_t p0, const mpz_t q0, const mpz_t d0,
                      char p_i, char q_i, char d_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_add_ui(temp, N, 1);
    mpz_mul(term, k, temp);         /* term = k·(N+1) */

    mpz_add(temp, p0, q0);
    mpz_mul(temp, k, temp);
    mpz_sub(term, term, temp);      /* term -= k·(p0+q0) */

    mpz_mul(temp, e, d0);
    mpz_sub(term, term, temp);      /* term -= e·d0 */

    mpz_add_ui(term, term, 1);     /* term += 1 */

    bool ok = (get_bit(term, i + tau_k) == (d_i ^ p_i ^ q_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

/* Equation (10): bit_{i+τ(kp)}(kp(p0−1) + 1 − e·dp0) = dp_i ⊕ p_i */
static bool check_eq10(const mpz_t e, const mpz_t kp, int tau_kp,
                       const mpz_t p0, const mpz_t dp0,
                       char p_i, char dp_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_sub_ui(temp, p0, 1);
    mpz_mul(term, kp, temp);       /* term = kp·(p0−1) */
    mpz_add_ui(term, term, 1);     /* term += 1 */
    mpz_mul(temp, e, dp0);
    mpz_sub(term, term, temp);     /* term -= e·dp0 */

    bool ok = (get_bit(term, i + tau_kp) == (dp_i ^ p_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

/* Equation (11): bit_{i+τ(kq)}(kq(q0−1) + 1 − e·dq0) = dq_i ⊕ q_i */
static bool check_eq11(const mpz_t e, const mpz_t kq, int tau_kq,
                       const mpz_t q0, const mpz_t dq0,
                       char q_i, char dq_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_sub_ui(temp, q0, 1);
    mpz_mul(term, kq, temp);       /* term = kq·(q0−1) */
    mpz_add_ui(term, term, 1);     /* term += 1 */
    mpz_mul(temp, e, dq0);
    mpz_sub(term, term, temp);     /* term -= e·dq0 */

    bool ok = (get_bit(term, i + tau_kq) == (dq_i ^ q_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

/* ---- Termination check ---- */

/*
 * Try to recover p, q from current partial dp, dq:
 *   p = (e·dp − 1) / kp + 1
 *   q = (e·dq − 1) / kq + 1
 * Check if p·q = N.
 */
static bool try_termination(const mpz_t e, const mpz_t kp, const mpz_t kq,
                            const mpz_t N, const mpz_t dp0, const mpz_t dq0,
                            mpz_t out_p, mpz_t out_q)
{
    mpz_t p_cand, q_cand, product, temp;
    mpz_inits(p_cand, q_cand, product, temp, NULL);

    /* p = (e·dp − 1) / kp + 1 */
    mpz_mul(p_cand, e, dp0);
    mpz_sub_ui(p_cand, p_cand, 1);

    /* Check divisibility before dividing */
    if (!mpz_divisible_p(p_cand, kp)) {
        mpz_clears(p_cand, q_cand, product, temp, NULL);
        return false;
    }
    mpz_divexact(p_cand, p_cand, kp);
    mpz_add_ui(p_cand, p_cand, 1);

    /* q = (e·dq − 1) / kq + 1 */
    mpz_mul(q_cand, e, dq0);
    mpz_sub_ui(q_cand, q_cand, 1);

    if (!mpz_divisible_p(q_cand, kq)) {
        mpz_clears(p_cand, q_cand, product, temp, NULL);
        return false;
    }
    mpz_divexact(q_cand, q_cand, kq);
    mpz_add_ui(q_cand, q_cand, 1);

    /* Check p·q = N */
    mpz_mul(product, p_cand, q_cand);
    bool found = (mpz_cmp(product, N) == 0);

    if (found) {
        mpz_set(out_p, p_cand);
        mpz_set(out_q, q_cand);
    }

    mpz_clears(p_cand, q_cand, product, temp, NULL);
    return found;
}

/* ---- Recursive Branch & Prune ---- */

/* Valid slice: up to 2 solutions per level */
typedef struct {
    char bits[5];   /* p_i, q_i, d_i, dp_i, dq_i */
} slice_t;

/*
 * Core recursive function.
 *
 * Parameters are passed via a context struct to keep the signature manageable
 * and avoid passing 20+ arguments through each recursive call.
 */
typedef struct {
    /* Constants (shared across all recursion levels) */
    const mpz_t *N;
    const mpz_t *e;
    const mpz_t *k;
    const mpz_t *kp;
    const mpz_t *kq;
    int tau_k, tau_kp, tau_kq;
    int half_bits;          /* n/2: number of bits in p, q */

    /* Known bits arrays */
    const int *known_p;
    const int *known_q;
    const int *known_d;
    const int *known_dp;
    const int *known_dq;

    /* Mutable state */
    hs_result_t *result;
    bool verbose;
} hs_context_t;

static void branch_and_prune(
    hs_context_t *ctx,
    mpz_t my_p, mpz_t my_q, mpz_t my_d, mpz_t my_dp, mpz_t my_dq,
    int depth)
{
    /* Early exit if solution already found */
    if (ctx->result->success) return;

    ctx->result->branches_explored++;

    if (depth > ctx->result->max_depth_reached)
        ctx->result->max_depth_reached = depth;

    /* Progress reporting every 50 slices */
    if (ctx->verbose && depth % 50 == 0 && depth > 0) {
        printf("  [HS] Depth %d/%d, branches explored: %lld\n",
               depth, ctx->half_bits, ctx->result->branches_explored);
    }

    /* ---- Try termination at each level ---- */
    mpz_t cand_p, cand_q;
    mpz_inits(cand_p, cand_q, NULL);

    if (try_termination(*ctx->e, *ctx->kp, *ctx->kq, *ctx->N,
                        my_dp, my_dq, cand_p, cand_q))
    {
        /* SUCCESS — factorization found! */
        ctx->result->success = true;

        mpz_set(ctx->result->recovered_key.p, cand_p);
        mpz_set(ctx->result->recovered_key.q, cand_q);
        mpz_set(ctx->result->recovered_key.N, *ctx->N);
        mpz_set(ctx->result->recovered_key.e, *ctx->e);

        /* Compute d = e⁻¹ mod φ(N) */
        mpz_t p1, q1, phi;
        mpz_inits(p1, q1, phi, NULL);
        mpz_sub_ui(p1, cand_p, 1);
        mpz_sub_ui(q1, cand_q, 1);
        mpz_mul(phi, p1, q1);
        mpz_invert(ctx->result->recovered_key.d, *ctx->e, phi);

        /* Compute dp, dq, qinv */
        mpz_mod(ctx->result->recovered_key.dp, ctx->result->recovered_key.d, p1);
        mpz_mod(ctx->result->recovered_key.dq, ctx->result->recovered_key.d, q1);
        mpz_invert(ctx->result->recovered_key.qinv, cand_q, cand_p);

        ctx->result->recovered_key.bits = ctx->half_bits * 2;

        mpz_clears(p1, q1, phi, NULL);
        mpz_clears(cand_p, cand_q, NULL);
        return;
    }

    mpz_clears(cand_p, cand_q, NULL);

    /* If we've gone past the half-bit boundary, stop */
    if (depth >= ctx->half_bits) return;

    /* ---- Test all 32 candidate slices ---- */
    slice_t valid[2];
    int valid_count = 0;

    for (int s = 0; s < 32; s++) {
        char p_i  = POSSIBILITIES[s][0];
        char q_i  = POSSIBILITIES[s][1];
        char d_i  = POSSIBILITIES[s][2];
        char dp_i = POSSIBILITIES[s][3];
        char dq_i = POSSIBILITIES[s][4];

        /* Check all 4 equations */
        if (!check_eq8(*ctx->N, my_p, my_q, p_i, q_i, depth)) {
            ctx->result->branches_pruned++;
            continue;
        }
        if (!check_eq9(*ctx->N, *ctx->e, *ctx->k, ctx->tau_k,
                       my_p, my_q, my_d, p_i, q_i, d_i, depth)) {
            ctx->result->branches_pruned++;
            continue;
        }
        if (!check_eq10(*ctx->e, *ctx->kp, ctx->tau_kp,
                        my_p, my_dp, p_i, dp_i, depth)) {
            ctx->result->branches_pruned++;
            continue;
        }
        if (!check_eq11(*ctx->e, *ctx->kq, ctx->tau_kq,
                        my_q, my_dq, q_i, dq_i, depth)) {
            ctx->result->branches_pruned++;
            continue;
        }

        /* Check against known bits */
        bool matches_known = true;

        if (ctx->known_p[depth] != BIT_UNKNOWN &&
            ctx->known_p[depth] != p_i)
            matches_known = false;

        if (matches_known &&
            ctx->known_q[depth] != BIT_UNKNOWN &&
            ctx->known_q[depth] != q_i)
            matches_known = false;

        if (matches_known &&
            (depth + ctx->tau_k) < ctx->half_bits * 2 &&
            ctx->known_d[depth + ctx->tau_k] != BIT_UNKNOWN &&
            ctx->known_d[depth + ctx->tau_k] != d_i)
            matches_known = false;

        if (matches_known &&
            (depth + ctx->tau_kp) < ctx->half_bits &&
            ctx->known_dp[depth + ctx->tau_kp] != BIT_UNKNOWN &&
            ctx->known_dp[depth + ctx->tau_kp] != dp_i)
            matches_known = false;

        if (matches_known &&
            (depth + ctx->tau_kq) < ctx->half_bits &&
            ctx->known_dq[depth + ctx->tau_kq] != BIT_UNKNOWN &&
            ctx->known_dq[depth + ctx->tau_kq] != dq_i)
            matches_known = false;

        if (!matches_known) {
            ctx->result->branches_known_pruned++;
            continue;
        }

        /* This slice is valid — store it */
        if (valid_count < 2) {
            memcpy(valid[valid_count].bits, POSSIBILITIES[s], 5);
            valid_count++;
        }
        /* Theoretically at most 2 valid slices from equations alone;
         * with known-bit pruning it can be 0, 1, or 2 */
    }

    /* ---- Recurse on valid slices ---- */
    for (int v = 0; v < valid_count; v++) {
        if (ctx->result->success) return;

        /* Clone partial values for this branch */
        mpz_t cl_p, cl_q, cl_d, cl_dp, cl_dq;
        mpz_init_set(cl_p,  my_p);
        mpz_init_set(cl_q,  my_q);
        mpz_init_set(cl_d,  my_d);
        mpz_init_set(cl_dp, my_dp);
        mpz_init_set(cl_dq, my_dq);

        /* Set the new bits */
        if (valid[v].bits[0]) mpz_setbit(cl_p,  depth);
        if (valid[v].bits[1]) mpz_setbit(cl_q,  depth);
        if (valid[v].bits[2]) mpz_setbit(cl_d,  depth + ctx->tau_k);
        if (valid[v].bits[3]) mpz_setbit(cl_dp, depth + ctx->tau_kp);
        if (valid[v].bits[4]) mpz_setbit(cl_dq, depth + ctx->tau_kq);

        /* Recurse */
        branch_and_prune(ctx, cl_p, cl_q, cl_d, cl_dp, cl_dq, depth + 1);

        mpz_clears(cl_p, cl_q, cl_d, cl_dp, cl_dq, NULL);
    }
}

/* ---- Single run with a given kp/kq assignment ---- */

static int hs_run_once(const degraded_key_t *dk,
                       const mpz_t k, const mpz_t kp, const mpz_t kq,
                       int tau_k, int tau_kp, int tau_kq,
                       hs_result_t *result,
                       bool verbose)
{
    int half_bits = dk->half_bits;

    /* Build initial partial values */
    mpz_t p0, q0, d0, dp0, dq0;
    mpz_inits(p0, q0, d0, dp0, dq0, NULL);

    /* p[0] = q[0] = 1 (primes are odd) */
    mpz_set_ui(p0, 0);
    mpz_set_ui(q0, 0);
    mpz_setbit(p0, 0);
    mpz_setbit(q0, 0);

    /* Correct LSBs of d, dp, dq */
    /* d:  e⁻¹ mod 2^(τ(k)+2)  */
    /* dp: e⁻¹ mod 2^(τ(kp)+1) */
    /* dq: e⁻¹ mod 2^(τ(kq)+1) */
    mpz_t modulus;
    mpz_init(modulus);

    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_k + 2));
    mpz_invert(d0, dk->key.e, modulus);

    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kp + 1));
    mpz_invert(dp0, dk->key.e, modulus);

    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kq + 1));
    mpz_invert(dq0, dk->key.e, modulus);

    mpz_clear(modulus);

    if (verbose) {
        gmp_printf("[HS] Starting branch & prune (half_bits=%d)\n", half_bits);
        printf("[HS] Slice(0): p[0]=%d, q[0]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n",
               get_bit(p0, 0), get_bit(q0, 0),
               tau_k, get_bit(d0, tau_k),
               tau_kp, get_bit(dp0, tau_kp),
               tau_kq, get_bit(dq0, tau_kq));
    }

    /* Set up context */
    hs_context_t ctx;
    ctx.N = &dk->key.N;
    ctx.e = &dk->key.e;
    ctx.k = &k;
    ctx.kp = &kp;
    ctx.kq = &kq;
    ctx.tau_k = tau_k;
    ctx.tau_kp = tau_kp;
    ctx.tau_kq = tau_kq;
    ctx.half_bits = half_bits;
    ctx.known_p  = dk->known_p;
    ctx.known_q  = dk->known_q;
    ctx.known_d  = dk->known_d;
    ctx.known_dp = dk->known_dp;
    ctx.known_dq = dk->known_dq;
    ctx.result = result;
    ctx.verbose = verbose;

    /* Launch recursion starting at bit 1 (bit 0 already set) */
    branch_and_prune(&ctx, p0, q0, d0, dp0, dq0, 1);

    mpz_clears(p0, q0, d0, dp0, dq0, NULL);

    return result->success ? 0 : -1;
}

/* ---- Public entry point ---- */

int run_heninger_shacham(const degraded_key_t *dk,
                         const init_result_t *init,
                         hs_result_t *result,
                         bool try_swap,
                         bool verbose)
{
    printf("╔══════════════════════════════════════════════════╗\n");
    printf("║   Heninger-Shacham Branch & Prune (HS09)         ║\n");
    printf("╚══════════════════════════════════════════════════╝\n\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    /* ---- Attempt 1: kp, kq as given ---- */
    gmp_printf("[HS] Attempt 1: kp=%Zd, kq=%Zd\n", init->kp, init->kq);
    printf("[HS] τ(k)=%d, τ(kp)=%d, τ(kq)=%d\n\n",
           init->tau_k, init->tau_kp, init->tau_kq);

    int ret = hs_run_once(dk, init->k, init->kp, init->kq,
                          init->tau_k, init->tau_kp, init->tau_kq,
                          result, verbose);

    if (ret == 0) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                                + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("\n[HS] SUCCESS on attempt 1!\n");
        hs_result_print(result);
        return 0;
    }

    /* ---- Attempt 2: swap kp and kq ---- */
    if (try_swap) {
        printf("\n[HS] Attempt 1 failed (%lld branches explored).\n",
               result->branches_explored);
        gmp_printf("[HS] Attempt 2: swapping kp=%Zd, kq=%Zd\n",
                   init->kq, init->kp);

        /* Recompute τ values for swapped assignment */
        int tau_kp_swap = tau(init->kq);
        int tau_kq_swap = tau(init->kp);

        /* Reset counters but keep cumulative stats */
        long long prev_explored = result->branches_explored;
        long long prev_pruned = result->branches_pruned;
        long long prev_known = result->branches_known_pruned;
        result->branches_explored = 0;
        result->branches_pruned = 0;
        result->branches_known_pruned = 0;
        result->max_depth_reached = 0;

        ret = hs_run_once(dk, init->k, init->kq, init->kp,
                          init->tau_k, tau_kp_swap, tau_kq_swap,
                          result, verbose);

        /* Accumulate stats */
        result->branches_explored += prev_explored;
        result->branches_pruned += prev_pruned;
        result->branches_known_pruned += prev_known;

        if (ret == 0) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                                    + (end.tv_nsec - start.tv_nsec) / 1e9;
            printf("\n[HS] SUCCESS on attempt 2 (swapped kp/kq)!\n");
            hs_result_print(result);
            return 0;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                            + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("\n[HS] FAILED — key not recovered.\n");
    printf("  Total branches explored: %lld\n", result->branches_explored);
    printf("  Time: %.3f s\n", result->elapsed_seconds);
    printf("  (This can happen with too much decay — δ must be ≥ ~0.27)\n");

    return -1;
}
