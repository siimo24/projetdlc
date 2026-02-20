#define _POSIX_C_SOURCE 199309L
#include "heninger_shacham.h"
#include <time.h>

/* ========================================================================
 *  heninger_shacham.c — Branch & Prune RSA Key Reconstruction
 *
 *  ITERATIVE implementation using explicit candidate lists.
 *  (The recursive DFS approach caused stack overflow when kp/kq are wrong,
 *  because wrong kp/kq breaks equations (10)-(11) pruning, causing the
 *  search tree to explode.)
 *
 *  From HS09 paper, equations (8)-(11):
 *
 *  At bit position i, with current partial values p0, q0, d0, dp0, dq0,
 *  and candidate bits p_i, q_i, d_i, dp_i, dq_i:
 *
 *  (8)  bit_i(N - p0*q0)                          = p_i XOR q_i
 *  (9)  bit_{i+t(k)}(k(N+1) - k(p0+q0) - e*d0 + 1) = d_i XOR p_i XOR q_i
 *  (10) bit_{i+t(kp)}(kp(p0-1) + 1 - e*dp0)      = dp_i XOR p_i
 *  (11) bit_{i+t(kq)}(kq(q0-1) + 1 - e*dq0)      = dq_i XOR q_i
 *
 *  These 4 equations on 5 unknowns yield exactly 2 valid candidates
 *  per slice. Known bits further reduce this to <= 2, often to 1.
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

/* ---- Candidate management ---- */

typedef struct {
    mpz_t p, q, d, dp, dq;
} hs_candidate_t;

typedef struct {
    hs_candidate_t *items;
    int count;
    int capacity;
} hs_clist_t;

static void cand_init(hs_candidate_t *c) {
    mpz_inits(c->p, c->q, c->d, c->dp, c->dq, NULL);
}

static void cand_clear(hs_candidate_t *c) {
    mpz_clears(c->p, c->q, c->d, c->dp, c->dq, NULL);
}

static void cand_copy(hs_candidate_t *dst, const hs_candidate_t *src) {
    mpz_set(dst->p,  src->p);
    mpz_set(dst->q,  src->q);
    mpz_set(dst->d,  src->d);
    mpz_set(dst->dp, src->dp);
    mpz_set(dst->dq, src->dq);
}

static void clist_init(hs_clist_t *cl, int cap) {
    cl->capacity = cap;
    cl->count = 0;
    cl->items = (hs_candidate_t *)malloc(sizeof(hs_candidate_t) * cl->capacity);
    for (int i = 0; i < cl->capacity; i++)
        cand_init(&cl->items[i]);
}

static void clist_clear(hs_clist_t *cl) {
    for (int i = 0; i < cl->capacity; i++)
        cand_clear(&cl->items[i]);
    free(cl->items);
    cl->items = NULL;
    cl->count = 0;
    cl->capacity = 0;
}

static void clist_ensure(hs_clist_t *cl, int needed) {
    if (needed <= cl->capacity) return;
    int new_cap = cl->capacity * 2;
    if (new_cap < needed) new_cap = needed;
    cl->items = (hs_candidate_t *)realloc(cl->items,
                                           sizeof(hs_candidate_t) * new_cap);
    for (int i = cl->capacity; i < new_cap; i++)
        cand_init(&cl->items[i]);
    cl->capacity = new_cap;
}

static void clist_push(hs_clist_t *cl, const hs_candidate_t *c) {
    clist_ensure(cl, cl->count + 1);
    cand_copy(&cl->items[cl->count], c);
    cl->count++;
}

static void clist_reset(hs_clist_t *cl) {
    cl->count = 0;
}

/* ---- Equation checks ---- */

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

static bool check_eq9(const mpz_t N, const mpz_t e, const mpz_t k, int tau_k,
                      const mpz_t p0, const mpz_t q0, const mpz_t d0,
                      char p_i, char q_i, char d_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_add_ui(temp, N, 1);
    mpz_mul(term, k, temp);

    mpz_add(temp, p0, q0);
    mpz_mul(temp, k, temp);
    mpz_sub(term, term, temp);

    mpz_mul(temp, e, d0);
    mpz_sub(term, term, temp);
    mpz_add_ui(term, term, 1);

    bool ok = (get_bit(term, i + tau_k) == (d_i ^ p_i ^ q_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

static bool check_eq10(const mpz_t e, const mpz_t kp, int tau_kp,
                       const mpz_t p0, const mpz_t dp0,
                       char p_i, char dp_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_sub_ui(temp, p0, 1);
    mpz_mul(term, kp, temp);
    mpz_add_ui(term, term, 1);
    mpz_mul(temp, e, dp0);
    mpz_sub(term, term, temp);

    bool ok = (get_bit(term, i + tau_kp) == (dp_i ^ p_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

static bool check_eq11(const mpz_t e, const mpz_t kq, int tau_kq,
                       const mpz_t q0, const mpz_t dq0,
                       char q_i, char dq_i, int i)
{
    mpz_t term, temp;
    mpz_inits(term, temp, NULL);

    mpz_sub_ui(temp, q0, 1);
    mpz_mul(term, kq, temp);
    mpz_add_ui(term, term, 1);
    mpz_mul(temp, e, dq0);
    mpz_sub(term, term, temp);

    bool ok = (get_bit(term, i + tau_kq) == (dq_i ^ q_i));
    mpz_clears(term, temp, NULL);
    return ok;
}

/* ---- Termination check ---- */

static bool try_termination(const mpz_t e, const mpz_t kp, const mpz_t kq,
                            const mpz_t N, const mpz_t dp0, const mpz_t dq0,
                            mpz_t out_p, mpz_t out_q)
{
    mpz_t p_cand, q_cand, product;
    mpz_inits(p_cand, q_cand, product, NULL);

    mpz_mul(p_cand, e, dp0);
    mpz_sub_ui(p_cand, p_cand, 1);
    if (!mpz_divisible_p(p_cand, kp)) {
        mpz_clears(p_cand, q_cand, product, NULL);
        return false;
    }
    mpz_divexact(p_cand, p_cand, kp);
    mpz_add_ui(p_cand, p_cand, 1);

    mpz_mul(q_cand, e, dq0);
    mpz_sub_ui(q_cand, q_cand, 1);
    if (!mpz_divisible_p(q_cand, kq)) {
        mpz_clears(p_cand, q_cand, product, NULL);
        return false;
    }
    mpz_divexact(q_cand, q_cand, kq);
    mpz_add_ui(q_cand, q_cand, 1);

    mpz_mul(product, p_cand, q_cand);
    bool found = (mpz_cmp(product, N) == 0);

    if (found) {
        mpz_set(out_p, p_cand);
        mpz_set(out_q, q_cand);
    }

    mpz_clears(p_cand, q_cand, product, NULL);
    return found;
}

/* ========================================================================
 *  Iterative Branch & Prune
 *
 *  BFS with two candidate lists (current / next).
 *  At each bit position:
 *    1. Check termination on all current candidates
 *    2. For each candidate, test 32 slices through equations (8)-(11)
 *    3. Prune with known bits
 *    4. Add survivors to next list
 *    5. Swap current <-> next
 *
 *  Safety: if candidate count exceeds MAX_CANDIDATES, the attempt is
 *  aborted (wrong kp/kq assignment causes runaway growth).
 * ======================================================================== */

#define HS_MAX_CANDIDATES 100000

static int hs_run_once(const degraded_key_t *dk,
                       const mpz_t k, const mpz_t kp, const mpz_t kq,
                       int tau_k, int tau_kp, int tau_kq,
                       hs_result_t *result,
                       bool verbose)
{
    int half_bits = dk->half_bits;

    /* Build initial candidate */
    hs_candidate_t seed;
    cand_init(&seed);
    mpz_set_ui(seed.p, 0);  mpz_setbit(seed.p, 0);
    mpz_set_ui(seed.q, 0);  mpz_setbit(seed.q, 0);

    mpz_t modulus;
    mpz_init(modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_k + 2));
    mpz_invert(seed.d, dk->key.e, modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kp + 1));
    mpz_invert(seed.dp, dk->key.e, modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kq + 1));
    mpz_invert(seed.dq, dk->key.e, modulus);
    mpz_clear(modulus);

    if (verbose) {
        printf("[HS] Starting iterative branch & prune (half_bits=%d)\n", half_bits);
        printf("[HS] Slice(0): p[0]=%d, q[0]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n",
               get_bit(seed.p, 0), get_bit(seed.q, 0),
               tau_k, get_bit(seed.d, tau_k),
               tau_kp, get_bit(seed.dp, tau_kp),
               tau_kq, get_bit(seed.dq, tau_kq));
    }

    hs_clist_t current, next;
    clist_init(&current, 8);
    clist_init(&next, 16);
    clist_push(&current, &seed);
    cand_clear(&seed);

    const int *kn_p  = dk->known_p;
    const int *kn_q  = dk->known_q;
    const int *kn_d  = dk->known_d;
    const int *kn_dp = dk->known_dp;
    const int *kn_dq = dk->known_dq;

    for (int depth = 1; depth < half_bits; depth++) {

        if (depth > result->max_depth_reached)
            result->max_depth_reached = depth;

        if (verbose && (depth % 100 == 0 || depth == 1)) {
            printf("  [HS] Depth %d/%d, candidates: %d, explored: %lld\n",
                   depth, half_bits, current.count, result->branches_explored);
        }

        /* Check termination */
        mpz_t cp, cq;
        mpz_inits(cp, cq, NULL);
        for (int c = 0; c < current.count; c++) {
            if (try_termination(dk->key.e, kp, kq, dk->key.N,
                                current.items[c].dp, current.items[c].dq,
                                cp, cq))
            {
                result->success = true;
                mpz_set(result->recovered_key.p, cp);
                mpz_set(result->recovered_key.q, cq);
                mpz_set(result->recovered_key.N, dk->key.N);
                mpz_set(result->recovered_key.e, dk->key.e);

                mpz_t p1, q1, phi;
                mpz_inits(p1, q1, phi, NULL);
                mpz_sub_ui(p1, cp, 1);
                mpz_sub_ui(q1, cq, 1);
                mpz_mul(phi, p1, q1);
                mpz_invert(result->recovered_key.d, dk->key.e, phi);
                mpz_mod(result->recovered_key.dp, result->recovered_key.d, p1);
                mpz_mod(result->recovered_key.dq, result->recovered_key.d, q1);
                mpz_invert(result->recovered_key.qinv, cq, cp);
                result->recovered_key.bits = half_bits * 2;
                mpz_clears(p1, q1, phi, NULL);

                mpz_clears(cp, cq, NULL);
                clist_clear(&current);
                clist_clear(&next);
                return 0;
            }
        }
        mpz_clears(cp, cq, NULL);

        /* Expand candidates */
        clist_reset(&next);

        for (int c = 0; c < current.count; c++) {
            hs_candidate_t *cand = &current.items[c];

            for (int s = 0; s < 32; s++) {
                char p_i  = POSSIBILITIES[s][0];
                char q_i  = POSSIBILITIES[s][1];
                char d_i  = POSSIBILITIES[s][2];
                char dp_i = POSSIBILITIES[s][3];
                char dq_i = POSSIBILITIES[s][4];

                result->branches_explored++;

                if (!check_eq8(dk->key.N, cand->p, cand->q, p_i, q_i, depth)) {
                    result->branches_pruned++;
                    continue;
                }
                if (!check_eq9(dk->key.N, dk->key.e, k, tau_k,
                               cand->p, cand->q, cand->d,
                               p_i, q_i, d_i, depth)) {
                    result->branches_pruned++;
                    continue;
                }
                if (!check_eq10(dk->key.e, kp, tau_kp,
                                cand->p, cand->dp, p_i, dp_i, depth)) {
                    result->branches_pruned++;
                    continue;
                }
                if (!check_eq11(dk->key.e, kq, tau_kq,
                                cand->q, cand->dq, q_i, dq_i, depth)) {
                    result->branches_pruned++;
                    continue;
                }

                /* Check known bits */
                bool ok = true;

                if (kn_p[depth] != BIT_UNKNOWN && kn_p[depth] != p_i)
                    ok = false;

                if (ok && kn_q[depth] != BIT_UNKNOWN && kn_q[depth] != q_i)
                    ok = false;

                if (ok && (depth + tau_k) < half_bits * 2 &&
                    kn_d[depth + tau_k] != BIT_UNKNOWN &&
                    kn_d[depth + tau_k] != d_i)
                    ok = false;

                if (ok && (depth + tau_kp) < half_bits &&
                    kn_dp[depth + tau_kp] != BIT_UNKNOWN &&
                    kn_dp[depth + tau_kp] != dp_i)
                    ok = false;

                if (ok && (depth + tau_kq) < half_bits &&
                    kn_dq[depth + tau_kq] != BIT_UNKNOWN &&
                    kn_dq[depth + tau_kq] != dq_i)
                    ok = false;

                if (!ok) {
                    result->branches_known_pruned++;
                    continue;
                }

                /* Valid - create child */
                hs_candidate_t child;
                cand_init(&child);
                cand_copy(&child, cand);
                if (p_i)  mpz_setbit(child.p,  depth);
                if (q_i)  mpz_setbit(child.q,  depth);
                if (d_i)  mpz_setbit(child.d,  depth + tau_k);
                if (dp_i) mpz_setbit(child.dp, depth + tau_kp);
                if (dq_i) mpz_setbit(child.dq, depth + tau_kq);
                clist_push(&next, &child);
                cand_clear(&child);
            }
        }

        /* Safety: abort if tree explodes (wrong kp/kq) */
        if (next.count > HS_MAX_CANDIDATES) {
            if (verbose) {
                printf("  [HS] Candidates exploded to %d at depth %d - "
                       "aborting (likely wrong kp/kq)\n",
                       next.count, depth);
            }
            clist_clear(&current);
            clist_clear(&next);
            return -1;
        }

        if (next.count == 0) {
            if (verbose)
                printf("  [HS] No candidates at depth %d - aborting\n", depth);
            clist_clear(&current);
            clist_clear(&next);
            return -1;
        }

        /* Swap */
        hs_clist_t tmp = current;
        current = next;
        next = tmp;
        clist_reset(&next);
    }

    /* Final termination check */
    mpz_t fp, fq;
    mpz_inits(fp, fq, NULL);
    for (int c = 0; c < current.count; c++) {
        if (try_termination(dk->key.e, kp, kq, dk->key.N,
                            current.items[c].dp, current.items[c].dq,
                            fp, fq))
        {
            result->success = true;
            mpz_set(result->recovered_key.p, fp);
            mpz_set(result->recovered_key.q, fq);
            mpz_set(result->recovered_key.N, dk->key.N);
            mpz_set(result->recovered_key.e, dk->key.e);

            mpz_t p1, q1, phi;
            mpz_inits(p1, q1, phi, NULL);
            mpz_sub_ui(p1, fp, 1);
            mpz_sub_ui(q1, fq, 1);
            mpz_mul(phi, p1, q1);
            mpz_invert(result->recovered_key.d, dk->key.e, phi);
            mpz_mod(result->recovered_key.dp, result->recovered_key.d, p1);
            mpz_mod(result->recovered_key.dq, result->recovered_key.d, q1);
            mpz_invert(result->recovered_key.qinv, fq, fp);
            result->recovered_key.bits = half_bits * 2;
            mpz_clears(p1, q1, phi, NULL);
            break;
        }
    }
    mpz_clears(fp, fq, NULL);

    clist_clear(&current);
    clist_clear(&next);
    return result->success ? 0 : -1;
}

/* ---- Public entry point ---- */

int run_heninger_shacham(const degraded_key_t *dk,
                         const init_result_t *init,
                         hs_result_t *result,
                         bool try_swap,
                         bool verbose)
{
    printf("=== Heninger-Shacham Branch & Prune (HS09) ===\n\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    /* Attempt 1: kp, kq as given */
    gmp_printf("[HS] Attempt 1: kp=%Zd, kq=%Zd\n", init->kp, init->kq);
    printf("[HS] tau(k)=%d, tau(kp)=%d, tau(kq)=%d\n\n",
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

    /* Attempt 2: swap kp and kq */
    if (try_swap) {
        printf("\n[HS] Attempt 1 failed (%lld branches explored).\n",
               result->branches_explored);
        gmp_printf("[HS] Attempt 2: swapping kp=%Zd, kq=%Zd\n",
                   init->kq, init->kp);

        int tau_kp_swap = tau(init->kq);
        int tau_kq_swap = tau(init->kp);

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

    printf("\n[HS] FAILED - key not recovered.\n");
    printf("  Total branches explored: %lld\n", result->branches_explored);
    printf("  Time: %.3f s\n", result->elapsed_seconds);

    return -1;
}