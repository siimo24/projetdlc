#include "henecka_may.h"
#include <time.h>


// struct d'un candidat donné

typedef struct {
    mpz_t p, q, d, dp, dq;
} hmm_candidate_t;

static void candidate_init(hmm_candidate_t *c) {
    mpz_inits(c->p, c->q, c->d, c->dp, c->dq, NULL);
}

static void candidate_clear(hmm_candidate_t *c) {
    mpz_clears(c->p, c->q, c->d, c->dp, c->dq, NULL);
}

static void candidate_copy(hmm_candidate_t *dst, const hmm_candidate_t *src) {
    mpz_set(dst->p,  src->p);
    mpz_set(dst->q,  src->q);
    mpz_set(dst->d,  src->d);
    mpz_set(dst->dp, src->dp);
    mpz_set(dst->dq, src->dq);
}

// liste dynamique de candidats + fct utilitaires
typedef struct {
    hmm_candidate_t *items;
    int count;
    int capacity;
} candidate_list_t;

static void list_init(candidate_list_t *l, int initial_cap) {
    l->capacity = initial_cap > 0 ? initial_cap : 16;
    l->count = 0;
    l->items = (hmm_candidate_t *)malloc(l->capacity * sizeof(hmm_candidate_t));
    for (int i = 0; i < l->capacity; i++)
        candidate_init(&l->items[i]);
}

static void list_clear(candidate_list_t *l) {
    for (int i = 0; i < l->capacity; i++)
        candidate_clear(&l->items[i]);
    free(l->items);
    l->items = NULL;
    l->count = 0;
    l->capacity = 0;
}

// si la liste est pleine on lui  double la taille avec realloc
// c'est équivalent au std vector en cpp
static void list_ensure(candidate_list_t *l, int needed) {
    if (needed <= l->capacity) return;
    int new_cap = l->capacity;
    while (new_cap < needed) new_cap *= 2;
    l->items = (hmm_candidate_t *)realloc(l->items, new_cap * sizeof(hmm_candidate_t));
    for (int i = l->capacity; i < new_cap; i++)
        candidate_init(&l->items[i]);
    l->capacity = new_cap;
}

// ajout d'un candidat
static void list_push(candidate_list_t *l, const hmm_candidate_t *c) {
    list_ensure(l, l->count + 1);
    candidate_copy(&l->items[l->count], c);
    l->count++;
}

static void list_reset(candidate_list_t *l) {
    l->count = 0;
}

// liste statique qui permettra des comparaisons rapides

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

/* ---- Equation checks (same as HS09) ---- */

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

/* ---- Parameters ---- */

void hmm_params_default(hmm_params_t *params) {
    params->block_size_t = 0;       /* Auto */
    params->threshold_C = 0;        /* Auto */
    params->max_candidates = 100000;
    params->target_delta = 0.20;
}

/*
 * Auto-compute parameters from the paper:
 *   t ≈ ⌈10 · ln(2) / ε²⌉  where ε = 1/2 − δ   (simplified from paper)
 *   γ₀ = √((1 + 1/t) · ln(2) / 10)
 *   C = ⌊5t · (1/2 + γ₀)⌋
 *
 * The factor 5 comes from 5 components per slice.
 */
void hmm_params_auto(hmm_params_t *params, int key_bits, double error_rate) {
    (void)key_bits;  /* Not needed for current formula, but kept for flexibility */

    double delta = error_rate;
    if (delta <= 0.0) delta = 0.01;
    if (delta >= 0.5) delta = 0.49;

    double epsilon = 0.5 - delta;

    /* Block size t */
    if (params->block_size_t <= 0) {
        double t_float = 10.0 * log(2.0) / (epsilon * epsilon);
        params->block_size_t = (int)ceil(t_float);

        /* Clamp to practical range */
        if (params->block_size_t < 5)  params->block_size_t = 5;
        if (params->block_size_t > 15) params->block_size_t = 15;
    }

    int t = params->block_size_t;

    /* Threshold C */
    if (params->threshold_C <= 0) {
        double gamma_0 = sqrt((1.0 + 1.0 / t) * log(2.0) / 10.0);
        params->threshold_C = (int)floor(5.0 * t * (0.5 + gamma_0));
    }

    params->target_delta = delta;
}

void hmm_params_print(const hmm_params_t *params) {
    printf("=== HMM Parameters ===\n");
    printf("  Block size (t):    %d\n", params->block_size_t);
    printf("  Threshold (C):     %d  (max per block = %d)\n",
           params->threshold_C, 5 * params->block_size_t);
    printf("  Max candidates:    %d\n", params->max_candidates);
    printf("  Target δ:          %.4f\n", params->target_delta);

    /* Show expected match counts */
    int t = params->block_size_t;
    printf("  Expected correct:  ~%.1f / %d matches\n",
           5.0 * t * (1.0 - params->target_delta), 5 * t);
    printf("  Expected wrong:    ~%.1f / %d matches\n",
           5.0 * t * 0.5, 5 * t);
    printf("======================\n");
}

/* ---- Result init/clear ---- */

void hmm_result_init(hmm_result_t *res) {
    rsa_key_init(&res->recovered_key);
    res->success = false;
    res->total_candidates_generated = 0;
    res->total_candidates_pruned = 0;
    res->max_candidates_at_once = 0;
    res->blocks_processed = 0;
    res->elapsed_seconds = 0.0;
}

void hmm_result_clear(hmm_result_t *res) {
    rsa_key_clear(&res->recovered_key);
}

void hmm_result_print(const hmm_result_t *res) {
    printf("=== Henecka-May-Meurer Results ===\n");
    printf("  Success:              %s\n", res->success ? "YES" : "NO");
    printf("  Blocks processed:     %d\n", res->blocks_processed);
    printf("  Candidates generated: %lld\n", res->total_candidates_generated);
    printf("  Candidates pruned:    %lld\n", res->total_candidates_pruned);
    printf("  Max candidates:       %d\n", res->max_candidates_at_once);
    printf("  Time:                 %.3f s\n", res->elapsed_seconds);
    if (res->success) {
        gmp_printf("  p = %Zd\n", res->recovered_key.p);
        gmp_printf("  q = %Zd\n", res->recovered_key.q);
    }
    printf("==================================\n");
}

/* ========================================================================
 *  Block expansion: depth-first recursive expansion within a t-bit block.
 *
 *  For each bit position, we try all valid slices (from equations 8-11).
 *  For each valid slice, we count how many of its 5 bits match the
 *  degraded key bits at the corresponding positions. This match count
 *  is accumulated. At the end of the block (depth == t), we check if
 *  match_count >= C. If so, the candidate is added to the output list.
 *
 *  This avoids materializing all 2^t leaves — only survivors are stored.
 * ======================================================================== */

typedef struct {
    /* Constants */
    const mpz_t *N;
    const mpz_t *e;
    const mpz_t *k;
    const mpz_t *kp;
    const mpz_t *kq;
    int tau_k, tau_kp, tau_kq;
    int half_bits;

    /* Degraded key bit values (for matching) */
    const mpz_t *deg_p;
    const mpz_t *deg_q;
    const mpz_t *deg_d;
    const mpz_t *deg_dp;
    const mpz_t *deg_dq;

    /* Block parameters */
    int block_start;    /* Global bit index where this block starts */
    int block_size;     /* t */
    int threshold;      /* C */

    /* Output */
    candidate_list_t *output;
    int max_candidates;

    /* Stats */
    long long *candidates_generated;
} hmm_expand_ctx_t;

static void expand_block_recursive(
    hmm_expand_ctx_t *ctx,
    mpz_t cur_p, mpz_t cur_q, mpz_t cur_d, mpz_t cur_dp, mpz_t cur_dq,
    int depth_in_block,
    int match_count)
{
    /* NEW: Prevent massive CPU hang if kp/kq are swapped and the tree explodes.
     * If we generate way too many candidates in a single block, it means our
     * equations are not pruning properly (which happens when kp/kq are wrong). */
    if (*ctx->candidates_generated > ctx->max_candidates * 10) {
        return;
    }

    int global_pos = ctx->block_start + depth_in_block;

    /* ---- Block boundary: threshold check ---- */
    if (depth_in_block == ctx->block_size) {
        (*ctx->candidates_generated)++;

        if (match_count >= ctx->threshold) {
            /* Survivor — add to output if under limit */
            if (ctx->output->count < ctx->max_candidates) {
                hmm_candidate_t cand;
                candidate_init(&cand);
                mpz_set(cand.p,  cur_p);
                mpz_set(cand.q,  cur_q);
                mpz_set(cand.d,  cur_d);
                mpz_set(cand.dp, cur_dp);
                mpz_set(cand.dq, cur_dq);
                list_push(ctx->output, &cand);
                candidate_clear(&cand);
            }
        }
        return;
    }

    /* ---- Past the half-bit boundary: stop expanding ---- */
    if (global_pos >= ctx->half_bits) {
        /* Treat as end of block — check threshold with what we have */
        (*ctx->candidates_generated)++;
        if (match_count >= ctx->threshold ||
            depth_in_block < ctx->block_size / 2)  /* lenient at key end */
        {
            if (ctx->output->count < ctx->max_candidates) {
                hmm_candidate_t cand;
                candidate_init(&cand);
                mpz_set(cand.p,  cur_p);
                mpz_set(cand.q,  cur_q);
                mpz_set(cand.d,  cur_d);
                mpz_set(cand.dp, cur_dp);
                mpz_set(cand.dq, cur_dq);
                list_push(ctx->output, &cand);
                candidate_clear(&cand);
            }
        }
        return;
    }

    /* ---- Early termination: even with max remaining matches, can't reach C ---- */
    int remaining_bits = ctx->block_size - depth_in_block;
    int max_possible = match_count + 5 * remaining_bits;
    if (max_possible < ctx->threshold) {
        return;  /* Prune: impossible to reach threshold */
    }

    /* ---- Try all 32 slice candidates ---- */
    for (int s = 0; s < 32; s++) {
        char p_i  = POSSIBILITIES[s][0];
        char q_i  = POSSIBILITIES[s][1];
        char d_i  = POSSIBILITIES[s][2];
        char dp_i = POSSIBILITIES[s][3];
        char dq_i = POSSIBILITIES[s][4];

        /* Check equations (8)-(11) */
        if (!check_eq8(*ctx->N, cur_p, cur_q, p_i, q_i, global_pos))
            continue;
        if (!check_eq9(*ctx->N, *ctx->e, *ctx->k, ctx->tau_k,
                       cur_p, cur_q, cur_d, p_i, q_i, d_i, global_pos))
            continue;
        if (!check_eq10(*ctx->e, *ctx->kp, ctx->tau_kp,
                        cur_p, cur_dp, p_i, dp_i, global_pos))
            continue;
        if (!check_eq11(*ctx->e, *ctx->kq, ctx->tau_kq,
                        cur_q, cur_dq, q_i, dq_i, global_pos))
            continue;

        /* Count matches with degraded key */
        int new_matches = match_count;

        if (p_i == get_bit(*ctx->deg_p, global_pos))
            new_matches++;
        if (q_i == get_bit(*ctx->deg_q, global_pos))
            new_matches++;

        int d_pos = global_pos + ctx->tau_k;
        if (d_pos < ctx->half_bits * 2) {
            if (d_i == get_bit(*ctx->deg_d, d_pos))
                new_matches++;
        }

        int dp_pos = global_pos + ctx->tau_kp;
        if (dp_pos < ctx->half_bits) {
            if (dp_i == get_bit(*ctx->deg_dp, dp_pos))
                new_matches++;
        }

        int dq_pos = global_pos + ctx->tau_kq;
        if (dq_pos < ctx->half_bits) {
            if (dq_i == get_bit(*ctx->deg_dq, dq_pos))
                new_matches++;
        }

        /* Clone, set bits, and recurse */
        mpz_t cl_p, cl_q, cl_d, cl_dp, cl_dq;
        mpz_init_set(cl_p,  cur_p);
        mpz_init_set(cl_q,  cur_q);
        mpz_init_set(cl_d,  cur_d);
        mpz_init_set(cl_dp, cur_dp);
        mpz_init_set(cl_dq, cur_dq);

        if (p_i)  mpz_setbit(cl_p,  global_pos);
        if (q_i)  mpz_setbit(cl_q,  global_pos);
        if (d_i)  mpz_setbit(cl_d,  global_pos + ctx->tau_k);
        if (dp_i) mpz_setbit(cl_dp, global_pos + ctx->tau_kp);
        if (dq_i) mpz_setbit(cl_dq, global_pos + ctx->tau_kq);

        expand_block_recursive(ctx, cl_p, cl_q, cl_d, cl_dp, cl_dq,
                               depth_in_block + 1, new_matches);

        mpz_clears(cl_p, cl_q, cl_d, cl_dp, cl_dq, NULL);

        /* Early exit if output is full */
        /* Early exit if output is full or safety limit reached */
        if (ctx->output->count >= ctx->max_candidates ||
            *ctx->candidates_generated > ctx->max_candidates * 10) {
            return;
        }
    }
}

/* ========================================================================
 *  Main reconstruction driver
 * ======================================================================== */

static int hmm_run_once(const degraded_key_t *dk,
                        const mpz_t k, const mpz_t kp, const mpz_t kq,
                        int tau_k, int tau_kp, int tau_kq,
                        const hmm_params_t *params,
                        hmm_result_t *result,
                        bool verbose)
{
    int half_bits = dk->half_bits;
    int t = params->block_size_t;
    int C = params->threshold_C;
    int num_blocks = (half_bits + t - 1) / t;  /* ceiling division */

    /* ---- Initialize starting candidate ---- */
    candidate_list_t current, next;
    list_init(&current, 16);
    list_init(&next, 64);

    /* Initial candidate: p[0]=1, q[0]=1, corrected LSBs */
    hmm_candidate_t start;
    candidate_init(&start);
    mpz_set_ui(start.p, 0);
    mpz_set_ui(start.q, 0);
    mpz_setbit(start.p, 0);
    mpz_setbit(start.q, 0);

    /* Correct LSBs */
    mpz_t modulus;
    mpz_init(modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_k + 2));
    mpz_invert(start.d, dk->key.e, modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kp + 1));
    mpz_invert(start.dp, dk->key.e, modulus);
    mpz_ui_pow_ui(modulus, 2, (unsigned long)(tau_kq + 1));
    mpz_invert(start.dq, dk->key.e, modulus);
    mpz_clear(modulus);

    list_push(&current, &start);
    candidate_clear(&start);

    if (verbose) {
        printf("[HMM] Starting: %d blocks of %d bits (C=%d, half_bits=%d)\n",
               num_blocks, t, C, half_bits);
    }

    /* ---- Process blocks ---- */
    for (int block = 0; block < num_blocks; block++) {
        int block_start = 1 + block * t;  /* bit 0 already set */

        if (block_start >= half_bits) break;

        list_reset(&next);

        /* Check termination on each current candidate before expanding */
        mpz_t cand_p, cand_q;
        mpz_inits(cand_p, cand_q, NULL);

        for (int c = 0; c < current.count; c++) {
            if (try_termination(dk->key.e, kp, kq, dk->key.N,
                                current.items[c].dp, current.items[c].dq,
                                cand_p, cand_q))
            {
                /* Found it! */
                result->success = true;
                mpz_set(result->recovered_key.p, cand_p);
                mpz_set(result->recovered_key.q, cand_q);
                mpz_set(result->recovered_key.N, dk->key.N);
                mpz_set(result->recovered_key.e, dk->key.e);

                mpz_t p1, q1, phi;
                mpz_inits(p1, q1, phi, NULL);
                mpz_sub_ui(p1, cand_p, 1);
                mpz_sub_ui(q1, cand_q, 1);
                mpz_mul(phi, p1, q1);
                mpz_invert(result->recovered_key.d, dk->key.e, phi);
                mpz_mod(result->recovered_key.dp, result->recovered_key.d, p1);
                mpz_mod(result->recovered_key.dq, result->recovered_key.d, q1);
                mpz_invert(result->recovered_key.qinv, cand_q, cand_p);
                result->recovered_key.bits = half_bits * 2;
                mpz_clears(p1, q1, phi, NULL);

                mpz_clears(cand_p, cand_q, NULL);
                list_clear(&current);
                list_clear(&next);
                return 0;
            }
        }
        mpz_clears(cand_p, cand_q, NULL);

        /* Expand each candidate through the block */
        long long block_generated = 0;

        hmm_expand_ctx_t ctx;
        ctx.N = &dk->key.N;
        ctx.e = &dk->key.e;

        /* FIX: Cast the decayed pointers correctly to pass the heap arrays */
        ctx.k = (const mpz_t *)k;
        ctx.kp = (const mpz_t *)kp;
        ctx.kq = (const mpz_t *)kq;

        ctx.tau_k = tau_k;
        ctx.tau_kp = tau_kp;
        ctx.tau_kq = tau_kq;
        ctx.half_bits = half_bits;
        ctx.deg_p  = &dk->key.p;
        ctx.deg_q  = &dk->key.q;
        ctx.deg_d  = &dk->key.d;
        ctx.deg_dp = &dk->key.dp;
        ctx.deg_dq = &dk->key.dq;
        ctx.block_start = block_start;
        ctx.block_size = t;

        /* Adjust threshold for last block if it's shorter */
        int effective_t = t;
        if (block_start + t > half_bits)
            effective_t = half_bits - block_start;
        ctx.threshold = (effective_t < t) ?
            (int)floor((double)C * effective_t / t) : C;
        ctx.output = &next;
        ctx.max_candidates = params->max_candidates;
        ctx.candidates_generated = &block_generated;

        for (int c = 0; c < current.count; c++) {
            expand_block_recursive(
                &ctx,
                current.items[c].p, current.items[c].q,
                current.items[c].d, current.items[c].dp, current.items[c].dq,
                0, 0);
        }

        result->total_candidates_generated += block_generated;
        result->total_candidates_pruned += (block_generated - next.count);
        result->blocks_processed++;

        if (next.count > result->max_candidates_at_once)
            result->max_candidates_at_once = next.count;

        if (verbose) {
            printf("  [HMM] Block %d/%d (bits %d-%d): %d candidates → %d survivors"
                   " (%lld generated)\n",
                   block + 1, num_blocks, block_start,
                   block_start + effective_t - 1,
                   current.count, next.count, block_generated);
        }

        /* No survivors → algorithm failed */
        if (next.count == 0) {
            if (verbose)
                printf("  [HMM] No candidates survived! Aborting.\n");
            list_clear(&current);
            list_clear(&next);
            return -1;
        }

        /* Swap current and next */
        candidate_list_t tmp = current;
        current = next;
        next = tmp;
        list_reset(&next);
    }

    /* ---- Final termination check on remaining candidates ---- */
    mpz_t final_p, final_q;
    mpz_inits(final_p, final_q, NULL);

    for (int c = 0; c < current.count; c++) {
        if (try_termination(dk->key.e, kp, kq, dk->key.N,
                            current.items[c].dp, current.items[c].dq,
                            final_p, final_q))
        {
            result->success = true;
            mpz_set(result->recovered_key.p, final_p);
            mpz_set(result->recovered_key.q, final_q);
            mpz_set(result->recovered_key.N, dk->key.N);
            mpz_set(result->recovered_key.e, dk->key.e);

            mpz_t p1, q1, phi;
            mpz_inits(p1, q1, phi, NULL);
            mpz_sub_ui(p1, final_p, 1);
            mpz_sub_ui(q1, final_q, 1);
            mpz_mul(phi, p1, q1);
            mpz_invert(result->recovered_key.d, dk->key.e, phi);
            mpz_mod(result->recovered_key.dp, result->recovered_key.d, p1);
            mpz_mod(result->recovered_key.dq, result->recovered_key.d, q1);
            mpz_invert(result->recovered_key.qinv, final_q, final_p);
            result->recovered_key.bits = half_bits * 2;
            mpz_clears(p1, q1, phi, NULL);
            break;
        }
    }

    mpz_clears(final_p, final_q, NULL);
    list_clear(&current);
    list_clear(&next);

    return result->success ? 0 : -1;
}

/* ---- Public entry point ---- */

int run_henecka_may(const degraded_key_t *dk,
                    const init_result_t *init,
                    const hmm_params_t *params,
                    hmm_result_t *result,
                    bool try_swap,
                    bool verbose)
{
    printf("╔══════════════════════════════════════════════════╗\n");
    printf("║   Henecka-May-Meurer Block Threshold (HMM10)     ║\n");
    printf("╚══════════════════════════════════════════════════╝\n\n");

    hmm_params_print(params);
    printf("\n");

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    /* ---- Attempt 1: kp, kq as given ---- */
    gmp_printf("[HMM] Attempt 1: kp=%Zd, kq=%Zd\n", init->kp, init->kq);
    printf("[HMM] τ(k)=%d, τ(kp)=%d, τ(kq)=%d\n\n",
           init->tau_k, init->tau_kp, init->tau_kq);

    int ret = hmm_run_once(dk, init->k, init->kp, init->kq,
                           init->tau_k, init->tau_kp, init->tau_kq,
                           params, result, verbose);

    if (ret == 0) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                                + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("\n[HMM] SUCCESS on attempt 1!\n");
        hmm_result_print(result);
        return 0;
    }

    /* ---- Attempt 2: swap kp/kq ---- */
    if (try_swap) {
        printf("\n[HMM] Attempt 1 failed. Retrying with swapped kp/kq ...\n");

        int tau_kp_swap = tau(init->kq);
        int tau_kq_swap = tau(init->kp);

        /* Keep cumulative stats */
        long long prev_gen = result->total_candidates_generated;
        long long prev_pruned = result->total_candidates_pruned;
        int prev_max = result->max_candidates_at_once;
        result->total_candidates_generated = 0;
        result->total_candidates_pruned = 0;
        result->max_candidates_at_once = 0;
        result->blocks_processed = 0;

        ret = hmm_run_once(dk, init->k, init->kq, init->kp,
                           init->tau_k, tau_kp_swap, tau_kq_swap,
                           params, result, verbose);

        result->total_candidates_generated += prev_gen;
        result->total_candidates_pruned += prev_pruned;
        if (prev_max > result->max_candidates_at_once)
            result->max_candidates_at_once = prev_max;

        if (ret == 0) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                                    + (end.tv_nsec - start.tv_nsec) / 1e9;
            printf("\n[HMM] SUCCESS on attempt 2 (swapped kp/kq)!\n");
            hmm_result_print(result);
            return 0;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    result->elapsed_seconds = (end.tv_sec - start.tv_sec)
                            + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("\n[HMM] FAILED — key not recovered.\n");
    printf("  Total candidates generated: %lld\n", result->total_candidates_generated);
    printf("  Time: %.3f s\n", result->elapsed_seconds);
    printf("  (Try reducing δ or adjusting block_size/threshold)\n");

    return -1;
}
