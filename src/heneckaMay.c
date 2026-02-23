#include "henecka_may.h"
#include <time.h>

// ==========================================================
// Gestion de la liste des candidats
// ==========================================================

// structure d'un candidat partiel
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

// liste dynamique facon std::vector en C
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

// si on a pu de place, on double la capa (opti classique)
static void list_ensure(candidate_list_t *l, int needed) {
    if (needed <= l->capacity) return;
    int new_cap = l->capacity;
    while (new_cap < needed) new_cap *= 2;
    l->items = (hmm_candidate_t *)realloc(l->items, new_cap * sizeof(hmm_candidate_t));
    for (int i = l->capacity; i < new_cap; i++)
        candidate_init(&l->items[i]);
    l->capacity = new_cap;
}

static void list_push(candidate_list_t *l, const hmm_candidate_t *c) {
    list_ensure(l, l->count + 1);
    candidate_copy(&l->items[l->count], c);
    l->count++;
}

static void list_reset(candidate_list_t *l) {
    l->count = 0; // effacement "logique" ultra rapide, on garde la RAM allouée
}

// la table ultime des 32 combis de bits
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

// ==========================================================
// Filtres mathématiques (les fameuses equations 2-adiques)
// ==========================================================

static bool doesNEqualpq(const mpz_t N, const mpz_t p0, const mpz_t q0,
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

static bool eTimesd(const mpz_t N, const mpz_t e, const mpz_t k, int tau_k,
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

static bool eTimesdp(const mpz_t e, const mpz_t kp, int tau_kp,
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

static bool eTimesdq(const mpz_t e, const mpz_t kq, int tau_kq,
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

// ==========================================================
// Core Algorithm
// ==========================================================

// test de victoire: on derive p et q direct avec l'algebre pour voir si ca passe
static bool verify_success(const mpz_t e, const mpz_t kp, const mpz_t kq,
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

void hmm_params_default(hmm_params_t *params) {
    params->block_size_t = 0;
    params->threshold_C = 0;
    params->max_candidates = 100000;
    params->target_delta = 0.20;
}

// calcul auto des params selon l'equa de Hoeffding dans HMM10
void hmm_params_auto(hmm_params_t *params, int key_bits, double error_rate) {
    (void)key_bits;

    double delta = error_rate;
    if (delta <= 0.0) delta = 0.01;
    if (delta >= 0.5) delta = 0.49;

    double epsilon = 0.5 - delta;

    // fix taille de bloc
    if (params->block_size_t <= 0) {
        double t_float = 10.0 * log(2.0) / (epsilon * epsilon);
        params->block_size_t = (int)ceil(t_float);

        if (params->block_size_t < 5)  params->block_size_t = 5;
        if (params->block_size_t > 15) params->block_size_t = 15;
    }

    int t = params->block_size_t;

    // fix seuil C de tolérence
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

    int t = params->block_size_t;
    printf("  Expected correct:  ~%.1f / %d matches\n",
           5.0 * t * (1.0 - params->target_delta), 5 * t);
    printf("  Expected wrong:    ~%.1f / %d matches\n",
           5.0 * t * 0.5, 5 * t);
    printf("======================\n");
}


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


// context global pdt la fouille de l'arbre
typedef struct {
    const mpz_t *N;
    const mpz_t *e;
    const mpz_t *k;
    const mpz_t *kp;
    const mpz_t *kq;
    int tau_k, tau_kp, tau_kq;
    int half_bits;

    const mpz_t *deg_p;
    const mpz_t *deg_q;
    const mpz_t *deg_d;
    const mpz_t *deg_dp;
    const mpz_t *deg_dq;

    int block_start;
    int block_size;
    int threshold;

    candidate_list_t *output;
    int max_candidates;
    long long *candidates_generated;
} hmm_expand_ctx_t;

// fct recusive qui fouille toutes les combis possible pour un bloc donné
static void process_block_of_bits(
    hmm_expand_ctx_t *ctx,
    mpz_t cur_p, mpz_t cur_q, mpz_t cur_d, mpz_t cur_dp, mpz_t cur_dq,
    int depth_in_block,
    int match_count)
{
    // crash prev : si l'arbre explose (mauvais kp/kq), on tej pour pas freeze le CPU
    if (*ctx->candidates_generated > ctx->max_candidates * 10) {
        return;
    }

    int global_pos = ctx->block_start + depth_in_block;

    // fin du bloc t : c'est l'heure du couperet, on check le seuil C
    if (depth_in_block == ctx->block_size) {
        (*ctx->candidates_generated)++;

        if (match_count >= ctx->threshold) {
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

    // si on deborde de la moitie de N, on stop
    if (global_pos >= ctx->half_bits) {
        (*ctx->candidates_generated)++;
        if (match_count >= ctx->threshold || depth_in_block < ctx->block_size / 2) {
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

    // elagage prematuré : meme si TOUS les prochains bits sont parfaits, on atteindra pas C
    int remaining_bits = ctx->block_size - depth_in_block;
    int max_possible = match_count + 5 * remaining_bits;
    if (max_possible < ctx->threshold) {
        return; // rip le candidat
    }

    // OPTI MAJEURE : On declare et init les vars ici EN DEHORS du for.
    // Avant, elles etaient dans la boucle, ca faisait 2.5 millions de malloc/free par bloc !
    mpz_t cl_p, cl_q, cl_d, cl_dp, cl_dq;
    mpz_inits(cl_p, cl_q, cl_d, cl_dp, cl_dq, NULL);

    for (int s = 0; s < 32; s++) {
        char p_i  = POSSIBILITIES[s][0];
        char q_i  = POSSIBILITIES[s][1];
        char d_i  = POSSIBILITIES[s][2];
        char dp_i = POSSIBILITIES[s][3];
        char dq_i = POSSIBILITIES[s][4];

        // filtres mathématiques violents
        if (!doesNEqualpq(*ctx->N, cur_p, cur_q, p_i, q_i, global_pos)) continue;
        if (!eTimesd(*ctx->N, *ctx->e, *ctx->k, ctx->tau_k, cur_p, cur_q, cur_d, p_i, q_i, d_i, global_pos)) continue;
        if (!eTimesdp(*ctx->e, *ctx->kp, ctx->tau_kp, cur_p, cur_dp, p_i, dp_i, global_pos)) continue;
        if (!eTimesdq(*ctx->e, *ctx->kq, ctx->tau_kq, cur_q, cur_dq, q_i, dq_i, global_pos)) continue;

        // le chemin est legit niveau maths, on le note contre la RAM petee (erreur tolerée)
        int new_matches = match_count;

        if (p_i == get_bit(*ctx->deg_p, global_pos)) new_matches++;
        if (q_i == get_bit(*ctx->deg_q, global_pos)) new_matches++;

        int d_pos = global_pos + ctx->tau_k;
        if (d_pos < ctx->half_bits * 2) {
            if (d_i == get_bit(*ctx->deg_d, d_pos)) new_matches++;
        }

        int dp_pos = global_pos + ctx->tau_kp;
        if (dp_pos < ctx->half_bits) {
            if (dp_i == get_bit(*ctx->deg_dp, dp_pos)) new_matches++;
        }

        int dq_pos = global_pos + ctx->tau_kq;
        if (dq_pos < ctx->half_bits) {
            if (dq_i == get_bit(*ctx->deg_dq, dq_pos)) new_matches++;
        }

        // MAJ de l'etat des vars locales, on ecrase les valeurs par celles du parent
        mpz_set(cl_p,  cur_p);
        mpz_set(cl_q,  cur_q);
        mpz_set(cl_d,  cur_d);
        mpz_set(cl_dp, cur_dp);
        mpz_set(cl_dq, cur_dq);

        // FIX "Ghost Bit" : On set à 1, mais on force AUSSI à 0.
        // Sans le clrbit, des bits résiduels pouvaient corrompre les calculs 2-adiques
        if (p_i)  mpz_setbit(cl_p,  global_pos); else mpz_clrbit(cl_p, global_pos);
        if (q_i)  mpz_setbit(cl_q,  global_pos); else mpz_clrbit(cl_q, global_pos);
        if (d_i)  mpz_setbit(cl_d,  global_pos + ctx->tau_k); else mpz_clrbit(cl_d, global_pos + ctx->tau_k);
        if (dp_i) mpz_setbit(cl_dp, global_pos + ctx->tau_kp); else mpz_clrbit(cl_dp, global_pos + ctx->tau_kp);
        if (dq_i) mpz_setbit(cl_dq, global_pos + ctx->tau_kq); else mpz_clrbit(cl_dq, global_pos + ctx->tau_kq);

        // On s'enfonce plus profond
        process_block_of_bits(ctx, cl_p, cl_q, cl_d, cl_dp, cl_dq,
                               depth_in_block + 1, new_matches);

        // check de securite : on break au lieu de return pour bien free la RAM a la fin
        if (ctx->output->count >= ctx->max_candidates ||
            *ctx->candidates_generated > ctx->max_candidates * 10) {
            break;
        }
    }

    // menage, on nettoie TOUT le block d'un coup
    mpz_clears(cl_p, cl_q, cl_d, cl_dp, cl_dq, NULL);
}

// --- Moteur principal HMM10 ---

static int main_hmm_function(const degraded_key_t *dk,
                        const mpz_t k, const mpz_t kp, const mpz_t kq,
                        int tau_k, int tau_kp, int tau_kq,
                        const hmm_params_t *params,
                        hmm_result_t *result,
                        bool verbose)
{
    int half_bits = dk->half_bits;
    int t = params->block_size_t;
    int C = params->threshold_C;
    int num_blocks = (half_bits + t - 1) / t;

    candidate_list_t current, next;
    list_init(&current, 16);
    list_init(&next, 64);

    // root: p et q sont tjs impairs
    hmm_candidate_t start;
    candidate_init(&start);
    mpz_set_ui(start.p, 0);
    mpz_set_ui(start.q, 0);
    mpz_setbit(start.p, 0);
    mpz_setbit(start.q, 0);

    // correction lsb
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

    // on decoupe la clef en tranche et on avance par blocs
    for (int block = 0; block < num_blocks; block++) {
        int block_start = 1 + block * t;

        if (block_start >= half_bits) break;

        list_reset(&next);

        // check precoce de victoire pour tous les mecs en lice
        mpz_t cand_p, cand_q;
        mpz_inits(cand_p, cand_q, NULL);

        for (int c = 0; c < current.count; c++) {
            if (verify_success(dk->key.e, kp, kq, dk->key.N,
                                current.items[c].dp, current.items[c].dq,
                                cand_p, cand_q))
            {
                // BINGO !
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

                if (!verbose) printf("\n"); // clean l'affichage
                return 0;
            }
        }
        mpz_clears(cand_p, cand_q, NULL);

        long long block_generated = 0;

        hmm_expand_ctx_t ctx;
        ctx.N = &dk->key.N;
        ctx.e = &dk->key.e;

        // fix crade mais obligatoire pour passer les mpz_t sans segfault
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

        // petit ajustement du seuil si on coupe court sur le dernier bloc
        int effective_t = t;
        if (block_start + t > half_bits)
            effective_t = half_bits - block_start;
        ctx.threshold = (effective_t < t) ?
            (int)floor((double)C * effective_t / t) : C;
        ctx.output = &next;
        ctx.max_candidates = params->max_candidates;
        ctx.candidates_generated = &block_generated;

        for (int c = 0; c < current.count; c++) {
            process_block_of_bits(
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

        // feedback visuel pour savoir que ca avance
        if (verbose) {
            printf("  [HMM] Block %d/%d (bits %d-%d): %d candidates → %d survivors"
                   " (%lld generated)\n",
                   block + 1, num_blocks, block_start,
                   block_start + effective_t - 1,
                   current.count, next.count, block_generated);
        } else {
            // ptite anim de points pour pas faire croire au user que c'est freezé
            printf(".");
            fflush(stdout);
        }

        if (next.count == 0) {
            if (verbose) printf("  [HMM] Game over, tous les candidats sont morts.\n");
            else printf("\n");
            list_clear(&current);
            list_clear(&next);
            return -1;
        }

        candidate_list_t tmp = current;
        current = next;
        next = tmp;
        list_reset(&next);
    }

    // si ca tient jusque la, ultime check sur les survivants
    mpz_t final_p, final_q;
    mpz_inits(final_p, final_q, NULL);

    for (int c = 0; c < current.count; c++) {
        if (verify_success(dk->key.e, kp, kq, dk->key.N,
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

    if (!verbose) printf("\n");
    return result->success ? 0 : -1;
}

// ==========================================================
// Wrapper public (appelé par main.c)
// ==========================================================

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

    gmp_printf("[HMM] Attempt 1: kp=%Zd, kq=%Zd\n", init->kp, init->kq);
    printf("[HMM] τ(k)=%d, τ(kp)=%d, τ(kq)=%d\n", init->tau_k, init->tau_kp, init->tau_kq);

    int ret = main_hmm_function(dk, init->k, init->kp, init->kq,
                           init->tau_k, init->tau_kp, init->tau_kq,
                           params, result, verbose);

    if (ret == 0) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        result->elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("[HMM] SUCCESS on attempt 1!\n");
        hmm_result_print(result);
        return 0;
    }

    // si ca a fail, on part du principe que kp et kq ont ete lock sur le mauvais prime
    if (try_swap) {
        printf("[HMM] Attempt 1 failed. Retrying with swapped kp/kq ...\n");

        int tau_kp_swap = tau(init->kq);
        int tau_kq_swap = tau(init->kp);

        long long prev_gen = result->total_candidates_generated;
        long long prev_pruned = result->total_candidates_pruned;
        int prev_max = result->max_candidates_at_once;
        result->total_candidates_generated = 0;
        result->total_candidates_pruned = 0;
        result->max_candidates_at_once = 0;
        result->blocks_processed = 0;

        ret = main_hmm_function(dk, init->k, init->kq, init->kp,
                           init->tau_k, tau_kp_swap, tau_kq_swap,
                           params, result, verbose);

        result->total_candidates_generated += prev_gen;
        result->total_candidates_pruned += prev_pruned;
        if (prev_max > result->max_candidates_at_once)
            result->max_candidates_at_once = prev_max;

        if (ret == 0) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            result->elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
            printf("[HMM] SUCCESS on attempt 2 (swapped kp/kq)!\n");
            hmm_result_print(result);
            return 0;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    result->elapsed_seconds = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("[HMM] FAILED — key not recovered.\n");
    printf("  Total candidates generated: %lld\n", result->total_candidates_generated);
    printf("  Time: %.3f s\n", result->elapsed_seconds);
    printf("  (Try reducing δ or adjusting block_size/threshold)\n");

    return -1;
}