#include "decay_sim.h"


static unsigned long rng_state;

static void rng_seed(unsigned long seed) {
    rng_state = seed ? seed : (unsigned long)time(NULL);
}

// entier uniforme dans 0,1
static double rng_uniform(void) {
    /* xorshift64 — fast and reasonable quality */
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state & 0x7FFFFFFFFFFFFFFF) / (double)0x7FFFFFFFFFFFFFFF;
}

//on aura besoin du bernoulli pour décider si on va dégrader un  bit ou pas
static bool rng_bernoulli(double p) {
    return rng_uniform() < p;
}

/* Modélsiation du modèle de dégradation selon le papier :
 *                    k (rate)     t_mid (seconds)
 *   Room (25°C):     0.5          10
 *   Cold (-50°C):    0.01         600
 *   Azote liquide (-196°C):   0.0003       36000
 */

void decay_params_preset(decay_params_t *params, const char *preset) {
    memset(params, 0, sizeof(*params));

    params->wrong_dir_prob = 0.001;         // delta1
    params->ground_state_block_size = 256;  // taille du bloc ayant un ground identique
    params->random_ground_start = true;
    params->seed = 0;                       // 0 donc on utilise time

    if (strcmp(preset, "room") == 0) { // on verifie en fct de ce qui a été donné dans les flags
        params->temperature_celsius = 25.0;
        params->time_seconds = 10.0;
        params->decay_rate_k = 0.5;
        params->decay_midpoint = 10.0;
    }
    else if (strcmp(preset, "cold") == 0) {
        params->temperature_celsius = -50.0;
        params->time_seconds = 60.0;
        params->decay_rate_k = 0.01;
        params->decay_midpoint = 600.0;
    }
    else if (strcmp(preset, "frozen") == 0) {
        params->temperature_celsius = -196.0;
        params->time_seconds = 3600.0;  // 60 minutes
        params->decay_rate_k = 0.0003;
        params->decay_midpoint = 36000.0;
    }
    else {
        // autre que les cas classiques, ils sont passés via les flags
        params->temperature_celsius = 25.0;
        params->time_seconds = 0.0;
        params->decay_rate_k = 0.5;
        params->decay_midpoint = 10.0;
    }
}

double compute_decay_probability(const decay_params_t *params) {
    double exponent = -params->decay_rate_k * (params->time_seconds - params->decay_midpoint);

    /* Clamp to avoid overflow */
    if (exponent > 500.0)  return 0.0;
    if (exponent < -500.0) return 1.0;

    return 1.0 / (1.0 + exp(exponent));
}


// résultats de la simulation de dégradation
void decay_params_print(const decay_params_t *params) {
    double decay_prob = compute_decay_probability(params);

    printf("=== Decay Parameters ===\n");
    printf("  Temperature:   %.1f°C\n", params->temperature_celsius);
    printf("  Time:          %.1f s\n", params->time_seconds);
    printf("  Decay model:   logistic(k=%.4f, t_mid=%.1f)\n",
           params->decay_rate_k, params->decay_midpoint);
    printf("  δ₀ (decay prob):       %.6f\n", decay_prob);
    printf("  δ₁ (wrong-dir prob):   %.6f\n", params->wrong_dir_prob);
    printf("  Ground state blocks:   %d bits\n", params->ground_state_block_size);
    printf("  Expected known frac:   ~%.4f\n",
           0.5 * (1.0 - decay_prob) + 0.5);
    printf("  Expected error rate:   ~%.4f\n",
           0.5 * decay_prob + 0.5 * params->wrong_dir_prob);
    printf("========================\n");
}

// on décide de l'état ground à donner à un bit
static int ground_state_for_bit(int bit_index, int block_size, int first_block_polarity) {
    int block_num = bit_index / block_size;
    return (block_num + first_block_polarity) % 2;
}


/*
 * decay_prob      delta0
 * wrong_dir_prob  delta1
 * bit_offset      offset pour le calul du ground
 */

// la dégradation pour un seul nombre gmp se fait ici, decay_prob correspond à la proba
// qu'un 1 devienne 0, et wrong dir à ce qu'un bit 0 devienne 1

static void degrade_component(
    const mpz_t original, mpz_t degraded, int *known_bits, int num_bits,
    double decay_prob, double wrong_dir_prob,
    int block_size, int first_polarity, int bit_offset,
    int *stats_flipped, int *stats_known)
{
    mpz_set(degraded, original);
    *stats_flipped = 0;
    *stats_known = 0;

    for (int i = 0; i < num_bits; i++) {
        int orig_bit = get_bit(original, i);
        int gs = ground_state_for_bit(i + bit_offset, block_size, first_polarity);
        int degraded_bit = orig_bit;

        if (orig_bit != gs) {
            // si le bit != gnd, on le flippe avec une proba delta0 (qu'un 1 devienne 0)
            if (rng_bernoulli(decay_prob)) {
                degraded_bit = gs;
                (*stats_flipped)++;
            }
        } else {
            // si le bit = gnd, on le flippe avec une proba delta1 (qu'un 0 devienne 1)
            if (rng_bernoulli(wrong_dir_prob)) {
                degraded_bit = 1 - gs;
                (*stats_flipped)++;
            }
        }

        // on applique le bit qu'on a (probablement) flippé
        if (degraded_bit)
            mpz_setbit(degraded, i);
        else
            mpz_clrbit(degraded, i);

      // l'attaquant connait le ground state, il ne peur pas se fier aux bits qui sont égaux au gnd
        if (degraded_bit != gs) {
            // diff de ground -> forcément correct
            known_bits[i] = degraded_bit;
            (*stats_known)++;
        } else {
            // égal au ground, on n'est pas trop sûr
            known_bits[i] = BIT_UNKNOWN;
        }
    }
}

void apply_decay(const rsa_key_t *original, degraded_key_t *dk,
                 const decay_params_t *params)
{
    double decay_prob = compute_decay_probability(params);
    int block_size = params->ground_state_block_size;

    rng_seed(params->seed);

    /* Random first block polarity */
    int first_polarity = params->random_ground_start ? (rng_bernoulli(0.5) ? 1 : 0) : 0;

    printf("[DECAY] Applying decay: T=%.1f°C, t=%.1fs, δ₀=%.6f, δ₁=%.6f\n",
           params->temperature_celsius, params->time_seconds,
           decay_prob, params->wrong_dir_prob);

    dk->key.bits = original->bits;
    dk->half_bits = original->bits / 2;
    dk->full_bits = original->bits;
    mpz_set(dk->key.N, original->N);
    mpz_set(dk->key.e, original->e);

    int total_flipped = 0, total_known = 0, total_bits = 0;
    int flipped, known;

    // on applique différents offsets pour que les états ground soient indépendants et qu'on parvienne à simuler
    // des zones mémoires non contigues
    int offsets[5];
    for (int i = 0; i < 5; i++) {
        offsets[i] = (int)(rng_uniform() * 10000);
    }

    // degradation de p
    degrade_component(original->p, dk->key.p, dk->known_p, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[0], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  p:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    // pour q
    degrade_component(original->q, dk->key.q, dk->known_q, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[1], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  q:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    // pour d
    degrade_component(original->d, dk->key.d, dk->known_d, dk->full_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[2], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->full_bits;
    printf("  d:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->full_bits, 100.0*known/dk->full_bits, flipped);

    // pour dp
    degrade_component(original->dp, dk->key.dp, dk->known_dp, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[3], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  dp: %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    // pour dq
    degrade_component(original->dq, dk->key.dq, dk->known_dq, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[4], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  dq: %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    dk->fraction_known = (double)total_known / total_bits;
    dk->error_rate     = (double)total_flipped / total_bits;

    printf("[DECAY] Overall: %.1f%% known (δ=%.4f), %.2f%% error rate\n",
           100.0 * dk->fraction_known, dk->fraction_known,
           100.0 * dk->error_rate);
}

// fonction pour dégradation simple
static void simple_degrade_component(const mpz_t orig, mpz_t deg, int *known, int nbits,
                                     double known_fraction, double error_rate) {
    mpz_set(deg, orig);
    for (int i = 0; i < nbits; i++) {
        int orig_bit = get_bit(orig, i);

        if (rng_bernoulli(known_fraction)) {
            if (rng_bernoulli(error_rate)) {
                int wrong_bit = 1 - orig_bit;
                if (wrong_bit) mpz_setbit(deg, i);
                else           mpz_clrbit(deg, i);
                known[i] = wrong_bit;
            } else {
                known[i] = orig_bit;
            }
        } else {
            known[i] = BIT_UNKNOWN;
        }
    }
}

void apply_simple_decay(const rsa_key_t *original, degraded_key_t *dk,
                        double known_fraction, double error_rate,
                        unsigned long seed)
{
    rng_seed(seed);

    dk->key.bits = original->bits;
    dk->half_bits = original->bits / 2;
    dk->full_bits = original->bits;
    mpz_set(dk->key.N, original->N);
    mpz_set(dk->key.e, original->e);

    printf("[DECAY-SIMPLE] known_fraction=%.4f, error_rate=%.4f\n",
           known_fraction, error_rate);

    simple_degrade_component(original->p,  dk->key.p,  dk->known_p,  dk->half_bits,  known_fraction, error_rate);
    simple_degrade_component(original->q,  dk->key.q,  dk->known_q,  dk->half_bits,  known_fraction, error_rate);
    simple_degrade_component(original->d,  dk->key.d,  dk->known_d,  dk->full_bits,  known_fraction, error_rate);
    simple_degrade_component(original->dp, dk->key.dp, dk->known_dp, dk->half_bits,  known_fraction, error_rate);
    simple_degrade_component(original->dq, dk->key.dq, dk->known_dq, dk->half_bits,  known_fraction, error_rate);

    dk->fraction_known = known_fraction;
    dk->error_rate = error_rate;
}
