#include "decay_sim.h"

// rng local
static unsigned long rng_state;

static void rng_seed(unsigned long seed) {
    rng_state = seed ? seed : (unsigned long)time(NULL);
}

// sort un double entre 0 et 1 (algo xorshift64, basique mais rapide et fait le taff)
static double rng_uniform(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state & 0x7FFFFFFFFFFFFFFF) / (double)0x7FFFFFFFFFFFFFFF;
}

// un pile ou face truqué avec proba p (pour le % de decay)
static bool rng_bernoulli(double p) {
    return rng_uniform() < p;
}

// les presets tirés du papier HSH08 (température et tps)
void decay_params_preset(decay_params_t *params, const char *preset) {
    memset(params, 0, sizeof(*params));

    params->wrong_dir_prob = 0.001;         // le fameux delta1 (erreur thermique random)
    params->ground_state_block_size = 256;  // taille d'un bloc memoire avec le meme gnd
    params->random_ground_start = true;
    params->seed = 0;                       // 0 = on tape sur l'heure de l'OS

    if (strcmp(preset, "room") == 0) {
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
        params->time_seconds = 3600.0;  // 1h sous azote liquide
        params->decay_rate_k = 0.0003;
        params->decay_midpoint = 36000.0;
    }
    else {
        // fallback si le mec rentre nimp ou "custom"
        params->temperature_celsius = 25.0;
        params->time_seconds = 0.0;
        params->decay_rate_k = 0.5;
        params->decay_midpoint = 10.0;
    }
}




// calcul proba de decay basee sur la temperature et le tps (la fameuse courbe logistique)
double compute_decay_probability(const decay_params_t *params) {
    double exponent = -params->decay_rate_k * (params->time_seconds - params->decay_midpoint);

    // on clamp pcq les floats en C ca part vite en overflow si l'expo est enorme
    if (exponent > 500.0)  return 0.0;
    if (exponent < -500.0) return 1.0;

    return 1.0 / (1.0 + exp(exponent));
}

// ptit print des stats de simu
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

// check si la zone memoire a un etat de repos (gnd) à 0 ou à 1
static int ground_state_for_bit(int bit_index, int block_size, int first_block_polarity) {
    int block_num = bit_index / block_size;
    return (block_num + first_block_polarity) % 2;
}


// fct core de la simu: elle detruit un entier GMP bit par bit en fct de son gnd
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
            // bit opposé au ground -> il perd son jus avec proba delta0
            if (rng_bernoulli(decay_prob)) {
                degraded_bit = gs;
                (*stats_flipped)++;
            }
        } else {
            // bit dejà sur le ground -> rare bruit thermique (delta1)
            if (rng_bernoulli(wrong_dir_prob)) {
                degraded_bit = 1 - gs;
                (*stats_flipped)++;
            }
        }

        // on ecrit le bit (peut-etre) flippé
        if (degraded_bit)
            mpz_setbit(degraded, i);
        else
            mpz_clrbit(degraded, i);

        // POV de l'attaquant: s'il voit un bit opposé au gnd, il sait que c'est legit
        if (degraded_bit != gs) {
            known_bits[i] = degraded_bit; // bingo
            (*stats_known)++;
        } else {
            // bit confondu avec le gnd -> effacement
            known_bits[i] = BIT_UNKNOWN;
        }
    }
}

// applique la simu a toute la clef entiere
void apply_decay(const rsa_key_t *original, degraded_key_t *dk,
                 const decay_params_t *params)
{
    double decay_prob = compute_decay_probability(params);
    int block_size = params->ground_state_block_size;

    rng_seed(params->seed);

    // on tire au pif si la ram physique commence par un bloc de 0 ou de 1
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

    // on met des offsets random pcq les vars sont physiquement eparpillées ds la puce
    int offsets[5];
    for (int i = 0; i < 5; i++) {
        offsets[i] = (int)(rng_uniform() * 10000);
    }

    // on degrade les 5 morceaux un par un
    degrade_component(original->p, dk->key.p, dk->known_p, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[0], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  p:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    degrade_component(original->q, dk->key.q, dk->known_q, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[1], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  q:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    degrade_component(original->d, dk->key.d, dk->known_d, dk->full_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[2], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->full_bits;
    printf("  d:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->full_bits, 100.0*known/dk->full_bits, flipped);

    degrade_component(original->dp, dk->key.dp, dk->known_dp, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[3], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  dp: %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

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

// ==========================================================
// Mode simple : pour test d'algo brut (sans physique de ram)
// ==========================================================
static void simple_degrade_component(const mpz_t orig, mpz_t deg, int *known, int nbits,
                                     double known_fraction, double error_rate) {
    mpz_set(deg, orig);
    for (int i = 0; i < nbits; i++) {
        int orig_bit = get_bit(orig, i);

        // on le garde ?
        if (rng_bernoulli(known_fraction)) {
            // oui, mais est-ce qu'il est foiré par une erreur ?
            if (rng_bernoulli(error_rate)) {
                int wrong_bit = 1 - orig_bit;
                if (wrong_bit) mpz_setbit(deg, i);
                else           mpz_clrbit(deg, i);
                known[i] = wrong_bit;
            } else {
                known[i] = orig_bit;
            }
        } else {
            // effacé
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