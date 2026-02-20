#include "decay_sim.h"

/* ========================================================================
 *  decay_sim.c — DRAM Remanence Decay Simulation
 *
 *  Realistic model based on Halderman et al. (HSH08) experimental data:
 *
 *  Key observations from the paper:
 *   - Memory cells have a "ground state" (preferred value after power loss)
 *   - Ground states alternate in blocks across memory
 *   - Bits opposing ground state flip to ground state with probability δ₀
 *   - Bits matching ground state flip away with tiny probability δ₁ ≈ 0.001
 *   - δ₀ follows a logistic curve as a function of time
 *   - Temperature dramatically affects decay rate:
 *       Room temp (~25°C):   full decay in ~2.5–35 seconds
 *       Cooled (~-50°C):     < 0.001% error after 60 seconds
 *       LN₂ (~-196°C):      0.17% error after 60 minutes
 *
 *  Erasure view (for HS09):
 *   - Attacker detects ground state direction per region
 *   - Bits opposing ground state in degraded image → KNOWN (almost certainly original)
 *   - Bits matching ground state → UNKNOWN (could be original or decayed)
 *
 *  Error view (for HMM10):
 *   - All bits present, but some have wrong values
 *   - Error rate ≈ 0.5 * δ₀  (half the bits oppose ground state, fraction δ₀ flipped)
 * ======================================================================== */

/* ---- Random number generation (simple LCG for portability) ---- */

static unsigned long rng_state;

static void rng_seed(unsigned long seed) {
    rng_state = seed ? seed : (unsigned long)time(NULL);
}

/* Returns uniform random in [0, 1) */
static double rng_uniform(void) {
    /* xorshift64 — fast and reasonable quality */
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state & 0x7FFFFFFFFFFFFFFF) / (double)0x7FFFFFFFFFFFFFFF;
}

/* Returns true with given probability */
static bool rng_bernoulli(double p) {
    return rng_uniform() < p;
}

/* ---- Logistic decay model ---- */

/*
 * Logistic model parameters fitted to HSH08 experimental data.
 *
 * The decay probability (fraction of anti-ground-state bits that flip)
 * follows: δ₀(t) = 1 / (1 + exp(-k * (t - t_mid)))
 *
 * Parameters per temperature regime:
 *                    k (rate)     t_mid (seconds)
 *   Room (25°C):     0.5          10
 *   Cold (-50°C):    0.01         600
 *   LN₂ (-196°C):   0.0003       36000
 *
 * These are approximate fits to the paper's experimental curves.
 */

void decay_params_preset(decay_params_t *params, const char *preset) {
    memset(params, 0, sizeof(*params));

    /* Defaults for all presets */
    params->wrong_dir_prob = 0.001;         /* δ₁: wrong-direction flip rate */
    params->ground_state_block_size = 256;  /* Alternating block size (bits) */
    params->random_ground_start = true;
    params->seed = 0;                       /* 0 = use current time */

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
        params->time_seconds = 3600.0;  /* 60 minutes */
        params->decay_rate_k = 0.0003;
        params->decay_midpoint = 36000.0;
    }
    else {
        /* "custom" — zeroed for manual configuration */
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

/* ---- Ground state assignment ---- */

/*
 * Determine ground state for a given bit position.
 * Simulates alternating blocks of ground_state = 0 and ground_state = 1.
 * The first block's polarity is determined randomly if random_ground_start is set.
 */
static int ground_state_for_bit(int bit_index, int block_size, int first_block_polarity) {
    int block_num = bit_index / block_size;
    return (block_num + first_block_polarity) % 2;
}

/* ---- Core decay function for a single component ---- */

/*
 * Degrade a single GMP number and produce known-bits array.
 *
 * @param original        Original correct value
 * @param degraded        Output: degraded value (with flipped bits)
 * @param known_bits      Output: known-bits array (-1=unknown, 0/1=known)
 * @param num_bits        Number of bits in this component
 * @param decay_prob      Probability of anti-ground-state bit flipping (δ₀)
 * @param wrong_dir_prob  Probability of ground-state bit flipping (δ₁)
 * @param block_size      Ground state alternating block size
 * @param first_polarity  Polarity of the first ground state block
 * @param bit_offset      Offset for ground state calculation (to vary per component)
 * @param stats_flipped   Output: count of bits that were flipped
 * @param stats_known     Output: count of bits marked as known
 */
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
            /* Bit opposes ground state → flip to ground state with probability δ₀ */
            if (rng_bernoulli(decay_prob)) {
                degraded_bit = gs;
                (*stats_flipped)++;
            }
        } else {
            /* Bit matches ground state → wrong-direction flip with probability δ₁ */
            if (rng_bernoulli(wrong_dir_prob)) {
                degraded_bit = 1 - gs;
                (*stats_flipped)++;
            }
        }

        /* Apply the (possibly flipped) bit */
        if (degraded_bit)
            mpz_setbit(degraded, i);
        else
            mpz_clrbit(degraded, i);

        /*
         * Erasure view (attacker perspective):
         * The attacker determines the ground state direction for each region.
         * - If degraded bit OPPOSES ground state → it's almost certainly the
         *   original value (wrong-direction flips are negligible).
         *   Mark as KNOWN.
         * - If degraded bit MATCHES ground state → could be original or decayed.
         *   Mark as UNKNOWN.
         */
        if (degraded_bit != gs) {
            /* Anti-ground-state → known correct */
            known_bits[i] = degraded_bit;
            (*stats_known)++;
        } else {
            /* Matches ground state → ambiguous */
            known_bits[i] = BIT_UNKNOWN;
        }
    }
}

/* ---- Public API ---- */

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

    /* Each component gets a different bit_offset so ground state patterns
     * are independent (simulating different memory locations) */
    int offsets[5];
    for (int i = 0; i < 5; i++) {
        offsets[i] = (int)(rng_uniform() * 10000);
    }

    /* Degrade p */
    degrade_component(original->p, dk->key.p, dk->known_p, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[0], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  p:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    /* Degrade q */
    degrade_component(original->q, dk->key.q, dk->known_q, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[1], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  q:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    /* Degrade d */
    degrade_component(original->d, dk->key.d, dk->known_d, dk->full_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[2], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->full_bits;
    printf("  d:  %d/%d known (%.1f%%), %d flipped\n",
           known, dk->full_bits, 100.0*known/dk->full_bits, flipped);

    /* Degrade dp */
    degrade_component(original->dp, dk->key.dp, dk->known_dp, dk->half_bits,
                      decay_prob, params->wrong_dir_prob, block_size, first_polarity,
                      offsets[3], &flipped, &known);
    total_flipped += flipped; total_known += known; total_bits += dk->half_bits;
    printf("  dp: %d/%d known (%.1f%%), %d flipped\n",
           known, dk->half_bits, 100.0*known/dk->half_bits, flipped);

    /* Degrade dq */
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

/* Helper for simple flat-probability decay */
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
