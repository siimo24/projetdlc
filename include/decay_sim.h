#ifndef DECAY_SIM_H
#define DECAY_SIM_H

#include "common.h"

/* ========================================================================
 *  decay_sim.h — DRAM Remanence Decay Simulation
 *
 *  Models cold boot attack scenarios based on Halderman et al. (HSH08):
 *   - Binary asymmetric channel (ground state bias)
 *   - Temperature-dependent logistic decay curve
 *   - Alternating ground state block regions
 *   - Both erasure view (Heninger-Shacham) and error view (Henecka-May-Meurer)
 * ======================================================================== */

/*
 * Initialize default decay parameters for a given temperature scenario.
 *
 * Presets:
 *   "room"    — Room temperature (~25°C), fast decay
 *   "cold"    — Cooled with compressed air (~-50°C)
 *   "frozen"  — Liquid nitrogen (~-196°C), minimal decay
 *   "custom"  — Returns zeroed params for manual configuration
 */
void decay_params_preset(decay_params_t *params, const char *preset);

/*
 * Apply DRAM decay to an RSA key, producing a degraded key.
 *
 * The simulation:
 *   1. Assigns ground state regions (alternating blocks of 0s and 1s)
 *   2. For each bit in each component (p, q, d, dp, dq):
 *      - If bit disagrees with ground state: flip with probability decay_prob(t, T)
 *      - If bit agrees with ground state: flip with tiny probability (wrong_dir_prob)
 *   3. For erasure view: bits that agree with ground state are marked as "known"
 *      (since in a real attack, you can determine the decay direction)
 *   4. Tracks error rate (flipped bits) and known fraction statistics
 *
 * @param original   The original correct RSA key
 * @param dk         Output: degraded key with known-bits arrays
 * @param params     Decay parameters (temperature, time, etc.)
 */
void apply_decay(const rsa_key_t *original, degraded_key_t *dk,
                 const decay_params_t *params);

/*
 * Apply simple uniform decay (flat probability, no ground state model).
 * Useful for quick testing or comparison with the realistic model.
 *
 * @param original       The original correct RSA key
 * @param dk             Output: degraded key
 * @param known_fraction Fraction of bits to mark as known (for erasure view)
 * @param error_rate     Fraction of known bits that are wrong (for error view)
 * @param seed           RNG seed
 */
void apply_simple_decay(const rsa_key_t *original, degraded_key_t *dk,
                        double known_fraction, double error_rate,
                        unsigned long seed);

/*
 * Compute the decay probability for a given time and temperature.
 * Uses logistic model: P(decay) = 1 / (1 + exp(-k * (t - t_mid)))
 */
double compute_decay_probability(const decay_params_t *params);

/*
 * Print decay parameters summary.
 */
void decay_params_print(const decay_params_t *params);

#endif /* DECAY_SIM_H */
