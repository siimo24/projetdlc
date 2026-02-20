#ifndef COMMON_H
#define COMMON_H

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

/* ========================================================================
 *  RSA Key Reconstruction from Degraded Memory
 *  Based on: Heninger-Shacham (HS09), Henecka-May-Meurer (HMM10),
 *            Halderman et al. (HSH08)
 * ======================================================================== */

/* --- RSA-CRT Key Structure --- */
typedef struct {
    mpz_t N;        /* Modulus N = p * q */
    mpz_t e;        /* Public exponent */
    mpz_t d;        /* Private exponent */
    mpz_t p;        /* First prime factor */
    mpz_t q;        /* Second prime factor */
    mpz_t dp;       /* d mod (p-1)  (CRT exponent) */
    mpz_t dq;       /* d mod (q-1)  (CRT exponent) */
    mpz_t qinv;     /* q^{-1} mod p (CRT coefficient) */
    int bits;       /* Key size in bits (e.g., 512, 1024, 2048) */
} rsa_key_t;

/* --- Degraded Key with Bit Knowledge --- */
typedef enum {
    BIT_UNKNOWN = -1,   /* Erasure: bit value unknown */
    BIT_ZERO    =  0,   /* Known to be 0 */
    BIT_ONE     =  1    /* Known to be 1 */
} bit_status_t;

typedef struct {
    rsa_key_t key;          /* The degraded key values (may contain errors) */

    /* Erasure view (for Heninger-Shacham): -1 = unknown, 0/1 = known */
    int *known_p;
    int *known_q;
    int *known_d;
    int *known_dp;
    int *known_dq;

    /* Component sizes in bits */
    int half_bits;          /* n/2 bits for p, q, dp, dq */
    int full_bits;          /* n bits for d */

    /* Decay statistics */
    double fraction_known;  /* Fraction of bits with known values */
    double error_rate;      /* Fraction of bits with wrong values */
} degraded_key_t;

/* --- Decay Simulation Parameters --- */
typedef struct {
    double temperature_celsius;     /* Temperature after power loss */
    double time_seconds;            /* Time since power loss */
    double wrong_dir_prob;          /* Probability of wrong-direction flip (~0.001) */

    /* Logistic model parameters (fitted per temperature) */
    double decay_midpoint;          /* Time at which decay = 50% */
    double decay_rate_k;            /* Steepness of logistic curve */

    /* Ground state configuration */
    int ground_state_block_size;    /* Size of alternating ground state regions */
    bool random_ground_start;       /* Random first block polarity */

    /* RNG seed */
    unsigned long seed;
} decay_params_t;

/* --- common.c: Utility Functions --- */

/* Initialize an RSA key structure */
void rsa_key_init(rsa_key_t *key);

/* Free an RSA key structure */
void rsa_key_clear(rsa_key_t *key);

/* Initialize a degraded key structure */
void degraded_key_init(degraded_key_t *dk, int key_bits);

/* Free a degraded key structure */
void degraded_key_clear(degraded_key_t *dk);

/* Compute tau(x) = largest power of 2 dividing x */
int tau(const mpz_t x);

/* Get bit i of a GMP number (returns 0 or 1) */
int get_bit(const mpz_t num, int i);

/* Export RSA key to a text file */
int rsa_key_export(const rsa_key_t *key, const char *filename);

/* Import RSA key from a text file */
int rsa_key_import(rsa_key_t *key, const char *filename);

/* Export degraded key (known bits arrays) to files */
int degraded_key_export(const degraded_key_t *dk, const char *directory);

/* Import degraded key from files */
int degraded_key_import(degraded_key_t *dk, const char *directory);

/* Print key summary */
void rsa_key_print_summary(const rsa_key_t *key, const char *label);

/* Print degraded key statistics */
void degraded_key_print_stats(const degraded_key_t *dk);

/* Verify RSA key consistency (ed ≡ 1 mod φ(N), N = pq, etc.) */
bool rsa_key_verify(const rsa_key_t *key);

#endif /* COMMON_H */
