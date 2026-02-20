#include "keygen.h"

#include <openssl/evp.h>
#include <openssl/core_names.h>
#include <openssl/param_build.h>
#include <openssl/bn.h>
#include <openssl/err.h>

/* ========================================================================
 *  keygen.c — RSA-CRT Key Generation via OpenSSL 3.0
 *
 *  Uses EVP_PKEY API (modern OpenSSL 3.x), extracts all CRT components,
 *  and converts them to GMP mpz_t via hex string intermediary.
 * ======================================================================== */

/* Convert an OpenSSL BIGNUM to a GMP mpz_t */
static int bn_to_mpz(const BIGNUM *bn, mpz_t out) {
    if (!bn) return -1;

    char *hex = BN_bn2hex(bn);
    if (!hex) return -1;

    /* GMP can parse hex with 0x prefix or directly */
    if (mpz_set_str(out, hex, 16) != 0) {
        OPENSSL_free(hex);
        return -1;
    }

    OPENSSL_free(hex);
    return 0;
}

/* Extract a BIGNUM parameter from an EVP_PKEY and store in mpz_t */
static int extract_param(EVP_PKEY *pkey, const char *param_name, mpz_t out) {
    BIGNUM *bn = NULL;

    if (EVP_PKEY_get_bn_param(pkey, param_name, &bn) != 1) {
        fprintf(stderr, "[KEYGEN] Failed to extract parameter '%s'\n", param_name);
        ERR_print_errors_fp(stderr);
        return -1;
    }

    int ret = bn_to_mpz(bn, out);
    BN_free(bn);
    return ret;
}

int rsa_keygen(rsa_key_t *key, int bits, unsigned long pub_exp) {
    EVP_PKEY_CTX *ctx = NULL;
    EVP_PKEY *pkey = NULL;
    BIGNUM *bn_e = NULL;
    int ret = -1;

    /* --- 1. Create key generation context --- */
    ctx = EVP_PKEY_CTX_new_from_name(NULL, "RSA", NULL);
    if (!ctx) {
        fprintf(stderr, "[KEYGEN] Failed to create EVP_PKEY_CTX\n");
        goto cleanup;
    }

    if (EVP_PKEY_keygen_init(ctx) <= 0) {
        fprintf(stderr, "[KEYGEN] Failed keygen init\n");
        goto cleanup;
    }

    /* --- 2. Set key parameters --- */
    if (EVP_PKEY_CTX_set_rsa_keygen_bits(ctx, bits) <= 0) {
        fprintf(stderr, "[KEYGEN] Failed to set key size\n");
        goto cleanup;
    }

    bn_e = BN_new();
    BN_set_word(bn_e, pub_exp);
    if (EVP_PKEY_CTX_set1_rsa_keygen_pubexp(ctx, bn_e) <= 0) {
        fprintf(stderr, "[KEYGEN] Failed to set public exponent\n");
        goto cleanup;
    }

    /* --- 3. Generate the key --- */
    if (EVP_PKEY_keygen(ctx, &pkey) <= 0) {
        fprintf(stderr, "[KEYGEN] Key generation failed\n");
        ERR_print_errors_fp(stderr);
        goto cleanup;
    }

    printf("[KEYGEN] Generated %d-bit RSA key with e = %lu\n", bits, pub_exp);

    /* --- 4. Extract all CRT components into GMP --- */
    key->bits = bits;

    /* OpenSSL 3.0 RSA parameter names */
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_N,           key->N)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_E,           key->e)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_D,           key->d)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_FACTOR1,     key->p)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_FACTOR2,     key->q)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_EXPONENT1,   key->dp)   != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_EXPONENT2,   key->dq)   != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_COEFFICIENT1, key->qinv) != 0) goto cleanup;

    /* Ensure p > q (convention: some implementations swap them) */
    if (mpz_cmp(key->p, key->q) < 0) {
        mpz_swap(key->p, key->q);
        /* Recompute dp, dq, qinv for swapped primes */
        mpz_t p1, q1;
        mpz_inits(p1, q1, NULL);
        mpz_sub_ui(p1, key->p, 1);
        mpz_sub_ui(q1, key->q, 1);
        mpz_mod(key->dp, key->d, p1);
        mpz_mod(key->dq, key->d, q1);
        mpz_invert(key->qinv, key->q, key->p);
        mpz_clears(p1, q1, NULL);
    }

    ret = 0;

cleanup:
    if (bn_e) BN_free(bn_e);
    if (pkey) EVP_PKEY_free(pkey);
    if (ctx)  EVP_PKEY_CTX_free(ctx);
    return ret;
}

int rsa_keygen_and_export(int bits, unsigned long pub_exp, const char *outdir) {
    rsa_key_t key;
    rsa_key_init(&key);

    if (rsa_keygen(&key, bits, pub_exp) != 0) {
        rsa_key_clear(&key);
        return -1;
    }

    /* Verify key before saving */
    if (!rsa_key_verify(&key)) {
        fprintf(stderr, "[KEYGEN] Generated key failed verification!\n");
        rsa_key_clear(&key);
        return -1;
    }
    printf("[KEYGEN] Key verification: OK\n");

    /* Save to file */
    char path[512];
    snprintf(path, sizeof(path), "%s/original_key.txt", outdir);
    int ret = rsa_key_export(&key, path);
    if (ret == 0) {
        printf("[KEYGEN] Key exported to %s\n", path);
    }

    rsa_key_print_summary(&key, "Generated Key");
    rsa_key_clear(&key);
    return ret;
}
