#include "keygen.h"

#include <openssl/evp.h>
#include <openssl/core_names.h>
#include <openssl/param_build.h>
#include <openssl/bn.h>
#include <openssl/err.h>

// gen de clefs RSA-CRT via l'API moderne d'OpenSSL (v3)

// passe d'un BIGNUM openssl a un mpz_t de gmp en passant par une string hexa
static int bn_to_mpz(const BIGNUM *bn, mpz_t out) {
    if (!bn) return -1;

    char *hex = BN_bn2hex(bn);
    if (!hex) return -1;

    // gmp se debrouille tout seul avec l'hexa
    if (mpz_set_str(out, hex, 16) != 0) {
        OPENSSL_free(hex);
        return -1;
    }

    OPENSSL_free(hex);
    return 0;
}

// helper pour recup un param precis de la clef openssl
static int extract_param(EVP_PKEY *pkey, const char *param_name, mpz_t out) {
    BIGNUM *bn = NULL;

    if (EVP_PKEY_get_bn_param(pkey, param_name, &bn) != 1) {
        fprintf(stderr, "[KEYGEN] aie, impossible d'extraire '%s'\n", param_name);
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

    // 1. init du context openssl
    ctx = EVP_PKEY_CTX_new_from_name(NULL, "RSA", NULL);
    if (!ctx) {
        fprintf(stderr, "[KEYGEN] echec creation ctx EVP_PKEY\n");
        goto cleanup;
    }

    if (EVP_PKEY_keygen_init(ctx) <= 0) {
        fprintf(stderr, "[KEYGEN] pb init keygen\n");
        goto cleanup;
    }

    // 2. config de la taille et de l'exposant public e
    if (EVP_PKEY_CTX_set_rsa_keygen_bits(ctx, bits) <= 0) {
        fprintf(stderr, "[KEYGEN] echec setup taille clef\n");
        goto cleanup;
    }

    bn_e = BN_new();
    BN_set_word(bn_e, pub_exp);
    if (EVP_PKEY_CTX_set1_rsa_keygen_pubexp(ctx, bn_e) <= 0) {
        fprintf(stderr, "[KEYGEN] echec setup exposant e\n");
        goto cleanup;
    }

    // 3. magie: on genere la clef
    if (EVP_PKEY_keygen(ctx, &pkey) <= 0) {
        fprintf(stderr, "[KEYGEN] generation hs\n");
        ERR_print_errors_fp(stderr);
        goto cleanup;
    }

    printf("[KEYGEN] clef RSA %d-bit generee avec e = %lu\n", bits, pub_exp);

    // 4. on extrait tout le bordel CRT vers nos var GMP
    key->bits = bits;

    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_N,           key->N)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_E,           key->e)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_D,           key->d)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_FACTOR1,     key->p)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_FACTOR2,     key->q)    != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_EXPONENT1,   key->dp)   != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_EXPONENT2,   key->dq)   != 0) goto cleanup;
    if (extract_param(pkey, OSSL_PKEY_PARAM_RSA_COEFFICIENT1, key->qinv) != 0) goto cleanup;

    // FIX: on s'assure que p > q sinon ca casse l'algo de recup plus tard
    if (mpz_cmp(key->p, key->q) < 0) {
        mpz_swap(key->p, key->q);

        // du coup faut recalculer dp, dq et qinv a la mano
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
    // menage std
    if (bn_e) BN_free(bn_e);
    if (pkey) EVP_PKEY_free(pkey);
    if (ctx)  EVP_PKEY_CTX_free(ctx);
    return ret;
}

// wrapper pour tout gerer d'un coup (gen, verif, ecriture)
int rsa_keygen_and_export(int bits, unsigned long pub_exp, const char *outdir) {
    rsa_key_t key;
    rsa_key_init(&key);

    if (rsa_keygen(&key, bits, pub_exp) != 0) {
        rsa_key_clear(&key);
        return -1;
    }

    // check rapide que la clef marche avnt de sauv
    if (!rsa_key_verify(&key)) {
        fprintf(stderr, "[KEYGEN] aie, echec verif clef!\n");
        rsa_key_clear(&key);
        return -1;
    }
    printf("[KEYGEN] verif clef: OK\n");

    // sauvegarde sur le disque
    char path[512];
    snprintf(path, sizeof(path), "%s/original_key.txt", outdir);
    int ret = rsa_key_export(&key, path);
    if (ret == 0) {
        printf("[KEYGEN] clef exportee ds %s\n", path);
    }

    rsa_key_print_summary(&key, "Clef Generee");
    rsa_key_clear(&key);
    return ret;
}