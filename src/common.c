#include "common.h"

// fichier avec les fonctions utilitaires




//init d'objets mpz a partir d'une struct de clé rsa générée
void rsa_key_init(rsa_key_t *key) {
    mpz_inits(key->N, key->e, key->d, key->p, key->q,
              key->dp, key->dq, key->qinv, NULL);
    key->bits = 0;
}

// liberation de la mémoire
void rsa_key_clear(rsa_key_t *key) {
    mpz_clears(key->N, key->e, key->d, key->p, key->q,
               key->dp, key->dq, key->qinv, NULL);
}

// init d'une clé dégradée
void degraded_key_init(degraded_key_t *dk, int key_bits) {
    rsa_key_init(&dk->key);
    dk->key.bits = key_bits;
    dk->half_bits = key_bits / 2;
    dk->full_bits = key_bits;

    /* Allocate known-bits arrays, initialize all to UNKNOWN */
    dk->known_p  = (int *)malloc(dk->half_bits * sizeof(int));
    dk->known_q  = (int *)malloc(dk->half_bits * sizeof(int));
    dk->known_d  = (int *)malloc(dk->full_bits * sizeof(int));
    dk->known_dp = (int *)malloc(dk->half_bits * sizeof(int));
    dk->known_dq = (int *)malloc(dk->half_bits * sizeof(int));

    for (int i = 0; i < dk->half_bits; i++) {
        dk->known_p[i]  = BIT_UNKNOWN;
        dk->known_q[i]  = BIT_UNKNOWN;
        dk->known_dp[i] = BIT_UNKNOWN;
        dk->known_dq[i] = BIT_UNKNOWN;
    }
    for (int i = 0; i < dk->full_bits; i++) {
        dk->known_d[i] = BIT_UNKNOWN;
    }

    dk->fraction_known = 0.0;
    dk->error_rate = 0.0;
}

void degraded_key_clear(degraded_key_t *dk) {
    rsa_key_clear(&dk->key);
    free(dk->known_p);
    free(dk->known_q);
    free(dk->known_d);
    free(dk->known_dp);
    free(dk->known_dq);
}

// recupere la valuation 2-adique (cb de 0 trainent à la fin pour savoir quelle puissance de 10 divise ce nombre)
int tau(const mpz_t x) {
    if (mpz_cmp_ui(x, 0) == 0) return 0;
    return (int)mpz_scan1(x, 0);   /* GMP: index of lowest set bit */
}


//wrapper the test bit pour avoir un int
int get_bit(const mpz_t num, int i) {
    return (int)mpz_tstbit(num, (mp_bitcnt_t)i);
}

/* --- RSA Key Verification --- */

bool rsa_key_verify(const rsa_key_t *key) {
    bool valid = true;
    mpz_t tmp, phi, p1, q1;
    mpz_inits(tmp, phi, p1, q1, NULL);

    /* Check N = p * q */
    mpz_mul(tmp, key->p, key->q);
    if (mpz_cmp(tmp, key->N) != 0) {
        fprintf(stderr, "[VERIFY] FAIL: N != p * q\n");
        valid = false;
    }

    // On verifie que e * d est congru à 1 mod phi n
    mpz_sub_ui(p1, key->p, 1);
    mpz_sub_ui(q1, key->q, 1);
    mpz_mul(phi, p1, q1);
    mpz_mul(tmp, key->e, key->d);
    mpz_mod(tmp, tmp, phi);
    if (mpz_cmp_ui(tmp, 1) != 0) {
        fprintf(stderr, "[VERIFY] FAIL: e*d != 1 mod phi(N)\n");
        valid = false;
    }


    // on verif que dp congru à d mod (p-1)
    mpz_mod(tmp, key->d, p1);
    if (mpz_cmp(tmp, key->dp) != 0) {
        fprintf(stderr, "[VERIFY] FAIL: dp != d mod (p-1)\n");
        valid = false;
    }

    // dq  congru à d mod q-1 ?
    mpz_mod(tmp, key->d, q1);
    if (mpz_cmp(tmp, key->dq) != 0) {
        fprintf(stderr, "[VERIFY] FAIL: dq != d mod (q-1)\n");
        valid = false;
    }

    // q inv * q congru à 1 mod p?
    mpz_mul(tmp, key->qinv, key->q);
    mpz_mod(tmp, tmp, key->p);
    if (mpz_cmp_ui(tmp, 1) != 0) {
        fprintf(stderr, "[VERIFY] FAIL: qinv*q != 1 mod p\n");
        valid = false;
    }

    mpz_clears(tmp, phi, p1, q1, NULL);
    return valid;
}

/* --- I/O Functions --- */

int rsa_key_export(const rsa_key_t *key, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("rsa_key_export: fopen");
        return -1;
    }

    fprintf(f, "bits=%d\n", key->bits);
    gmp_fprintf(f, "N=%Zd\n", key->N);
    gmp_fprintf(f, "e=%Zd\n", key->e);
    gmp_fprintf(f, "d=%Zd\n", key->d);
    gmp_fprintf(f, "p=%Zd\n", key->p);
    gmp_fprintf(f, "q=%Zd\n", key->q);
    gmp_fprintf(f, "dp=%Zd\n", key->dp);
    gmp_fprintf(f, "dq=%Zd\n", key->dq);
    gmp_fprintf(f, "qinv=%Zd\n", key->qinv);

    fclose(f);
    return 0;
}


//lire une clé
int rsa_key_import(rsa_key_t *key, const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("rsa_key_import: fopen");
        return -1;
    }

    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        /* Remove trailing newline */
        line[strcspn(line, "\n")] = '\0';

        char *eq = strchr(line, '=');
        if (!eq) continue;

        char *name = line;
        *eq = '\0';
        char *value = eq + 1;

        if (strcmp(name, "bits") == 0) {
            key->bits = atoi(value);
        } else if (strcmp(name, "N") == 0) {
            mpz_set_str(key->N, value, 10);
        } else if (strcmp(name, "e") == 0) {
            mpz_set_str(key->e, value, 10);
        } else if (strcmp(name, "d") == 0) {
            mpz_set_str(key->d, value, 10);
        } else if (strcmp(name, "p") == 0) {
            mpz_set_str(key->p, value, 10);
        } else if (strcmp(name, "q") == 0) {
            mpz_set_str(key->q, value, 10);
        } else if (strcmp(name, "dp") == 0) {
            mpz_set_str(key->dp, value, 10);
        } else if (strcmp(name, "dq") == 0) {
            mpz_set_str(key->dq, value, 10);
        } else if (strcmp(name, "qinv") == 0) {
            mpz_set_str(key->qinv, value, 10);
        }
    }

    fclose(f);
    return 0;
}

// exporter les bits qu'on connait à un fichier donné
static int export_known_bits(const int *bits, int count, const char *filepath) {
    FILE *f = fopen(filepath, "w");
    if (!f) { perror("export_known_bits"); return -1; }

    for (int i = 0; i < count; i++) {
        if (bits[i] != BIT_UNKNOWN) {
            fprintf(f, "%d %d\n", i, bits[i]);
        }
    }
    fclose(f);
    return 0;
}


//utilise la fct au dessus pour exporter toutes les sous infos de la clé
int degraded_key_export(const degraded_key_t *dk, const char *directory) {
    char path[512];

    snprintf(path, sizeof(path), "%s/degraded_key.txt", directory);
    rsa_key_export(&dk->key, path);


    snprintf(path, sizeof(path), "%s/known_bits_p.txt", directory);
    export_known_bits(dk->known_p, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_q.txt", directory);
    export_known_bits(dk->known_q, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_d.txt", directory);
    export_known_bits(dk->known_d, dk->full_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_dp.txt", directory);
    export_known_bits(dk->known_dp, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_dq.txt", directory);
    export_known_bits(dk->known_dq, dk->half_bits, path);

    // on exporte les statistiques
    snprintf(path, sizeof(path), "%s/decay_stats.txt", directory);
    FILE *f = fopen(path, "w");
    if (f) {
        fprintf(f, "fraction_known=%.6f\n", dk->fraction_known);
        fprintf(f, "error_rate=%.6f\n", dk->error_rate);
        fclose(f);
    }

    return 0;
}

// recup des bits depuis un fichier
static void import_known_bits_from_file(int *bits, int count, const char *filepath) {
    for (int i = 0; i < count; i++) bits[i] = BIT_UNKNOWN;
    FILE *f = fopen(filepath, "r");
    if (!f) return;
    int idx, val;
    while (fscanf(f, "%d %d", &idx, &val) == 2) {
        if (idx >= 0 && idx < count)
            bits[idx] = val;
    }
    fclose(f);
}

int degraded_key_import(degraded_key_t *dk, const char *directory) {
    char path[512];

    /* Import degraded key values */
    snprintf(path, sizeof(path), "%s/degraded_key.txt", directory);
    if (rsa_key_import(&dk->key, path) != 0) return -1;

    /* Re-init arrays based on imported key size */
    dk->half_bits = dk->key.bits / 2;
    dk->full_bits = dk->key.bits;

    snprintf(path, sizeof(path), "%s/known_bits_p.txt", directory);
    import_known_bits_from_file(dk->known_p, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_q.txt", directory);
    import_known_bits_from_file(dk->known_q, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_d.txt", directory);
    import_known_bits_from_file(dk->known_d, dk->full_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_dp.txt", directory);
    import_known_bits_from_file(dk->known_dp, dk->half_bits, path);

    snprintf(path, sizeof(path), "%s/known_bits_dq.txt", directory);
    import_known_bits_from_file(dk->known_dq, dk->half_bits, path);

    return 0;
}

//affichage

void rsa_key_print_summary(const rsa_key_t *key, const char *label) {
    printf("=== %s (%d-bit RSA-CRT Key) ===\n", label, key->bits);
    gmp_printf("  N    = %Zd\n", key->N);
    gmp_printf("  e    = %Zd\n", key->e);
    gmp_printf("  d    = %Zd\n", key->d);
    gmp_printf("  p    = %Zd\n", key->p);
    gmp_printf("  q    = %Zd\n", key->q);
    gmp_printf("  dp   = %Zd\n", key->dp);
    gmp_printf("  dq   = %Zd\n", key->dq);
    gmp_printf("  qinv = %Zd\n", key->qinv);
    printf("=================================\n");
}

void degraded_key_print_stats(const degraded_key_t *dk) {
    // on compte le nb de bits connus par composant
    int known_p = 0, known_q = 0, known_d = 0, known_dp = 0, known_dq = 0;
    int total_known = 0, total_bits = 0;

    for (int i = 0; i < dk->half_bits; i++) {
        if (dk->known_p[i] != BIT_UNKNOWN)  known_p++;
        if (dk->known_q[i] != BIT_UNKNOWN)  known_q++;
        if (dk->known_dp[i] != BIT_UNKNOWN) known_dp++;
        if (dk->known_dq[i] != BIT_UNKNOWN) known_dq++;
    }
    for (int i = 0; i < dk->full_bits; i++) {
        if (dk->known_d[i] != BIT_UNKNOWN) known_d++;
    }

    total_known = known_p + known_q + known_d + known_dp + known_dq;
    total_bits  = 4 * dk->half_bits + dk->full_bits;

    printf("=== Degraded Key Statistics ===\n");
    printf("  Key size:     %d bits\n", dk->key.bits);
    printf("  Known bits:\n");
    printf("    p:   %d / %d  (%.1f%%)\n", known_p,  dk->half_bits, 100.0*known_p/dk->half_bits);
    printf("    q:   %d / %d  (%.1f%%)\n", known_q,  dk->half_bits, 100.0*known_q/dk->half_bits);
    printf("    d:   %d / %d  (%.1f%%)\n", known_d,  dk->full_bits, 100.0*known_d/dk->full_bits);
    printf("    dp:  %d / %d  (%.1f%%)\n", known_dp, dk->half_bits, 100.0*known_dp/dk->half_bits);
    printf("    dq:  %d / %d  (%.1f%%)\n", known_dq, dk->half_bits, 100.0*known_dq/dk->half_bits);
    printf("  Total: %d / %d  (%.1f%%)  [delta = %.4f]\n",
           total_known, total_bits, 100.0*total_known/total_bits,
           (double)total_known / total_bits);
    printf("  Error rate:   %.4f\n", dk->error_rate);
    printf("===============================\n");
}
