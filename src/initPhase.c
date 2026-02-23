#include "init_phase.h"

// Phase d'init de l'attaque. On recup k, kp et kq avant de lancer l'arbre de recherche.

/* --- Fonctions d'init basiques --- */

void init_result_init(init_result_t *res) {
    // init des var GMP, classique
    mpz_inits(res->k, res->kp, res->kq, NULL);
    mpz_inits(res->d_msb_corrected, NULL);
    mpz_inits(res->p0, res->q0, res->d0, res->dp0, res->dq0, NULL);
    res->tau_k = res->tau_kp = res->tau_kq = 0;
    res->valid = false;
    res->kp_kq_swapped = false; // au cas ou on les inverse plus tard
}

void init_result_clear(init_result_t *res) {
    mpz_clears(res->k, res->kp, res->kq, NULL);
    mpz_clears(res->d_msb_corrected, NULL);
    mpz_clears(res->p0, res->q0, res->d0, res->dp0, res->dq0, NULL);
}

void init_result_print(const init_result_t *res) {
    // affichage de debug, tkt c'est juste du print
    printf("=== Init Phase Results ===\n");
    gmp_printf("  k    = %Zd\n", res->k);
    gmp_printf("  kp   = %Zd\n", res->kp);
    gmp_printf("  kq   = %Zd\n", res->kq);
    printf("  τ(k) = %d,  τ(kp) = %d,  τ(kq) = %d\n",
           res->tau_k, res->tau_kp, res->tau_kq);
    printf("  Valid: %s\n", res->valid ? "YES" : "NO");
    if (res->kp_kq_swapped)
        printf("  (kp/kq assignment swapped)\n");

    printf("\n  Initial partial values:\n");
    gmp_printf("    p0  = %Zd\n", res->p0);
    gmp_printf("    q0  = %Zd\n", res->q0);
    gmp_printf("    d0  = %Zd\n", res->d0);
    gmp_printf("    dp0 = %Zd\n", res->dp0);
    gmp_printf("    dq0 = %Zd\n", res->dq0);
    printf("  Slice(0): p[0]=%d, q[0]=%d, d[%d]=%d, dp[%d]=%d, dq[%d]=%d\n",
           get_bit(res->p0, 0), get_bit(res->q0, 0),
           res->tau_k, get_bit(res->d0, res->tau_k),
           res->tau_kp, get_bit(res->dp0, res->tau_kp),
           res->tau_kq, get_bit(res->dq0, res->tau_kq));
    printf("==========================\n");
}


// on bruteforce k. e est petit (genre 65537) donc ca va super vite.
// On compare les MSB (poids fort) calculés avec ceux du dump
int find_k(const mpz_t N, const mpz_t e, const mpz_t degraded_d,
           mpz_t out_k, mpz_t out_d_tilde, bool verbose)
{
    mpz_t k_candidate, d_tilde, temp, N_plus_1;
    mpz_inits(k_candidate, d_tilde, temp, N_plus_1, NULL);

    mpz_add_ui(N_plus_1, N, 1);

    size_t total_bits = mpz_sizeinbase(degraded_d, 2);
    size_t msb_count = (total_bits / 2) - 2; // on check que la moitié haute

    int best_match = 0;
    int best_match_k_count = 0;

    if (verbose) {
        printf("[INIT] Finding k: comparing top %zu bits of d̃(k') with degraded d\n",
               msb_count);
        printf("[INIT] Searching k' in [1, e-1] ...\n");
    }

    // boucle sur tous les k possibles
    for (mpz_set_ui(k_candidate, 1);
         mpz_cmp(k_candidate, e) < 0;
         mpz_add_ui(k_candidate, k_candidate, 1))
    {
        // calcul de l'approximation de d
        mpz_mul(temp, k_candidate, N_plus_1);
        mpz_add_ui(temp, temp, 1);
        mpz_fdiv_q(d_tilde, temp, e);

        // on compte cb de bits matchent
        int matching = 0;
        for (size_t i = 0; i < msb_count; i++) {
            size_t bit_idx = total_bits - 1 - i;
            if (mpz_tstbit(d_tilde, bit_idx) == mpz_tstbit(degraded_d, bit_idx))
                matching++;
        }

        // maj du meilleur candidat
        if (matching > best_match) {
            best_match = matching;
            mpz_set(out_k, k_candidate);
            mpz_set(out_d_tilde, d_tilde);
            best_match_k_count = 1;
        } else if (matching == best_match) {
            best_match_k_count++;
        }
    }

    if (verbose) {
        gmp_printf("[INIT] Best k = %Zd  (%d / %zu MSBs match, %.1f%%)\n",
                   out_k, best_match, msb_count,
                   100.0 * best_match / msb_count);
        if (best_match_k_count > 1) {
            printf("[INIT] WARNING: %d candidates tied for best match\n",
                   best_match_k_count);
        }
    }

    mpz_clears(k_candidate, d_tilde, temp, N_plus_1, NULL);
    return best_match;
}


//  fusionne les bons MSB qu'on vient de calculer avec les LSB tout pétés de la RAM
void correct_d_msb(const mpz_t d_tilde, const mpz_t degraded_d,
                   mpz_t out_d, int n_bits)
{
    mpz_t msb_mask, lsb_mask, msb_part, lsb_part;
    mpz_inits(msb_mask, lsb_mask, msb_part, lsb_part, NULL);

    size_t total_bits = (size_t)n_bits;
    if (total_bits == 0)
        total_bits = mpz_sizeinbase(d_tilde, 2);

    size_t msb_bits = (total_bits / 2) - 2;
    size_t lsb_bits = total_bits - msb_bits;

    // creation des masques pour decouper les bits
    mpz_set_ui(msb_mask, 1);
    mpz_mul_2exp(msb_mask, msb_mask, msb_bits);
    mpz_sub_ui(msb_mask, msb_mask, 1);
    mpz_mul_2exp(msb_mask, msb_mask, lsb_bits);

    mpz_set_ui(lsb_mask, 1);
    mpz_mul_2exp(lsb_mask, lsb_mask, lsb_bits);
    mpz_sub_ui(lsb_mask, lsb_mask, 1);

    // on assemble le frankenstein
    mpz_and(msb_part, d_tilde, msb_mask);
    mpz_and(lsb_part, degraded_d, lsb_mask);
    mpz_ior(out_d, msb_part, lsb_part);

    mpz_clears(msb_mask, lsb_mask, msb_part, lsb_part, NULL);
}


// Algo de Tonelli-Shanks. relou à coder mais obligatoire pour les racines carr mod p
int tonelli_shanks(const mpz_t n, const mpz_t p, mpz_t root1, mpz_t root2) {
    // check si n est un residu quadratique sinon c'est mort d'avance
    if (mpz_legendre(n, p) != 1) {
        return -1;
    }

    mpz_t Q, z, M, c, t, R, temp, b, exponent;
    mpz_inits(Q, z, M, c, t, R, temp, b, exponent, NULL);

    mpz_sub_ui(Q, p, 1);
    unsigned long S = 0;
    while (mpz_even_p(Q)) {
        mpz_fdiv_q_2exp(Q, Q, 1);
        S++;
    }

    // recherche du non-residu
    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) {
        mpz_add_ui(z, z, 1);
    }

    // setup
    mpz_powm(c, z, Q, p);
    mpz_powm(t, n, Q, p);
    mpz_add_ui(temp, Q, 1);
    mpz_fdiv_q_2exp(temp, temp, 1);
    mpz_powm(R, n, temp, p);
    mpz_set_ui(M, S);

     // boucle principale
    while (mpz_cmp_ui(t, 0) != 0 && mpz_cmp_ui(t, 1) != 0) {
        unsigned long i = 0;
        mpz_set(temp, t);
        while (mpz_cmp_ui(temp, 1) != 0) {
            mpz_powm_ui(temp, temp, 2, p);
            i++;
        }

        if (i == mpz_get_ui(M)) {
            mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
            return -1; // erreur chelou
        }

        unsigned long exp_val = mpz_get_ui(M) - i - 1;
        mpz_ui_pow_ui(exponent, 2, exp_val);
        mpz_powm(b, c, exponent, p);

        mpz_powm_ui(c, b, 2, p);
        mpz_mul(t, t, c);
        mpz_mod(t, t, p);
        mpz_mul(R, R, b);
        mpz_mod(R, R, p);
        mpz_set_ui(M, i);
    }

    if (mpz_cmp_ui(t, 1) == 0) {
        mpz_set(root1, R);
        mpz_sub(root2, p, R); // root2 = p - R
    } else {
        mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
        return -1;
    }

    mpz_clears(Q, z, M, c, t, R, temp, b, exponent, NULL);
    return 0;
}


// on choppe kp et kq en resolvant l'eq du second degré
int solve_kp_kq(const mpz_t N, const mpz_t e, const mpz_t k,
                mpz_t kp, mpz_t kq)
{
    mpz_t A, B, sqrt1, sqrt2;
    mpz_inits(A, B, sqrt1, sqrt2, NULL);

    // A = [k(N-1) + 1 + e] / 2
    mpz_sub_ui(A, N, 1);
    mpz_mul(A, A, k);
    mpz_add_ui(A, A, 1);

    // FIX astuce: e est tjs impair, donc on l'ajoute pour rendre A pair avant la div
    if (mpz_odd_p(A)) {
        mpz_add(A, A, e);
    }

    mpz_divexact_ui(A, A, 2);
    mpz_mod(A, A, e);

    // B = A^2 + k
    mpz_mul(B, A, A);
    mpz_add(B, B, k);
    mpz_mod(B, B, e);

    // on recup les racines via tonelli shanks
    if (tonelli_shanks(B, e, sqrt1, sqrt2) != 0) {
        fprintf(stderr, "[INIT] oups: pas de soluce pour l'eq de kp/kq\n");
        mpz_clears(A, B, sqrt1, sqrt2, NULL);
        return -1;
    }

    // kp = sqrt + A mod e (et kq c'est l'autre rep)
    mpz_add(kp, sqrt1, A);
    mpz_mod(kp, kp, e);

    mpz_add(kq, sqrt2, A);
    mpz_mod(kq, kq, e);

    mpz_clears(A, B, sqrt1, sqrt2, NULL);
    return 0;
}


// fix des bits LSB. c'est deterministe pcq ed = 1 mod ...
void correct_lsb(mpz_t out, const mpz_t e, const mpz_t k_val,
                 int tau_val, bool is_d)
{
    mpz_t modulus;
    mpz_init(modulus);

    // decallage diff pour d par rapport a dp/dq
    int power = is_d ? (tau_val + 2) : (tau_val + 1);

    mpz_ui_pow_ui(modulus, 2, (unsigned long)power);
    mpz_invert(out, e, modulus); // e^-1 mod 2^power

    mpz_clear(modulus);
}


// Fct de debug pour verif si on a bon en trichant avec la clef d'origine
bool verify_init_result(const rsa_key_t *original, const init_result_t *result) {
    bool ok = true;
    mpz_t tmp, phi, p1, q1;
    mpz_inits(tmp, phi, p1, q1, NULL);

    // verif k
    mpz_sub_ui(p1, original->p, 1);
    mpz_sub_ui(q1, original->q, 1);
    mpz_mul(phi, p1, q1);

    mpz_mul(tmp, original->e, original->d);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, phi);

    if (mpz_cmp(tmp, result->k) != 0) {
        gmp_printf("[VERIFY] pb sur k. attendu: %Zd, eu: %Zd\n", tmp, result->k);
        ok = false;
    } else {
        printf("[VERIFY] k: ok\n");
    }

    // verif kp
    mpz_mul(tmp, original->e, original->dp);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, p1);

    // kp et kq peuvent etre inversés, c normal
    if (mpz_cmp(tmp, result->kp) == 0) {
        printf("[VERIFY] kp: ok\n");
    } else if (mpz_cmp(tmp, result->kq) == 0) {
        printf("[VERIFY] kp: ok (mais inversé avec kq)\n");
    } else {
        ok = false;
    }

    // pareil pour kq
    mpz_mul(tmp, original->e, original->dq);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_divexact(tmp, tmp, q1);

    if (mpz_cmp(tmp, result->kp) == 0 || mpz_cmp(tmp, result->kq) == 0) {
        printf("[VERIFY] kq: ok\n");
    } else {
        ok = false;
    }

    // check rapide des tau()
    int expected_tau_k = tau(result->k);
    if (result->tau_k != expected_tau_k) {
        printf("[VERIFY] nope: τ(k) = %d au lieu de %d\n", result->tau_k, expected_tau_k);
        ok = false;
    }

    mpz_clears(tmp, phi, p1, q1, NULL);
    return ok;
}


// MAIN de la phase d'init. c'est ça qu'on appelle depuis le prg principal
int run_init_phase(const degraded_key_t *dk, init_result_t *result, bool verbose) {

    printf("╔══════════════════════════════════════╗\n");
    printf("║   Init Phase: recup de k, kp, kq     ║\n");
    printf("╚══════════════════════════════════════╝\n\n");

    // -- Etape 1 --
    printf("[INIT] Etape 1/4: recup k...\n");

    mpz_t d_tilde;
    mpz_init(d_tilde);

    int match_count = find_k(dk->key.N, dk->key.e, dk->key.d,
                             result->k, d_tilde, verbose);

    if (match_count == 0) {
        fprintf(stderr, "[INIT] erreur critique: impossible de trouver k\n");
        mpz_clear(d_tilde);
        return -1;
    }
    gmp_printf("[INIT] k = %Zd\n\n", result->k);

    // -- Etape 2 --
    printf("[INIT] Etape 2/4: fix des MSB de d...\n");

    correct_d_msb(d_tilde, dk->key.d, result->d_msb_corrected, dk->full_bits);
    printf("[INIT] MSB fixés\n\n");
    mpz_clear(d_tilde);

    // -- Etape 3 --
    printf("[INIT] Etape 3/4: resolution de l'eq pour kp et kq...\n");

    if (solve_kp_kq(dk->key.N, dk->key.e, result->k,
                    result->kp, result->kq) != 0) {
        fprintf(stderr, "[INIT] crash sur la resolution kp/kq\n");
        return -1;
    }
    gmp_printf("[INIT] kp = %Zd\n", result->kp);
    gmp_printf("[INIT] kq = %Zd\n\n", result->kq);

    // -- Etape 4 --
    printf("[INIT] Etape 4/4: calcul des tau et LSB...\n");

    result->tau_k  = tau(result->k);
    result->tau_kp = tau(result->kp);
    result->tau_kq = tau(result->kq);

    // on set p[0] et q[0] à 1 (un nombre premier est tjs impair)
    mpz_set_ui(result->p0, 0);
    mpz_set_ui(result->q0, 0);
    mpz_setbit(result->p0, 0);
    mpz_setbit(result->q0, 0);

    // fix deterministe des lsb
    correct_lsb(result->d0,  dk->key.e, result->k,  result->tau_k,  true);
    correct_lsb(result->dp0, dk->key.e, result->kp, result->tau_kp, false);
    correct_lsb(result->dq0, dk->key.e, result->kq, result->tau_kq, false);

    result->valid = true;
    printf("\n[INIT] c'est good l'init est fini.\n");

    if (verbose) {
        init_result_print(result);
    }

    return 0;
}