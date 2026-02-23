#include "common.h"
#include "keygen.h"
#include "decay_sim.h"
#include "init_phase.h"
#include "heninger_shacham.h"
#include "henecka_may.h"

#include <getopt.h>
#include <sys/stat.h>

/* ========================================================================
 * main.c — Recup de clef RSA post Cold-Boot Attack
 *
 * déroulement:
 * Phase 1: Gen de clef & Simu de la RAM qui freeze (decay)
 * Phase 2: Init — calcul des k, kp, kq depuis la clef pétée
 * Phase 3: Attaque Heninger-Shacham (HS09) -> rapide mais sensible au bruit
 * Phase 4: Attaque Henecka-May-Meurer (HMM10) -> stat thresholding
 *
 * flags importants pour les tests
 * ./rsa_recon --bits 1024 --preset cold --time 120
 * ./rsa_recon --bits 512 --preset room --time 8
 * ./rsa_recon --bits 1024 --simple --known 0.30 --error 0.0
 * ./rsa_recon --import keys/original_key.txt --preset frozen --time 3600
 * ======================================================================== */

// --- Helpers pour sortir des stats de debug ---

// compte cb de bits ont flip
static void compare_bits(const char *name, const mpz_t orig, const mpz_t deg, int nbits) {
    int differ = 0;
    for (int i = 0; i < nbits; i++) {
        if (get_bit(orig, i) != get_bit(deg, i))
            differ++;
    }
    printf("  %-4s: %d / %d bits differ (%.2f%%)\n",
           name, differ, nbits, 100.0 * differ / nbits);
}

// check si les bits marqués "connus" sont pas des mythos
static void check_known_accuracy(const char *name, const mpz_t orig, const int *known, int nbits) {
    int total_known = 0, correct = 0;
    for (int i = 0; i < nbits; i++) {
        if (known[i] != BIT_UNKNOWN) {
            total_known++;
            if (known[i] == get_bit(orig, i))
                correct++;
        }
    }
    if (total_known > 0) {
        printf("  %-4s: %d / %d known bits correct (%.4f%% accuracy)\n",
               name, correct, total_known, 100.0 * correct / total_known);
    }
}

// affichage de l'aide CLI
static void print_usage(const char *prog) {
    printf("Usage: %s [OPTIONS]\n\n", prog);
    printf("Key generation:\n");
    printf("  --bits N         Key size in bits (default: 1024)\n");
    printf("  --exp E          Public exponent (default: 65537)\n");
    printf("  --import FILE    Import existing key instead of generating\n");
    printf("\nDecay simulation:\n");
    printf("  --preset P       Temperature preset: room, cold, frozen, custom\n");
    printf("  --time T         Time since power loss in seconds\n");
    printf("  --temp C         Temperature in Celsius (overrides preset)\n");
    printf("  --wrong-dir P    Wrong-direction flip probability (default: 0.001)\n");
    printf("  --block-size B   Ground state block size in bits (default: 256)\n");
    printf("  --seed S         RNG seed (default: current time)\n");
    printf("\nSimple decay (flat probability, no ground state model):\n");
    printf("  --simple         Use simple decay model\n");
    printf("  --known F        Fraction of bits known (0.0 - 1.0)\n");
    printf("  --error F        Fraction of bits with wrong values (0.0 - 0.5)\n");
    printf("\nOutput:\n");
    printf("  --outdir DIR     Output directory (default: keys/)\n");
    printf("  -v, --verbose    Verbose output\n");
    printf("  -h, --help       Show this help\n");
}

int main(int argc, char *argv[]) {
    // config par defaut
    int key_bits = 1024;
    unsigned long pub_exp = 65537;
    char *import_file = NULL;
    char *preset = "cold";
    char *outdir = "keys";
    bool use_simple = false;
    bool verbose = false;
    double custom_time = -1;
    double custom_temp = -999;
    double custom_wrong_dir = -1;
    double simple_known = 0.30;
    double simple_error = 0.0;
    int block_size = 256;
    unsigned long seed = 0;

    // parsing des args facon linux
    static struct option long_opts[] = {
        {"bits",       required_argument, NULL, 'b'},
        {"exp",        required_argument, NULL, 'E'},
        {"import",     required_argument, NULL, 'i'},
        {"preset",     required_argument, NULL, 'p'},
        {"time",       required_argument, NULL, 't'},
        {"temp",       required_argument, NULL, 'T'},
        {"wrong-dir",  required_argument, NULL, 'w'},
        {"block-size", required_argument, NULL, 'B'},
        {"seed",       required_argument, NULL, 's'},
        {"simple",     no_argument,       NULL, 'S'},
        {"known",      required_argument, NULL, 'k'},
        {"error",      required_argument, NULL, 'r'},
        {"outdir",     required_argument, NULL, 'o'},
        {"verbose",    no_argument,       NULL, 'v'},
        {"help",       no_argument,       NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "b:i:p:t:o:vh", long_opts, NULL)) != -1) {
        switch (opt) {
            case 'b': key_bits = atoi(optarg); break;
            case 'E': pub_exp = strtoul(optarg, NULL, 10); break;
            case 'i': import_file = optarg; break;
            case 'p': preset = optarg; break;
            case 't': custom_time = atof(optarg); break;
            case 'T': custom_temp = atof(optarg); break;
            case 'w': custom_wrong_dir = atof(optarg); break;
            case 'B': block_size = atoi(optarg); break;
            case 's': seed = strtoul(optarg, NULL, 10); break;
            case 'S': use_simple = true; break;
            case 'k': simple_known = atof(optarg); break;
            case 'r': simple_error = atof(optarg); break;
            case 'o': outdir = optarg; break;
            case 'v': verbose = true; break;
            case 'h': print_usage(argv[0]); return 0;
            default:  print_usage(argv[0]); return 1;
        }
    }

    mkdir(outdir, 0755);

    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║   RSA Key Reconstruction — Cold Boot Attack Simulator    ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    // --- Etape 1: La clef cible ---
    // soit on la charge depuis un fichier, soit on en pond une nouvelle toute belle
    rsa_key_t original;
    rsa_key_init(&original);

    if (import_file) {
        printf("[STEP 1] Import clef cible depuis %s\n", import_file);
        if (rsa_key_import(&original, import_file) != 0) {
            fprintf(stderr, "Erreur: echec import\n");
            return 1;
        }
        key_bits = original.bits;
    } else {
        printf("[STEP 1] Generation clef RSA %d-bit (e = %lu)\n", key_bits, pub_exp);
        if (rsa_keygen(&original, key_bits, pub_exp) != 0) {
            fprintf(stderr, "Erreur: pb pdt la generation de la clef\n");
            return 1;
        }
    }

    // ptit check de sante
    if (rsa_key_verify(&original)) {
        printf("[STEP 1] Verif clef: OK ✓\n");
    } else {
        fprintf(stderr, "[STEP 1] Verif clef: FAILED ✗\n");
        rsa_key_clear(&original);
        return 1;
    }

    if (verbose) {
        rsa_key_print_summary(&original, "Clef Originale");
    }

    // backup de la clef propre pour comparer apres
    char path[512];
    snprintf(path, sizeof(path), "%s/original_key.txt", outdir);
    rsa_key_export(&original, path);
    printf("[STEP 1] Clef originale sauvée ds %s\n\n", path);

    // --- Etape 2: On pète la RAM ---
    // application du decay physique (courbe logistique)
    degraded_key_t dk;
    degraded_key_init(&dk, key_bits);

    if (use_simple) {
        printf("[STEP 2] Decay mode simple (known=%.2f, error=%.4f)\n", simple_known, simple_error);
        apply_simple_decay(&original, &dk, simple_known, simple_error, seed);
    } else {
        decay_params_t params;
        decay_params_preset(&params, preset);

        // overrides si besoin
        if (custom_time >= 0)    params.time_seconds = custom_time;
        if (custom_temp > -998)  params.temperature_celsius = custom_temp;
        if (custom_wrong_dir >= 0) params.wrong_dir_prob = custom_wrong_dir;
        params.ground_state_block_size = block_size;
        params.seed = seed;

        printf("[STEP 2] Lancement simu physique DRAM\n");
        decay_params_print(&params);
        printf("\n");

        apply_decay(&original, &dk, &params);
    }

    printf("\n");
    degraded_key_print_stats(&dk);

    // on dump la clef petée (ce que l'attaquant recupere en vrai)
    printf("\n[STEP 2] Export clef petée vers %s/\n", outdir);
    degraded_key_export(&dk, outdir);

    if (verbose) {
        rsa_key_print_summary(&dk.key, "Clef Degradee");
    }

    // --- Etape 3: Verif des stats ---
    printf("\n[STEP 3] Check des degats (original vs degrade)\n");

    compare_bits("p",  original.p,  dk.key.p,  dk.half_bits);
    compare_bits("q",  original.q,  dk.key.q,  dk.half_bits);
    compare_bits("d",  original.d,  dk.key.d,  dk.full_bits);
    compare_bits("dp", original.dp, dk.key.dp, dk.half_bits);
    compare_bits("dq", original.dq, dk.key.dq, dk.half_bits);

    printf("\n[STEP 3] Check si nos bits 'connus' sont legits\n");

    check_known_accuracy("p",  original.p,  dk.known_p,  dk.half_bits);
    check_known_accuracy("q",  original.q,  dk.known_q,  dk.half_bits);
    check_known_accuracy("d",  original.d,  dk.known_d,  dk.full_bits);
    check_known_accuracy("dp", original.dp, dk.known_dp, dk.half_bits);
    check_known_accuracy("dq", original.dq, dk.known_dq, dk.half_bits);

    printf("\n[DONE] Phase 1 & 2 terminees. Fichiers ds %s/\n", outdir);
    printf("  original_key.txt    — Clef d'origine\n");
    printf("  degraded_key.txt    — Clef post-freeze\n");
    printf("  known_bits_*.txt    — Arrays des bits qui ont survecus\n");
    printf("  decay_stats.txt     — Stats du decay\n");

    // --- Etape 4: Init ---
    // faut recup les variables mathematiques avant de lancer la reconstruction
    printf("\n");
    init_result_t init_res;
    init_result_init(&init_res);

    if (run_init_phase(&dk, &init_res, verbose) != 0) {
        fprintf(stderr, "[STEP 4] Crash phase d'init\n");
        init_result_clear(&init_res);
        rsa_key_clear(&original);
        degraded_key_clear(&dk);
        return 1;
    }

    printf("\n[STEP 4] Verif de l'init avec la vrai clef (cheat de dev)...\n");
    if (verify_init_result(&original, &init_res)) {
        printf("[STEP 4] Verif init: TOUT EST BON ✓\n");
    } else {
        printf("[STEP 4] Verif init: Y'A DES PB ✗\n");
        printf("  (Normal si la memoire est trop détruite)\n");
    }

    printf("\n[DONE] Go pour la phase 3: Reconstruction.\n");

    // --- Etape 5: Tentative via Heninger-Shacham (HS09) ---
    // parfait si ya des bits effacés mais 0 erreurs.
    printf("\n");
    hs_result_t hs_res;
    hs_result_init(&hs_res);

    int hs_ret = run_heninger_shacham(&dk, &init_res, &hs_res,
                                       true,  /* try_swap (astuce kp/kq) */
                                       verbose);

    if (hs_ret == 0) {
        printf("\n[STEP 5] Verif clef recuperée par HS09...\n");

        bool p_ok = (mpz_cmp(hs_res.recovered_key.p, original.p) == 0 ||
                     mpz_cmp(hs_res.recovered_key.p, original.q) == 0);
        bool q_ok = (mpz_cmp(hs_res.recovered_key.q, original.p) == 0 ||
                     mpz_cmp(hs_res.recovered_key.q, original.q) == 0);
        bool d_ok = (mpz_cmp(hs_res.recovered_key.d, original.d) == 0);

        printf("  p match:  %s\n", p_ok ? "OK ✓" : "MISMATCH ✗");
        printf("  q match:  %s\n", q_ok ? "OK ✓" : "MISMATCH ✗");
        printf("  d match:  %s\n", d_ok ? "OK ✓" : "MISMATCH ✗");

        if (rsa_key_verify(&hs_res.recovered_key)) {
            printf("  Verif crypto complete: OK ✓\n");
        } else {
            printf("  Verif crypto complete: FAILED ✗\n");
        }

        char rpath[512];
        snprintf(rpath, sizeof(rpath), "%s/recovered_key.txt", outdir);
        rsa_key_export(&hs_res.recovered_key, rpath);
        printf("\n[STEP 5] Clef recupérée sauvée dans %s\n", rpath);

        if (verbose) {
            rsa_key_print_summary(&hs_res.recovered_key, "Clef Recup HS");
        }
    } else {
        printf("\n[STEP 5] HS a fail — Il y a surement des bits inversés par erreur thermique.\n");
    }

    // --- Etape 6: Tentative via Henecka-May-Meurer (HMM10) ---
    // si HS a fail, on sort l'artillerie statistique lourde pour absorber les erreurs.
    printf("\n");
    hmm_params_t hmm_params;
    hmm_params_default(&hmm_params);
    hmm_params_auto(&hmm_params, key_bits, dk.error_rate); // calcul auto du seuil

    hmm_result_t hmm_res;
    hmm_result_init(&hmm_res);

    int hmm_ret = run_henecka_may(&dk, &init_res, &hmm_params, &hmm_res,
                                   true,  /* try_swap */
                                   verbose);

    if (hmm_ret == 0) {
        printf("\n[STEP 6] Verif clef recuperee par HMM10...\n");

        bool p_ok = (mpz_cmp(hmm_res.recovered_key.p, original.p) == 0 ||
                     mpz_cmp(hmm_res.recovered_key.p, original.q) == 0);
        bool q_ok = (mpz_cmp(hmm_res.recovered_key.q, original.p) == 0 ||
                     mpz_cmp(hmm_res.recovered_key.q, original.q) == 0);

        printf("  p match:  %s\n", p_ok ? "OK ✓" : "MISMATCH ✗");
        printf("  q match:  %s\n", q_ok ? "OK ✓" : "MISMATCH ✗");

        if (rsa_key_verify(&hmm_res.recovered_key)) {
            printf("  Verif crypto complete: OK ✓\n");
        } else {
            printf("  Verif crypto complete: FAILED ✗\n");
        }

        char hpath[512];
        snprintf(hpath, sizeof(hpath), "%s/recovered_key_hmm.txt", outdir);
        rsa_key_export(&hmm_res.recovered_key, hpath);
        printf("\n[STEP 6] Clef recupérée HMM sauvée ds %s\n", hpath);
    } else {
        printf("\n[STEP 6] HMM a fail aussi — la memoire est surement trop detruite.\n");
    }

    // --- Fin: Recap global ---
    printf("\n╔══════════════════════════════════════════════════╗\n");
    printf("║                  Bilan du run                    ║\n");
    printf("╚══════════════════════════════════════════════════╝\n");
    printf("  Taille clef:           %d bits\n", key_bits);
    printf("  Fraction de bits ok:   %.4f (%.1f%%)\n",
           dk.fraction_known, dk.fraction_known * 100.0);
    printf("  Tx d'erreur:           %.4f (%.1f%%)\n",
           dk.error_rate, dk.error_rate * 100.0);

    printf("\n  Heninger-Shacham (algo pur effacements):\n");
    printf("    Resultat:            %s\n", hs_ret == 0 ? "SUCCESS" : "FAILED");
    if (hs_ret == 0) {
        printf("    Temps d'exec:        %.3f s\n", hs_res.elapsed_seconds);
        printf("    Branches checkees:   %lld\n", hs_res.branches_explored);
    }

    printf("\n  Henecka-May-Meurer (algo stats tolerances erreurs):\n");
    printf("    Resultat:            %s\n", hmm_ret == 0 ? "SUCCESS" : "FAILED");
    if (hmm_ret == 0) {
        printf("    Temps d'exec:        %.3f s\n", hmm_res.elapsed_seconds);
        printf("    Canditats testes:    %lld\n", hmm_res.total_candidates_generated);
    }
    printf("    Taille bloc (t):     %d\n", hmm_params.block_size_t);
    printf("    Seuil (C):           %d\n", hmm_params.threshold_C);

    // menage final
    hmm_result_clear(&hmm_res);
    hs_result_clear(&hs_res);
    init_result_clear(&init_res);
    rsa_key_clear(&original);
    degraded_key_clear(&dk);

    return 0;
}