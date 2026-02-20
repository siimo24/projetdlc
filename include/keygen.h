#ifndef KEYGEN_H
#define KEYGEN_H

#include "common.h"

/* ========================================================================
 *  keygen.h — RSA-CRT Key Generation via OpenSSL
 * ======================================================================== */

/*
 * Generate an RSA-CRT key pair using OpenSSL.
 * Extracts all components (N, e, d, p, q, dp, dq, qinv) into GMP format.
 *
 * @param key       Output: initialized rsa_key_t to fill
 * @param bits      Key size in bits (512, 1024, 2048, 4096)
 * @param pub_exp   Public exponent (typically 65537)
 * @return          0 on success, -1 on failure
 */
int rsa_keygen(rsa_key_t *key, int bits, unsigned long pub_exp);

/*
 * Generate a key and save it to files.
 *
 * @param bits      Key size
 * @param pub_exp   Public exponent
 * @param outdir    Output directory for key files
 * @return          0 on success
 */
int rsa_keygen_and_export(int bits, unsigned long pub_exp, const char *outdir);

#endif /* KEYGEN_H */
