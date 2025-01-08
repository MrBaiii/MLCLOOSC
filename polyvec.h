#ifndef POLYVEC_H
#define POLYVEC_H
#include "params.h"
#include "poly.h"
#include <stdint.h>

/* Vectors of polynomials of length L */
typedef struct {
    poly vec[L];
} polyvecl;

typedef struct {
    poly vec[K];
} polyveck;

void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);
void polyveck_uniform_eta(polyveck *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);

void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v);

void polyvecl_ntt(polyvecl *v);
void polyvecl_invntt(polyvecl *v);
void polyveck_invntt(polyveck *v);
void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a, const polyvecl *v);
void polyvecl_pointwise_acc_montgomery(poly *w,
        const polyvecl *u,
        const polyvecl *v);

void polyvec_matrix_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]);

void polyvec_matrix_pointwise_montgomery(polyveck *t, const polyvecl mat[K], const polyvecl *v);

#endif
