#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include <stdint.h>

/*************************************************
* Name:        expand_mat
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|j|i)
*              or AES256CTR(rho,j|i).
*
* Arguments:   - polyvecl mat[K]: output matrix
*              - const uint8_t rho[]: byte array containing seed rho
**************************************************/
void polyvec_matrix_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]) {
    unsigned int i, j;

    for (i = 0; i < K; ++i) {
        for (j = 0; j < L; ++j) {
            poly_uniform(&mat[i].vec[j], rho, (uint16_t) ((i << 8) + j));
        }
    }
}

void polyvec_matrix_pointwise_montgomery(polyveck *t, const polyvecl mat[K], const polyvecl *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        polyvecl_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
    }
}

/**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/

void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
    unsigned int i;

    for (i = 0; i < L; ++i) {
        poly_uniform_eta(&v->vec[i], seed, nonce++);
    }
}

void polyveck_uniform_eta(polyveck *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        poly_uniform_eta(&v->vec[i], seed, nonce++);
    }
}

/*************************************************
* Name:        polyvecl_add
*
* Description: Add vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first summand
*              - const polyvecl *v: pointer to second summand
**************************************************/
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
    unsigned int i;

    for (i = 0; i < L; ++i) {
        poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
    }
}

/*************************************************
* Name:        polyvecl_ntt
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_ntt(polyvecl *v) {
    unsigned int i;

    for (i = 0; i < L; ++i) {
        poly_ntt(&v->vec[i]);
    }
}

void polyvecl_invntt(polyvecl *v) {
    unsigned int i;

    for (i = 0; i < L; ++i) {
        poly_invntt(&v->vec[i]);
    }
}

void polyveck_invntt(polyveck *v) {
    unsigned int i;

    for (i = 0; i < K; ++i) {
        poly_invntt(&v->vec[i]);
    }
}

void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a, const polyvecl *v) {
    unsigned int i;

    for (i = 0; i < L; ++i) {
        poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
    }
}

/*************************************************
* Name:        polyvecl_pointwise_acc_montgomery
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_pointwise_acc_montgomery(poly *w,
        const polyvecl *u,
        const polyvecl *v) {
    unsigned int i;
    poly t;
    poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
    for (i = 1; i < L; ++i) {
        poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
        poly_add(w, w, &t);
    }
}
