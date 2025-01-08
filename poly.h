#ifndef POLY_H
#define POLY_H
#include "params.h"
#include <stdint.h>

typedef struct {
    int32_t coeffs[N] __attribute__((aligned(16)));
} poly;

void poly_add(poly *c, const poly *a, const poly *b);

void poly_ntt(poly *a);
void poly_invntt(poly *a);
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

void poly_uniform(poly *a,
        const uint8_t seed[SEEDBYTES],
        uint16_t nonce);
void poly_uniform_eta(poly *a,
        const uint8_t seed[SEEDBYTES],
        uint16_t nonce);
#endif
