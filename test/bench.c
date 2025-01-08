#include <stdint.h>
#include "../params.h"
#include "../polyvec.h"
#include "../poly.h"
#include "../randombytes.h"
#include "../symmetric.h"
#include "../fips202.h"
#include "speed_print.h"

#define NTESTS 1000
#define USE_SECOND

#ifdef uSE_CPUCYCLE

#include "cpucycles.h"

int main() {
    int i;
    uint64_t t[NTESTS];
    uint8_t seedbuf[2*SEEDBYTES];
    polyvecl mat[K];
    polyvecl s1, s2, s1hat;
    polyveck b;
    const uint8_t *rho, *rhoprime;

    randombytes(seedbuf, SEEDBYTES);
    shake256(seedbuf, 2*SEEDBYTES, seedbuf, SEEDBYTES);
    rho = seedbuf;
    rhoprime = seedbuf + SEEDBYTES;

    /* Expand matrix */
    polyvec_matrix_expand(mat, rho);
    /* Sample short vectors s1 and s2 */
    polyvecl_uniform_eta(&s1, rhoprime, 0);
    polyvecl_uniform_eta(&s2, rhoprime, L);

    for(i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        polyvecl_add(&s1hat, &s1, &s2);
    }
    print_results("polyvecl_add:", t, NTESTS);

    s1hat = s1;
    for(i = 0; i < NTESTS; ++i) {
        t[i] = cpucycles();
        polyvecl_ntt(&s1hat);
        polyvec_matrix_pointwise_montgomery(&b, mat, &s1hat);
        polyveck_invntt(&b);
    }
    print_results("polyvec_mul:", t, NTESTS);

    return 0;
}

#else

#include <time.h>
#include <stdio.h>

int main() {
    clock_t start, finish;

    int i;
    double t;
    uint8_t seedbuf[2*SEEDBYTES];
    polyvecl mat[K];
    polyvecl s1, s2, s1hat;
    polyveck b;
    const uint8_t *rho, *rhoprime;

    randombytes(seedbuf, SEEDBYTES);
    shake256(seedbuf, 2*SEEDBYTES, seedbuf, SEEDBYTES);
    rho = seedbuf;
    rhoprime = seedbuf + SEEDBYTES;

    /* Expand matrix */
    polyvec_matrix_expand(mat, rho);
    /* Sample short vectors s1 and s2 */
    polyvecl_uniform_eta(&s1, rhoprime, 0);
    polyvecl_uniform_eta(&s2, rhoprime, L);

    start = clock();
    for(i = 0; i < NTESTS; ++i) {
        polyvecl_add(&s1hat, &s1, &s2);
    }
    finish = clock();
    t = (double)(finish - start) * 1000 / CLOCKS_PER_SEC / NTESTS;
    printf("The average performance of polyvecl_add= %f milliseconds\n", t);

    s1hat = s1;
    start = clock();
    for(i = 0; i < NTESTS; ++i) {
        polyvecl_ntt(&s1hat);
        polyvec_matrix_pointwise_montgomery(&b, mat, &s1hat);
        polyveck_invntt(&b);
    }
    finish = clock();
    t = (double)(finish - start) * 1000 / CLOCKS_PER_SEC / NTESTS;
    printf("The average performance of polyvec_mul= %f milliseconds\n", t);

    return 0;
}

#endif