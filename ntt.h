#ifndef NTT_H
#define NTT_H
#include "params.h"
#include <stdint.h>

void ntt(int32_t a[N]);

void invntt(int32_t a[N]);

#endif
