#ifndef REDUCE_H
#define REDUCE_H
#include "params.h"
#include <stdint.h>

#define MONT (84473) // 2^32 % Q
#define QINV 540932097 // q^(-1) mod 2^32

int32_t montgomery_reduce(int64_t a);

#endif
