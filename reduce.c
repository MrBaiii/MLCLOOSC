#include "params.h"
#include "reduce.h"
#include <stdint.h>

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a with -2^{31}Q <= a <= Q*2^31,
*              compute r \equiv a*2^{-32} (mod Q) such that -Q < r < Q.
*
* Arguments:   - int64_t: finite field element a
*
* Returns r.
**************************************************/
int32_t montgomery_reduce(int64_t a) {
    int32_t t;
    t = (int32_t)a * QINV;
    t = (a - (int64_t)t * Q) >> 32;
    return t;
}
