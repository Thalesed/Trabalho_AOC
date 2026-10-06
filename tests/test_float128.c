/* Independent decoder and binary-search quantizer: no tk_extract/from_log. */
#include "../takum.h"
#include <math.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned long checks;
static uint32_t state = 0x98765432;
static uint32_t rng(void) {
    state ^= state << 13; state ^= state >> 17; state ^= state << 5;
    return state;
}
static __float128 positive_log(uint32_t w, unsigned n) {
    uint32_t full = w << (32-n);
    unsigned d = (full >> 30) & 1, R = (full >> 27) & 7;
    unsigned r = d ? R : 7-R, p = 27-r;
    int bias = d ? (1 << r)-1 : -(1 << (r+1))+1;
    uint32_t payload = full & 0x07ffffff;
    return bias + (__float128)payload / (UINT32_C(1) << p);
}
static __float128 decode(uint32_t w, unsigned n) {
    uint32_t sign = UINT32_C(1) << (n-1), mask = UINT32_MAX >> (32-n);
    w &= mask;
    if (w == sign) return NAN;
    if (!w) return 0;
    int neg = (w & sign) != 0;
    if (neg) w = (0-w) & mask;
    return (neg ? -1 : 1) * expq(positive_log(w,n)/2);
}
static uint32_t encode(__float128 x, unsigned n) {
    uint32_t sign = UINT32_C(1) << (n-1), mask = UINT32_MAX >> (32-n);
    if (!isfinite(x)) return sign;
    if (!x) return 0;
    int neg = x < 0;
    __float128 l = 2*logq(fabsq(x));
    uint32_t lo = 1, hi = sign-1;
    while (lo < hi) {
        uint32_t m = lo+(hi-lo)/2;
        __float128 boundary = (positive_log(m,n)+positive_log(m+1,n))/2;
        if (l > boundary || (l == boundary && (m & 1))) lo = m+1;
        else hi = m;
    }
    return neg ? (0-lo) & mask : lo;
}
static void equal(uint32_t a, uint32_t b, const char *label, unsigned n,
                  uint32_t x, uint32_t y) {
    ++checks;
    if (a != b) {
        fprintf(stderr,"%s T%u x=%08x y=%08x ours=%08x ref=%08x\n",label,n,x,y,a,b);
        exit(EXIT_FAILURE);
    }
}
static void operations(uint32_t a, uint32_t b, unsigned n) {
    __float128 x = decode(a,n), y = decode(b,n);
    equal(tk_add(a,b,n),encode(x+y,n),"add128",n,a,b);
    equal(tk_sub(a,b,n),encode(x-y,n),"sub128",n,a,b);
    /* Product logarithms are exact dyadic rationals; exp/log in the oracle
       can perturb an exact tie, so multiplication uses exact log midpoints. */
}
int main(void) {
    for (unsigned n = 8; n <= 32; n *= 2) {
        uint32_t sign = UINT32_C(1) << (n-1);
        unsigned count = n < 32 ? sign-2 : 100000;
        for (unsigned i = 0; i < count; ++i) {
            uint32_t w = n < 32 ? i+1 : 1+rng()%(sign-2);
            __float128 mid = expq((positive_log(w,n)+positive_log(w+1,n))/4);
            long double center = (long double)mid;
            long double xs[] = {nextafterl(center,0),center,nextafterl(center,INFINITY)};
            for (unsigned j = 0; j < 3; ++j) {
                equal(tk_from_ld(xs[j],n),encode(xs[j],n),"boundary+",n,w,j);
                equal(tk_from_ld(-xs[j],n),encode(-xs[j],n),"boundary-",n,w,j);
            }
        }
        if (n == 8) {
            for (uint32_t a = 0; a < 256; ++a)
                for (uint32_t b = 0; b < 256; ++b) operations(a,b,n);
        } else {
            for (unsigned i = 0; i < 100000; ++i)
                operations(rng() & tk_mask(n),rng() & tk_mask(n),n);
        }
    }
    printf("OK: %lu comparacoes com referencia Float128 (fronteiras e soma/subtracao).\n",checks);
    return 0;
}
