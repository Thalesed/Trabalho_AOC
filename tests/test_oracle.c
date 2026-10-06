/* Validacao diferencial opcional contra libtakum, variante takum_log*. */
#include "../takum.h"
#include LIBTAKUM_HEADER
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long checks;
static uint32_t state = UINT32_C(0x27182818);
static uint32_t rng(void) { state ^= state << 13; state ^= state >> 17; state ^= state << 5; return state; }
static takum_log8 t8(uint32_t a) { uint8_t x = (uint8_t)a; takum_log8 t; memcpy(&t,&x,1); return t; }
static takum_log16 t16(uint32_t a) { uint16_t x = (uint16_t)a; takum_log16 t; memcpy(&t,&x,2); return t; }
static takum_log32 t32(uint32_t a) { takum_log32 t; memcpy(&t,&a,4); return t; }
static void equal(uint32_t ours, uint32_t ref, unsigned n, const char *label, uint32_t a, uint32_t b) {
    ++checks; ref &= tk_mask(n);
    if (ours != ref) {
        fprintf(stderr,"%s T%u a=%08x b=%08x ours=%08x ref=%08x\n",label,n,a,b,ours,ref);
        exit(EXIT_FAILURE);
    }
}
static double reference_decode(uint32_t a, unsigned n) {
    return n == 8 ? takum_log8_to_float64(t8(a)) : n == 16 ?
           takum_log16_to_float64(t16(a)) : takum_log32_to_float64(t32(a));
}
static uint32_t reference_encode(double x, unsigned n) {
    return n == 8 ? (uint8_t)takum_log8_from_float64(x) : n == 16 ?
           (uint16_t)takum_log16_from_float64(x) : (uint32_t)takum_log32_from_float64(x);
}
static uint32_t reference_op(uint32_t a, uint32_t b, unsigned n, unsigned op) {
    if (n == 8) {
        takum_log8 x = t8(a), y = t8(b);
        return (uint8_t)(op == 0 ? takum_log8_addition(x,y) : op == 1 ?
                         takum_log8_subtraction(x,y) : takum_log8_multiplication(x,y));
    }
    if (n == 16) {
        takum_log16 x = t16(a), y = t16(b);
        return (uint16_t)(op == 0 ? takum_log16_addition(x,y) : op == 1 ?
                          takum_log16_subtraction(x,y) : takum_log16_multiplication(x,y));
    }
    takum_log32 x = t32(a), y = t32(b);
    return (uint32_t)(op == 0 ? takum_log32_addition(x,y) : op == 1 ?
                      takum_log32_subtraction(x,y) : takum_log32_multiplication(x,y));
}
static void compare_operations(uint32_t a, uint32_t b, unsigned n) {
    equal(tk_add(a,b,n),reference_op(a,b,n,0),n,"add",a,b);
    equal(tk_sub(a,b,n),reference_op(a,b,n,1),n,"sub",a,b);
    equal(tk_mul(a,b,n),reference_op(a,b,n,2),n,"mul",a,b);
}
static void compare_decode(uint32_t word, unsigned n) {
    long double x = tk_to_ld(word,n);
    double y = reference_decode(word,n);
    ++checks;
    if (!((isnan(x) && isnan(y)) || x == y ||
          (y != 0 && fabsl((x-y)/y) < 3e-14L))) {
        fprintf(stderr,"decode T%u word=%08x x=%.20Le ref=%.20e rel=%.3Le\n",n,word,x,y,fabsl((x-y)/y)); exit(EXIT_FAILURE);
    }
}
int main(void) {
    for (unsigned n = 8; n <= 32; n *= 2) {
        uint32_t count = n < 32 ? 1u << n : 100000;
        for (uint32_t i = 0; i < count; ++i) {
            uint32_t word = n < 32 ? i : rng();
            compare_decode(word,n);
            if (n > 8) {
                uint32_t ref = n == 16 ? (uint8_t)takum_log8_from_takum_log16(t16(word)) :
                                        (uint16_t)takum_log16_from_takum_log32(t32(word));
                equal(tk_convert(word,n,n/2),ref,n/2,"narrow",word,0);
            }
        }
        const double cases[] = {0,-0.0,1,-1,INFINITY,-INFINITY,NAN,1e-300,-1e-300,1e300,-1e300};
        for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i)
            equal(tk_from_ld(cases[i],n),reference_encode(cases[i],n),n,"special",0,0);
        for (unsigned i = 0; i < 100000; ++i) {
            double l = ((double)rng()/UINT32_MAX)*600.0-300.0;
            double x = exp(l/2);
            if (rng() & 1) x = -x;
            equal(tk_from_ld(x,n),reference_encode(x,n),n,"encode",i,0);
        }
        if (n == 8) {
            for (uint32_t a = 0; a < 256; ++a)
                for (uint32_t b = 0; b < 256; ++b) compare_operations(a,b,n);
        } else {
            for (unsigned i = 0; i < 100000; ++i) {
                uint32_t a = rng() & tk_mask(n), b = rng() & tk_mask(n);
                compare_operations(a,b,n);
            }
        }
    }
    printf("OK: %lu comparacoes com libtakum (takum_log8/16/32).\n",checks);
    return EXIT_SUCCESS;
}
