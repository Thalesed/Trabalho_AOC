/* Codec Takum logarítmico T8/T16/T32.
 * x = (-1)^S exp(l/2), l = (-1)^S (c+m). */
#include "takum.h"
#include "src/tk_wide.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Indice formado pelos bits DR armazenados. */
static const int C_BIAS[16] = {
    -255, -127, -63, -31, -15, -7, -3, -1,
    0, 1, 3, 7, 15, 31, 63, 127
};
bool tk_width_valid(unsigned n) { return n == 8 || n == 16 || n == 32; }
uint32_t tk_mask(unsigned n) {
    return tk_width_valid(n) ? (UINT32_MAX >> (32 - n)) : 0;
}
uint32_t tk_nar(unsigned n) {
    return tk_width_valid(n) ? UINT32_C(1) << (n - 1) : 0;
}
bool tk_extract(uint32_t word, unsigned n, takum_fields *f) {
    if (!tk_width_valid(n) || f == NULL) return false;
    word &= tk_mask(n);
    *f = (takum_fields){0};
    f->S = word >> (n - 1);
    f->D = (word >> (n - 2)) & 1;
    f->R = (word >> (n - 5)) & 7;
    f->r = f->D ? f->R : 7 - f->R;
    f->zero = word == 0;
    f->nar = word == tk_nar(n);
    unsigned available = n - 5;
    f->characteristic_bits = f->r < available ? f->r : available;
    f->fraction_bits = available - f->characteristic_bits;
    uint32_t low = word & ((UINT32_C(1) << available) - 1);
    f->C = low >> f->fraction_bits;
    f->F = low & ((UINT32_C(1) << f->fraction_bits) - 1);
    /* T8: completar os bits ausentes da caracteristica com zeros. */
    unsigned missing = f->r - f->characteristic_bits;
    f->c = C_BIAS[8 * f->D + f->R] + (int)(f->C << missing);
    f->m = ldexpl((long double)f->F, -(int)f->fraction_bits);
    f->l = f->zero ? -INFINITY : f->nar ? NAN :
           (f->S ? -1.0L : 1.0L) * (f->c + f->m);
    return true;
}
long double tk_to_log(uint32_t word, unsigned n) {
    takum_fields f;
    return tk_extract(word, n, &f) ? f.l : NAN;
}
long double tk_to_ld(uint32_t word, unsigned n) {
    takum_fields f;
    if (!tk_extract(word, n, &f) || f.nar) return NAN;
    if (f.zero) return 0.0L;
    return (f.S ? -1.0L : 1.0L) * expl(f.l / 2.0L);
}
/* Ao mais proximo; empate para palavra par. */
static uint64_t round_even(uint64_t base, tk_wide x) {
    tk_wide integral = tk_floor(x), fraction = x - integral;
    uint64_t q = base + (uint64_t)integral;
    return q + (fraction > 0.5L || (fraction == 0.5L && (q & 1)));
}
static uint32_t from_log_wide(bool negative, tk_wide l, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    uint32_t sign = tk_nar(n), mask = tk_mask(n);
    if (isnan(l) || (isinf(l) && l > 0)) return sign;
    if (isinf(l)) return 0;
    const long double bound = 254.99999904632568359375L;
    if (l > bound) l = bound;
    if (l < -bound) l = -bound;
    tk_wide cpm = negative ? -l : l;
    int c = (int)tk_floor(cpm);
    tk_wide m = cpm - c;
    unsigned dr = c >= 0 ? 8 : 0;
    while (dr < 15 && c >= C_BIAS[dr + 1]) ++dr;
    unsigned r = dr >= 8 ? dr - 8 : 7 - dr;
    unsigned p32 = 27 - r;
    uint32_t prefix = ((uint32_t)negative << 31) | (dr << 27) |
                      ((uint32_t)(c - C_BIAS[dr]) << p32);
    /* Separar parte inteira preserva a fracao em long double de 53 bits. */
    unsigned shift = 32 - n;
    uint32_t remainder = shift ? prefix & ((UINT32_C(1) << shift) - 1) : 0;
    tk_wide tail = tk_ldexp((tk_wide)remainder + tk_ldexp(m, (int)p32),
                             -(int)shift);
    uint64_t raw = round_even(prefix >> shift, tail);
    uint32_t result = (uint32_t)raw & mask;
    /* Saturar: numero finito nao nulo nunca vira zero ou NaR. */
    if (result == 0) return negative ? mask : 1;
    if (result == sign) return negative ? sign + 1 : sign - 1;
    return result;
}
uint32_t tk_from_log(bool negative, long double l, unsigned n) {
    return from_log_wide(negative, l, n);
}
static uint32_t from_wide(tk_wide x, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    if (!isfinite(x)) return tk_nar(n);
    if (x == 0.0L) return 0;
    return from_log_wide(x < 0, 2 * tk_log(tk_abs(x)), n);
}
uint32_t tk_from_ld(long double x, unsigned n) { return from_wide(x,n); }
static tk_wide to_wide(uint32_t a, unsigned n) {
    takum_fields f;
    if (!tk_extract(a,n,&f) || f.nar) return NAN;
    if (f.zero) return 0;
    return (f.S ? -1 : 1) * tk_exp((tk_wide)f.l / 2);
}
uint32_t tk_add(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) + to_wide(b,n), n);
}
uint32_t tk_sub(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) - to_wide(b,n), n);
}
uint32_t tk_mul(uint32_t a, uint32_t b, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    a &= tk_mask(n); b &= tk_mask(n);
    if (a == tk_nar(n) || b == tk_nar(n)) return tk_nar(n);
    if (a == 0 || b == 0) return 0;
    /* Multiplicacao logaritmica: soma dos logaritmos, sem exp/log extras. */
    return tk_from_log(((a ^ b) & tk_nar(n)) != 0,
                       tk_to_log(a,n) + tk_to_log(b,n), n);
}
uint32_t tk_convert(uint32_t word, unsigned from, unsigned to) {
    if (!tk_width_valid(from) || !tk_width_valid(to)) return 0;
    word &= tk_mask(from);
    if (to >= from) return word << (to - from);
    unsigned shift = from - to;
    uint32_t q = word >> shift;
    uint32_t rem = word & ((UINT32_C(1) << shift) - 1);
    uint32_t half = UINT32_C(1) << (shift - 1);
    q = (q + (rem > half || (rem == half && (q & 1)))) & tk_mask(to);
    bool negative = (word & tk_nar(from)) != 0;
    if (q == 0 && word != 0) q = negative ? tk_mask(to) : 1;
    if (q == tk_nar(to) && word != tk_nar(from))
        q = negative ? q + 1 : q - 1;
    return q;
}
float takum_to_float(takum16_t n) { return (float)tk_to_ld(n,16); }
takum16_t float_to_takum(float x) { return (takum16_t)tk_from_ld(x,16); }
takum16_t takum_add(takum16_t a, takum16_t b) { return (takum16_t)tk_add(a,b,16); }
takum16_t takum_sub(takum16_t a, takum16_t b) { return (takum16_t)tk_sub(a,b,16); }
takum16_t takum_multiply(takum16_t a, takum16_t b) { return (takum16_t)tk_mul(a,b,16); }
void printBinary(int num) {
    for (int i = 31; i >= 0; --i) putchar('0' + (((uint32_t)num >> i) & 1));
    putchar('\n');
}
#ifndef TAKUM_NO_MAIN
int main(void) {
    const long double values[] = {0,1,-1,2,-2,0.5L,3.14159L,100,1e-10L,1e10L};
    for (unsigned n = 8; n <= 32; n *= 2) {
        printf("\nTakum logaritmico T%u\n", n);
        for (size_t i = 0; i < sizeof(values)/sizeof(values[0]); ++i) {
            uint32_t t = tk_from_ld(values[i],n);
            printf("% .7Le -> 0x%0*X -> % .7Le\n", values[i], (int)n/4, t, tk_to_ld(t,n));
        }
        uint32_t a = tk_from_ld(2.5L,n), b = tk_from_ld(3,n);
        printf("Soma: %.7Lf; subtracao: %.7Lf; multiplicacao: %.7Lf\n",
               tk_to_ld(tk_add(a,b,n),n), tk_to_ld(tk_sub(a,b,n),n), tk_to_ld(tk_mul(a,b,n),n));
    }
    return 0;
}
#endif
/* Codec Takum logarítmico T8/T16/T32.
 * x = (-1)^S exp(l/2), l = (-1)^S (c+m). */
#include "takum.h"
#include "src/tk_wide.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Indice formado pelos bits DR armazenados. */
static const int C_BIAS[16] = {
    -255, -127, -63, -31, -15, -7, -3, -1,
    0, 1, 3, 7, 15, 31, 63, 127
};
bool tk_width_valid(unsigned n) { return n == 8 || n == 16 || n == 32; }
uint32_t tk_mask(unsigned n) {
    return tk_width_valid(n) ? (UINT32_MAX >> (32 - n)) : 0;
}
uint32_t tk_nar(unsigned n) {
    return tk_width_valid(n) ? UINT32_C(1) << (n - 1) : 0;
}
bool tk_extract(uint32_t word, unsigned n, takum_fields *f) {
    if (!tk_width_valid(n) || f == NULL) return false;
    word &= tk_mask(n);
    *f = (takum_fields){0};
    f->S = word >> (n - 1);
    f->D = (word >> (n - 2)) & 1;
    f->R = (word >> (n - 5)) & 7;
    f->r = f->D ? f->R : 7 - f->R;
    f->zero = word == 0;
    f->nar = word == tk_nar(n);
    unsigned available = n - 5;
    f->characteristic_bits = f->r < available ? f->r : available;
    f->fraction_bits = available - f->characteristic_bits;
    uint32_t low = word & ((UINT32_C(1) << available) - 1);
    f->C = low >> f->fraction_bits;
    f->F = low & ((UINT32_C(1) << f->fraction_bits) - 1);
    /* T8: completar os bits ausentes da caracteristica com zeros. */
    unsigned missing = f->r - f->characteristic_bits;
    f->c = C_BIAS[8 * f->D + f->R] + (int)(f->C << missing);
    f->m = ldexpl((long double)f->F, -(int)f->fraction_bits);
    f->l = f->zero ? -INFINITY : f->nar ? NAN :
           (f->S ? -1.0L : 1.0L) * (f->c + f->m);
    return true;
}
long double tk_to_log(uint32_t word, unsigned n) {
    takum_fields f;
    return tk_extract(word, n, &f) ? f.l : NAN;
}
long double tk_to_ld(uint32_t word, unsigned n) {
    takum_fields f;
    if (!tk_extract(word, n, &f) || f.nar) return NAN;
    if (f.zero) return 0.0L;
    return (f.S ? -1.0L : 1.0L) * expl(f.l / 2.0L);
}
/* Ao mais proximo; empate para palavra par. */
static uint64_t round_even(uint64_t base, tk_wide x) {
    tk_wide integral = tk_floor(x), fraction = x - integral;
    uint64_t q = base + (uint64_t)integral;
    return q + (fraction > 0.5L || (fraction == 0.5L && (q & 1)));
}
static uint32_t from_log_wide(bool negative, tk_wide l, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    uint32_t sign = tk_nar(n), mask = tk_mask(n);
    if (isnan(l) || (isinf(l) && l > 0)) return sign;
    if (isinf(l)) return 0;
    const long double bound = 254.99999904632568359375L;
    if (l > bound) l = bound;
    if (l < -bound) l = -bound;
    tk_wide cpm = negative ? -l : l;
    int c = (int)tk_floor(cpm);
    tk_wide m = cpm - c;
    unsigned dr = c >= 0 ? 8 : 0;
    while (dr < 15 && c >= C_BIAS[dr + 1]) ++dr;
    unsigned r = dr >= 8 ? dr - 8 : 7 - dr;
    unsigned p32 = 27 - r;
    uint32_t prefix = ((uint32_t)negative << 31) | (dr << 27) |
                      ((uint32_t)(c - C_BIAS[dr]) << p32);
    /* Separar parte inteira preserva a fracao em long double de 53 bits. */
    unsigned shift = 32 - n;
    uint32_t remainder = shift ? prefix & ((UINT32_C(1) << shift) - 1) : 0;
    tk_wide tail = tk_ldexp((tk_wide)remainder + tk_ldexp(m, (int)p32),
                             -(int)shift);
    uint64_t raw = round_even(prefix >> shift, tail);
    uint32_t result = (uint32_t)raw & mask;
    /* Saturar: numero finito nao nulo nunca vira zero ou NaR. */
    if (result == 0) return negative ? mask : 1;
    if (result == sign) return negative ? sign + 1 : sign - 1;
    return result;
}
uint32_t tk_from_log(bool negative, long double l, unsigned n) {
    return from_log_wide(negative, l, n);
}
static uint32_t from_wide(tk_wide x, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    if (!isfinite(x)) return tk_nar(n);
    if (x == 0.0L) return 0;
    return from_log_wide(x < 0, 2 * tk_log(tk_abs(x)), n);
}
uint32_t tk_from_ld(long double x, unsigned n) { return from_wide(x,n); }
static tk_wide to_wide(uint32_t a, unsigned n) {
    takum_fields f;
    if (!tk_extract(a,n,&f) || f.nar) return NAN;
    if (f.zero) return 0;
    return (f.S ? -1 : 1) * tk_exp((tk_wide)f.l / 2);
}
uint32_t tk_add(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) + to_wide(b,n), n);
}
uint32_t tk_sub(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) - to_wide(b,n), n);
}
uint32_t tk_mul(uint32_t a, uint32_t b, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    a &= tk_mask(n); b &= tk_mask(n);
    if (a == tk_nar(n) || b == tk_nar(n)) return tk_nar(n);
    if (a == 0 || b == 0) return 0;
    /* Multiplicacao logaritmica: soma dos logaritmos, sem exp/log extras. */
    return tk_from_log(((a ^ b) & tk_nar(n)) != 0,
                       tk_to_log(a,n) + tk_to_log(b,n), n);
}
uint32_t tk_convert(uint32_t word, unsigned from, unsigned to) {
    if (!tk_width_valid(from) || !tk_width_valid(to)) return 0;
    word &= tk_mask(from);
    if (to >= from) return word << (to - from);
    unsigned shift = from - to;
    uint32_t q = word >> shift;
    uint32_t rem = word & ((UINT32_C(1) << shift) - 1);
    uint32_t half = UINT32_C(1) << (shift - 1);
    q = (q + (rem > half || (rem == half && (q & 1)))) & tk_mask(to);
    bool negative = (word & tk_nar(from)) != 0;
    if (q == 0 && word != 0) q = negative ? tk_mask(to) : 1;
    if (q == tk_nar(to) && word != tk_nar(from))
        q = negative ? q + 1 : q - 1;
    return q;
}
float takum_to_float(takum16_t n) { return (float)tk_to_ld(n,16); }
takum16_t float_to_takum(float x) { return (takum16_t)tk_from_ld(x,16); }
takum16_t takum_add(takum16_t a, takum16_t b) { return (takum16_t)tk_add(a,b,16); }
takum16_t takum_sub(takum16_t a, takum16_t b) { return (takum16_t)tk_sub(a,b,16); }
takum16_t takum_multiply(takum16_t a, takum16_t b) { return (takum16_t)tk_mul(a,b,16); }
void printBinary(int num) {
    for (int i = 31; i >= 0; --i) putchar('0' + (((uint32_t)num >> i) & 1));
    putchar('\n');
}
#ifndef TAKUM_NO_MAIN
int main(void) {
    const long double values[] = {0,1,-1,2,-2,0.5L,3.14159L,100,1e-10L,1e10L};
    for (unsigned n = 8; n <= 32; n *= 2) {
        printf("\nTakum logaritmico T%u\n", n);
        for (size_t i = 0; i < sizeof(values)/sizeof(values[0]); ++i) {
            uint32_t t = tk_from_ld(values[i],n);
            printf("% .7Le -> 0x%0*X -> % .7Le\n", values[i], (int)n/4, t, tk_to_ld(t,n));
        }
        uint32_t a = tk_from_ld(2.5L,n), b = tk_from_ld(3,n);
        printf("Soma: %.7Lf; subtracao: %.7Lf; multiplicacao: %.7Lf\n",
               tk_to_ld(tk_add(a,b,n),n), tk_to_ld(tk_sub(a,b,n),n), tk_to_ld(tk_mul(a,b,n),n));
    }
    return 0;
}
#endif
/* Continuacao do prototipo T16 de Thales: codec comum T8/T16/T32.
 * Variante logaritmica: x = (-1)^S exp(l/2), l = (-1)^S (c+m).
 * Consulte docs/arquitetura.md para convencoes e limites do simulador. */
#include "takum.h"
#include "src/tk_wide.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Indice formado pelos bits DR armazenados. */
static const int C_BIAS[16] = {
    -255, -127, -63, -31, -15, -7, -3, -1,
    0, 1, 3, 7, 15, 31, 63, 127
};
bool tk_width_valid(unsigned n) { return n == 8 || n == 16 || n == 32; }
uint32_t tk_mask(unsigned n) {
    return tk_width_valid(n) ? (UINT32_MAX >> (32 - n)) : 0;
}
uint32_t tk_nar(unsigned n) {
    return tk_width_valid(n) ? UINT32_C(1) << (n - 1) : 0;
}
bool tk_extract(uint32_t word, unsigned n, takum_fields *f) {
    if (!tk_width_valid(n) || f == NULL) return false;
    word &= tk_mask(n);
    *f = (takum_fields){0};
    f->S = word >> (n - 1);
    f->D = (word >> (n - 2)) & 1;
    f->R = (word >> (n - 5)) & 7;
    f->r = f->D ? f->R : 7 - f->R;
    f->zero = word == 0;
    f->nar = word == tk_nar(n);
    unsigned available = n - 5;
    f->characteristic_bits = f->r < available ? f->r : available;
    f->fraction_bits = available - f->characteristic_bits;
    uint32_t low = word & ((UINT32_C(1) << available) - 1);
    f->C = low >> f->fraction_bits;
    f->F = low & ((UINT32_C(1) << f->fraction_bits) - 1);
    /* T8: completar os bits ausentes da caracteristica com zeros. */
    unsigned missing = f->r - f->characteristic_bits;
    f->c = C_BIAS[8 * f->D + f->R] + (int)(f->C << missing);
    f->m = ldexpl((long double)f->F, -(int)f->fraction_bits);
    f->l = f->zero ? -INFINITY : f->nar ? NAN :
           (f->S ? -1.0L : 1.0L) * (f->c + f->m);
    return true;
}
long double tk_to_log(uint32_t word, unsigned n) {
    takum_fields f;
    return tk_extract(word, n, &f) ? f.l : NAN;
}
long double tk_to_ld(uint32_t word, unsigned n) {
    takum_fields f;
    if (!tk_extract(word, n, &f) || f.nar) return NAN;
    if (f.zero) return 0.0L;
    return (f.S ? -1.0L : 1.0L) * expl(f.l / 2.0L);
}
/* Ao mais proximo; empate para palavra par. */
static uint64_t round_even(uint64_t base, tk_wide x) {
    tk_wide integral = tk_floor(x), fraction = x - integral;
    uint64_t q = base + (uint64_t)integral;
    return q + (fraction > 0.5L || (fraction == 0.5L && (q & 1)));
}
static uint32_t from_log_wide(bool negative, tk_wide l, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    uint32_t sign = tk_nar(n), mask = tk_mask(n);
    if (isnan(l) || (isinf(l) && l > 0)) return sign;
    if (isinf(l)) return 0;
    const long double bound = 254.99999904632568359375L;
    if (l > bound) l = bound;
    if (l < -bound) l = -bound;
    tk_wide cpm = negative ? -l : l;
    int c = (int)tk_floor(cpm);
    tk_wide m = cpm - c;
    unsigned dr = c >= 0 ? 8 : 0;
    while (dr < 15 && c >= C_BIAS[dr + 1]) ++dr;
    unsigned r = dr >= 8 ? dr - 8 : 7 - dr;
    unsigned p32 = 27 - r;
    uint32_t prefix = ((uint32_t)negative << 31) | (dr << 27) |
                      ((uint32_t)(c - C_BIAS[dr]) << p32);
    /* Separar parte inteira preserva a fracao em long double de 53 bits. */
    unsigned shift = 32 - n;
    uint32_t remainder = shift ? prefix & ((UINT32_C(1) << shift) - 1) : 0;
    tk_wide tail = tk_ldexp((tk_wide)remainder + tk_ldexp(m, (int)p32),
                             -(int)shift);
    uint64_t raw = round_even(prefix >> shift, tail);
    uint32_t result = (uint32_t)raw & mask;
    /* Saturar: numero finito nao nulo nunca vira zero ou NaR. */
    if (result == 0) return negative ? mask : 1;
    if (result == sign) return negative ? sign + 1 : sign - 1;
    return result;
}
uint32_t tk_from_log(bool negative, long double l, unsigned n) {
    return from_log_wide(negative, l, n);
}
static uint32_t from_wide(tk_wide x, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    if (!isfinite(x)) return tk_nar(n);
    if (x == 0.0L) return 0;
    return from_log_wide(x < 0, 2 * tk_log(tk_abs(x)), n);
}
uint32_t tk_from_ld(long double x, unsigned n) { return from_wide(x,n); }
static tk_wide to_wide(uint32_t a, unsigned n) {
    takum_fields f;
    if (!tk_extract(a,n,&f) || f.nar) return NAN;
    if (f.zero) return 0;
    return (f.S ? -1 : 1) * tk_exp((tk_wide)f.l / 2);
}
uint32_t tk_add(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) + to_wide(b,n), n);
}
uint32_t tk_sub(uint32_t a, uint32_t b, unsigned n) {
    return from_wide(to_wide(a,n) - to_wide(b,n), n);
}
uint32_t tk_mul(uint32_t a, uint32_t b, unsigned n) {
    if (!tk_width_valid(n)) return 0;
    a &= tk_mask(n); b &= tk_mask(n);
    if (a == tk_nar(n) || b == tk_nar(n)) return tk_nar(n);
    if (a == 0 || b == 0) return 0;
    /* Multiplicacao logaritmica: soma dos logaritmos, sem exp/log extras. */
    return tk_from_log(((a ^ b) & tk_nar(n)) != 0,
                       tk_to_log(a,n) + tk_to_log(b,n), n);
}
uint32_t tk_convert(uint32_t word, unsigned from, unsigned to) {
    if (!tk_width_valid(from) || !tk_width_valid(to)) return 0;
    word &= tk_mask(from);
    if (to >= from) return word << (to - from);
    unsigned shift = from - to;
    uint32_t q = word >> shift;
    uint32_t rem = word & ((UINT32_C(1) << shift) - 1);
    uint32_t half = UINT32_C(1) << (shift - 1);
    q = (q + (rem > half || (rem == half && (q & 1)))) & tk_mask(to);
    bool negative = (word & tk_nar(from)) != 0;
    if (q == 0 && word != 0) q = negative ? tk_mask(to) : 1;
    if (q == tk_nar(to) && word != tk_nar(from))
        q = negative ? q + 1 : q - 1;
    return q;
}
float takum_to_float(takum16_t n) { return (float)tk_to_ld(n,16); }
takum16_t float_to_takum(float x) { return (takum16_t)tk_from_ld(x,16); }
takum16_t takum_add(takum16_t a, takum16_t b) { return (takum16_t)tk_add(a,b,16); }
takum16_t takum_sub(takum16_t a, takum16_t b) { return (takum16_t)tk_sub(a,b,16); }
takum16_t takum_multiply(takum16_t a, takum16_t b) { return (takum16_t)tk_mul(a,b,16); }
void printBinary(int num) {
    for (int i = 31; i >= 0; --i) putchar('0' + (((uint32_t)num >> i) & 1));
    putchar('\n');
}
#ifndef TAKUM_NO_MAIN
int main(void) {
    const long double values[] = {0,1,-1,2,-2,0.5L,3.14159L,100,1e-10L,1e10L};
    for (unsigned n = 8; n <= 32; n *= 2) {
        printf("\nTakum logaritmico T%u\n", n);
        for (size_t i = 0; i < sizeof(values)/sizeof(values[0]); ++i) {
            uint32_t t = tk_from_ld(values[i],n);
            printf("% .7Le -> 0x%0*X -> % .7Le\n", values[i], (int)n/4, t, tk_to_ld(t,n));
        }
        uint32_t a = tk_from_ld(2.5L,n), b = tk_from_ld(3,n);
        printf("Soma: %.7Lf; subtracao: %.7Lf; multiplicacao: %.7Lf\n",
               tk_to_ld(tk_add(a,b,n),n), tk_to_ld(tk_sub(a,b,n),n), tk_to_ld(tk_mul(a,b,n),n));
    }
    return 0;
}
#endif
