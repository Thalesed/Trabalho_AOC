#ifndef AOC_TAKUM_H
#define AOC_TAKUM_H
#include <stdbool.h>
#include <stdint.h>

/* Palavras brutas sem sinal; os bits representam Takum LOGARITMICO. */
typedef uint8_t takum8_t;
typedef uint16_t takum16_t;
typedef uint32_t takum32_t;

typedef struct {
    unsigned S, D, R, r;
    unsigned characteristic_bits, fraction_bits;
    uint32_t C, F; /* Campos armazenados; C ausente em T8 e completado com zeros. */
    int c;
    long double m, l;
    bool zero, nar;
} takum_fields;

bool tk_width_valid(unsigned bits);
uint32_t tk_mask(unsigned bits);
uint32_t tk_nar(unsigned bits);
bool tk_extract(uint32_t word, unsigned bits, takum_fields *out);
long double tk_to_ld(uint32_t word, unsigned bits);
uint32_t tk_from_ld(long double value, unsigned bits);
long double tk_to_log(uint32_t word, unsigned bits);
uint32_t tk_from_log(bool negative, long double l, unsigned bits);
uint32_t tk_add(uint32_t a, uint32_t b, unsigned bits);
uint32_t tk_sub(uint32_t a, uint32_t b, unsigned bits);
uint32_t tk_mul(uint32_t a, uint32_t b, unsigned bits);
uint32_t tk_convert(uint32_t word, unsigned from, unsigned to);

/* API T16 mantida para compatibilidade com o código existente. */
float takum_to_float(takum16_t num);
takum16_t float_to_takum(float x);
takum16_t takum_add(takum16_t a, takum16_t b);
takum16_t takum_sub(takum16_t a, takum16_t b);
takum16_t takum_multiply(takum16_t a, takum16_t b);
void printBinary(int num);
#endif
