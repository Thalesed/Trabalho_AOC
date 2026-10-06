#ifndef TK_WIDE_H
#define TK_WIDE_H
/* Float128 in the default GCC build; portable long double fallback. */
#ifdef TAKUM_USE_FLOAT128
#include <quadmath.h>
typedef __float128 tk_wide;
#define tk_exp expq
#define tk_log logq
#define tk_floor floorq
#define tk_ldexp ldexpq
#define tk_abs fabsq
#else
#include <math.h>
typedef long double tk_wide;
#define tk_exp expl
#define tk_log logl
#define tk_floor floorl
#define tk_ldexp ldexpl
#define tk_abs fabsl
#endif
#endif
