/* benchmark.c */
#if !defined(_WIN32)
#define _POSIX_C_SOURCE 200809L
#endif
#include "takum.h"
#include "src/rvv_takum.h"
#include "src/tk_wide.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#ifdef _WIN32
#include <windows.h>
#endif

/* ============================================================
   Configuração
   ============================================================ */
static const unsigned NUM_MATRICES   = 30;   /* matrizes DIFERENTES por tamanho */
static const unsigned MATRIX_SIZES[] = {64, 128, 256, 512};
#define SPARSITY   0.90
#define RNG_SEED   0xDEADBEEFu


/* Registrador vetorial de 128 bits (RVV padrão) → VL por formato */
#define RVV_VLEN_BITS (TK_VLEN_BYTES * 8u)   
static size_t vl_for_bits(unsigned elem_bits) {
    return RVV_VLEN_BITS / elem_bits;
}

/* Lista canônica de formatos (ordem usada nos CSVs) */
static const char *FORMAT_NAMES[] = {
    "E4M3_N", "E4M3_W",
    "E5M2_N", "E5M2_W",
    "F16_N",  "F16_W",
    "F32",
    "T8", "T16", "T32",
    "T16W"   /* Takum T16 acumulando em T-width (vtmul+vtadd) */
};
#define NUM_FORMATS (sizeof(FORMAT_NAMES) / sizeof(FORMAT_NAMES[0]))

/* ============================================================
   RNG determinístico (xorshift32)
   ============================================================ */
static uint32_t rng32_state(uint32_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    return *state;
}
static uint32_t rng32(uint32_t *s) { return rng32_state(s); }
static long double rng_ld(uint32_t *s) {
    return ((long double)rng32(s) / (long double)UINT32_MAX) * 2.0L - 1.0L;
}

/* ============================================================
   Matriz esparsa em CSR
   ============================================================ */
typedef struct {
    unsigned n;
    size_t nnz;
    uint32_t *row_ptr;
    uint32_t *col_idx;
    long double *vals;
} csr_matrix;

static void csr_free(csr_matrix *m) {
    free(m->row_ptr); free(m->col_idx); free(m->vals);
    memset(m, 0, sizeof(*m));
}

/*
 * Gera uma matriz CSR com seed própria. A mesma seed usada aqui é
 * gravada no CSV para reprodutibilidade.
 */
static int csr_generate(csr_matrix *m, unsigned n, uint32_t seed) {
    memset(m, 0, sizeof(*m));
    m->n = n;
    m->row_ptr = calloc((size_t)n + 1, sizeof(uint32_t));
    if (!m->row_ptr) return -1;

    uint32_t s_dec = seed ? seed : 1u;
    size_t total = 0;
    for (unsigned i = 0; i < n; ++i) {
        m->row_ptr[i] = (uint32_t)total;
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32(&s_dec) / (double)UINT32_MAX >= SPARSITY) ++total;
        }
    }
    m->row_ptr[n] = (uint32_t)total;
    m->nnz = total;
    m->col_idx = malloc((total ? total : 1) * sizeof(*m->col_idx));
    m->vals    = malloc((total ? total : 1) * sizeof(*m->vals));
    if ((!m->col_idx || !m->vals) && total > 0) { csr_free(m); return -1; }

    uint32_t s_dec2 = seed ? seed : 1u;
    uint32_t s_val  = seed ^ 0xA5A5A5A5u;
    if (s_val == 0) s_val = 1;
    size_t k = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32(&s_dec2) / (double)UINT32_MAX >= SPARSITY) {
                m->col_idx[k] = j;
                m->vals[k] = rng_ld(&s_val);
                ++k;
            }
        }
    }
    return 0;
}

/* ============================================================
   Referência Float128 (com flag) ou long double (fallback)
   ============================================================ */
#if defined(TAKUM_USE_FLOAT128) && !defined(_MSC_VER) && \
    (defined(__GNUC__) || defined(__clang__))
typedef __float128 ref_t;
#define REF(x) ((ref_t)(x))
#define HAVE_FLOAT128 1
#else
typedef long double ref_t;
#define REF(x) ((ref_t)(x))
#define HAVE_FLOAT128 0
#endif

static void spmv_ref(const csr_matrix *A, const long double *x, ref_t *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        ref_t s = 0;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k)
            s += REF(A->vals[k]) * REF(x[A->col_idx[k]]);
        y[i] = s;
    }
}

/* ============================================================
   Quantizadores (software) para os formatos IEEE
   ============================================================ */

/* ---- FP8 E4M3 (1-4-3, bias 7, max 448, sem Inf; NaN→max) ---- */
static float quantize_f8_e4m3(float val) {
    uint32_t b; memcpy(&b, &val, 4);
    uint32_t sign = (b >> 31) & 1u;
    int32_t  exp  = ((b >> 23) & 0xFF) - 127 + 7;
    uint32_t mant = b & 0x7FFFFFu;
    uint8_t  f8;

    if (isnan(val)) { f8 = 0x7Fu; } /* E4M3: 0x7F é NaN */
    else if (exp >= 15) { f8 = 0x7Eu; } /* clamp para max normal 448 */
    else if (exp <= 0) {
        if (exp < -9) f8 = 0;
        else {
            mant |= 0x800000u;
            int shift = 23 - 3 + (1 - exp);
            if (shift >= 32) f8 = 0;
            else {
                mant += (1u << (shift - 1));
                f8 = (uint8_t)((mant >> shift) & 0x7u); /* subnormal: exp=0 */
            }
        }
    } else {
        int shift = 23 - 3;
        mant += (1u << (shift - 1));
        if (mant & 0x800000u) { mant >>= 1; exp++; }
        if (exp >= 15) f8 = 0x7Eu;
        else f8 = (uint8_t)((exp << 3) | ((mant >> shift) & 0x7u));
    }
    f8 = (uint8_t)(f8 | (sign << 7));

    /* Reverter para float */
    uint32_t r_sign = (f8 >> 7) & 1u;
    uint32_t r_exp  = (f8 >> 3) & 0xFu;
    uint32_t r_mant = f8 & 0x7u;
    uint32_t r_bits = (r_sign << 31);
    if (r_exp == 0) {
        if (r_mant == 0) { r_bits = r_sign ? 0x80000000u : 0u; }
        else {
            int e = -6;
            while (!(r_mant & 0x8u)) { r_mant <<= 1; e--; }
            r_mant &= 0x7u;
            r_bits |= ((uint32_t)(e + 127) << 23) | (r_mant << 20);
        }
    } else if (r_exp == 0xFu && r_mant == 0x7u) {
        r_bits |= 0x7FC00000u; /* NaN */
    } else {
        r_bits |= ((uint32_t)((int)r_exp - 7 + 127) << 23) | (r_mant << 20);
    }
    float res; memcpy(&res, &r_bits, 4);
    return res;
}

/* ---- FP8 E5M2 (1-5-2, bias 15, max 57344; NaN/Inf como IEEE) ---- */
static float quantize_f8_e5m2(float val) {
    uint32_t b; memcpy(&b, &val, 4);
    uint32_t sign = (b >> 31) & 1u;
    int32_t  exp  = ((b >> 23) & 0xFF) - 127 + 15;
    uint32_t mant = b & 0x7FFFFFu;
    uint8_t  f8;

    if (isnan(val)) { f8 = 0x7Fu; }
    else if (isinf(val)) { f8 = 0x7Cu; }
    else if (exp >= 31) { f8 = 0x7Bu; } /* clamp para max normal 57344 */
    else if (exp <= 0) {
        if (exp < -16) f8 = 0;
        else {
            mant |= 0x800000u;
            int shift = 23 - 2 + (1 - exp);
            if (shift >= 32) f8 = 0;
            else {
                mant += (1u << (shift - 1));
                f8 = (uint8_t)((mant >> shift) & 0x3u);
            }
        }
    } else {
        int shift = 23 - 2;
        mant += (1u << (shift - 1));
        if (mant & 0x800000u) { mant >>= 1; exp++; }
        if (exp >= 31) f8 = 0x7Bu;
        else f8 = (uint8_t)((exp << 2) | ((mant >> shift) & 0x3u));
    }
    f8 = (uint8_t)(f8 | (sign << 7));

    uint32_t r_sign = (f8 >> 7) & 1u;
    uint32_t r_exp  = (f8 >> 2) & 0x1Fu;
    uint32_t r_mant = f8 & 0x3u;
    uint32_t r_bits = (r_sign << 31);
    if (r_exp == 0) {
        if (r_mant == 0) { r_bits = r_sign ? 0x80000000u : 0u; }
        else {
            int e = -14;
            while (!(r_mant & 0x4u)) { r_mant <<= 1; e--; }
            r_mant &= 0x3u;
            r_bits |= ((uint32_t)(e + 127) << 23) | (r_mant << 21);
        }
    } else if (r_exp == 0x1Fu) {
        r_bits |= 0x7F800000u | (r_mant << 21); /* Inf/NaN */
    } else {
        r_bits |= ((uint32_t)((int)r_exp - 15 + 127) << 23) | (r_mant << 21);
    }
    float res; memcpy(&res, &r_bits, 4);
    return res;
}


/* Quantização manual para half (RNE) */
static float q16(float val) {
    uint32_t b; memcpy(&b, &val, 4);
    uint32_t sign = (b >> 31) & 1u;
    int32_t  exp  = ((b >> 23) & 0xFF) - 127 + 15;
    uint32_t mant = b & 0x7FFFFFu;

    if (isnan(val)) { uint32_t r = (sign << 31) | 0x7E000000u; float f; memcpy(&f,&r,4); return f; }
    if (isinf(val)) { uint32_t r = (sign << 31) | 0x7C000000u; float f; memcpy(&f,&r,4); return f; }

    uint16_t h;
    if (exp >= 31) h = 0x7BFFu;                       /* max half */
    else if (exp <= 0) {
        if (exp < -10) h = 0;
        else {
            mant |= 0x800000u;
            int shift = 23 - 10 + (1 - exp);
            if (shift >= 32) h = 0;
            else {
                mant += (1u << (shift - 1));
                h = (uint16_t)((mant >> shift) & 0x3FFu);
            }
        }
    } else {
        int shift = 23 - 10;
        mant += (1u << (shift - 1));
        if (mant & 0x800000u) { mant >>= 1; exp++; }
        if (exp >= 31) h = 0x7BFFu;
        else h = (uint16_t)((exp << 10) | ((mant >> shift) & 0x3FFu));
    }
    h = (uint16_t)(h | (sign << 15));

    uint32_t r_sign = (h >> 15) & 1u;
    uint32_t r_exp  = (h >> 10) & 0x1Fu;
    uint32_t r_mant = h & 0x3FFu;
    uint32_t r_bits = (r_sign << 31);
    if (r_exp == 0) {
        if (r_mant == 0) { r_bits = r_sign ? 0x80000000u : 0u; }
        else {
            int e = -14;
            while (!(r_mant & 0x400u)) { r_mant <<= 1; e--; }
            r_mant &= 0x3FFu;
            r_bits |= ((uint32_t)(e + 127) << 23) | (r_mant << 13);
        }
    } else if (r_exp == 0x1Fu) {
        r_bits |= 0x7F800000u | (r_mant << 13);
    } else {
        r_bits |= ((uint32_t)((int)r_exp - 15 + 127) << 23) | (r_mant << 13);
    }
    float res; memcpy(&res, &r_bits, 4);
    return res;
}


/* ============================================================
   SpMV — IEEE 754 (Narrow e Wide)
   ============================================================ */

/* F32 puro */
static void spmv_f32(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k)
            s += (float)A->vals[k] * (float)x[A->col_idx[k]];
        y[i] = s;
    }
}

/* F16 narrow: produto quantizado em half, acumulador half (forçado por quantização a cada soma) */
static void spmv_f16_narrow(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        volatile float s = 0.f;   /* impede fusão q16(q16(...)) pelo -O2 */
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = q16((float)A->vals[k]);
            float b = q16((float)x[A->col_idx[k]]);
            float prod = q16(a * b);
            s = q16(prod + s);
        }
        y[i] = s;
    }
}

/* F16 wide: produto quantizado em half, acumulador FP32 */
static void spmv_f16_wide(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = q16((float)A->vals[k]);
            float b = q16((float)x[A->col_idx[k]]);
            s += a * b;
        }
        y[i] = s;
    }
}

/* FP8 E4M3 narrow: produto arredondado p/ E4M3, acumulador E4M3 (forçado a cada passo) */
static void spmv_e4m3_narrow(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e4m3((float)A->vals[k]);
            float b = quantize_f8_e4m3((float)x[A->col_idx[k]]);
            s = quantize_f8_e4m3(quantize_f8_e4m3(a * b) + s);
        }
        y[i] = s;
    }
}

/* FP8 E4M3 wide: produto arredondado p/ E4M3, acumulador FP32 */
static void spmv_e4m3_wide(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e4m3((float)A->vals[k]);
            float b = quantize_f8_e4m3((float)x[A->col_idx[k]]);
            s += a * b;
        }
        y[i] = s;
    }
}

/* FP8 E5M2 narrow */
static void spmv_e5m2_narrow(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e5m2((float)A->vals[k]);
            float b = quantize_f8_e5m2((float)x[A->col_idx[k]]);
            s = quantize_f8_e5m2(quantize_f8_e5m2(a * b) + s);
        }
        y[i] = s;
    }
}

/* FP8 E5M2 wide */
static void spmv_e5m2_wide(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e5m2((float)A->vals[k]);
            float b = quantize_f8_e5m2((float)x[A->col_idx[k]]);
            s += a * b;
        }
        y[i] = s;
    }
}

/* ============================================================
   SpMV — Takum
   ============================================================ */

/* Takum padrão: vtdot (acúmulo em precisão ampliada, conforme especificação) */
static int spmv_takum_vtdot(const csr_matrix *A, const long double *x,
                            unsigned sew, uint32_t *y_out) {
    tk_sim sim; tk_sim_init(&sim);
    size_t vlmax = tk_sim_vlmax(sew);
    if (vlmax == 0) return -1;
    for (unsigned i = 0; i < A->n; ++i) {
        uint32_t start = A->row_ptr[i], end = A->row_ptr[i + 1];
        uint32_t nnz = end - start;
        if (nnz == 0) { y_out[i] = tk_from_ld(0, sew); continue; }
        long double acc = 0;
        uint32_t remain = nnz, off = start;
        while (remain) {
            uint32_t chunk = remain < vlmax ? remain : (uint32_t)vlmax;
            long double a_buf[TK_VLEN_BYTES], b_buf[TK_VLEN_BYTES];
            for (uint32_t t = 0; t < chunk; ++t) {
                a_buf[t] = A->vals[off + t];
                b_buf[t] = x[A->col_idx[off + t]];
            }
            tk_sim_setvl(&sim, chunk, sew);
            tk_sim_load(&sim, 1, a_buf, chunk);
            tk_sim_load(&sim, 2, b_buf, chunk);
            uint32_t insn = tk_instruction(TK_VTDOT, 3, 1, 2);
            tk_sim_execute(&sim, insn);
            acc += sim.dot_accumulator;
            remain -= chunk; off += chunk;
        }
        y_out[i] = tk_from_ld(acc, sew);
    }
    return 0;
}

/*
 * Takum T-width: acumula em T-width via vtmul + vtadd, um lane por vez
 * (setvl=1). Isso espelha o IEEE narrow: cada operação intermediária é
 * arredondada para a largura do elemento. A leitura do resultado da
 * operação usa apenas a interface pública de tk_sim (dot_accumulator),
 * evitando depender de membros internos (sim.vreg).
 */
static int spmv_takum_twidth(const csr_matrix *A, const long double *x,
                             unsigned sew, uint32_t *y_out) {
    tk_sim sim; tk_sim_init(&sim);
    if (tk_sim_vlmax(sew) == 0) return -1;

    for (unsigned i = 0; i < A->n; ++i) {
        uint32_t start = A->row_ptr[i], end = A->row_ptr[i + 1];
        uint32_t nnz = end - start;
        if (nnz == 0) { y_out[i] = tk_from_ld(0, sew); continue; }

        uint32_t acc_t = tk_from_ld(0.0L, sew);   /* acumulador em T-width */

        for (uint32_t k = start; k < end; ++k) {
            long double a = A->vals[k];
            long double b = x[A->col_idx[k]];

            /* vtmul: 1 lane, v3 <- v1 * v2 (arredonda para T-width) */
            tk_sim_setvl(&sim, 1, sew);
            tk_sim_load(&sim, 1, &a, 1);
            tk_sim_load(&sim, 2, &b, 1);
                        uint32_t insn_mul = tk_instruction(TK_VTMUL, 3, 1, 2);
            tk_sim_execute(&sim, insn_mul);

            /* Lê o produto de v[3] via API pública (VTADD/VTMUL não tocam dot_accumulator) */
            uint32_t prod_word = 0;
            tk_sim_read(&sim, 3, sew, 0, &prod_word);
            long double prod = tk_to_ld(prod_word, sew);

            /* vtadd: acc_t + prod, um lane */
            long double acc_ld = tk_to_ld(acc_t, sew);
            long double sum_ld = acc_ld + prod;
            tk_sim_setvl(&sim, 1, sew);
            tk_sim_load(&sim, 1, &sum_ld, 1);
            tk_sim_load(&sim, 2, &sum_ld, 1);   /* placeholder; sobrescrito abaixo */
            /* Carrega acc_ld e prod_ld nas lanes 1 e 2, não sum_ld duas vezes */
            tk_sim_load(&sim, 1, &acc_ld, 1);
            tk_sim_load(&sim, 2, &prod, 1);
            uint32_t insn_add = tk_instruction(TK_VTADD, 3, 1, 2);
            tk_sim_execute(&sim, insn_add);

            tk_sim_read(&sim, 3, sew, 0, &acc_t);
        }
        y_out[i] = acc_t;
    }
    return 0;
}

/* ============================================================
   Estatísticas de erro
   ============================================================ */
typedef struct {
    double mean, rms, max, p95, p99;
    double normwise;
    size_t count;
} err_stats;

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    if (isnan(x)) return isnan(y) ? 0 : 1;
    if (isnan(y)) return -1;
    return (x > y) - (x < y);
}

/*
 * Calcula estatísticas para um formato. A CDF por (size, format) é
 * preenchida em cdf_buf/idx se fornecidos.
 */
static err_stats compute_err(const ref_t *ref, const float *est, unsigned n,
                             double *cdf_buf, size_t *cdf_idx) {
    err_stats s = {0};
    if (n == 0) return s;
    double *errors = malloc((size_t)n * sizeof(*errors));
    if (!errors) return s;

    long double sum = 0, sumsq = 0, num_norm = 0, den_norm = 0;
    for (unsigned i = 0; i < n; ++i) {
        long double r = (long double)ref[i];
        long double v = (long double)est[i];
        long double denom = fmaxl(fabsl(r), 1.0e-30L);
        long double rel = fabsl(v - r) / denom;
        if (!isfinite(rel) || rel > (long double)DBL_MAX) rel = (long double)DBL_MAX;
        errors[i] = (double)rel;
        if (cdf_buf && cdf_idx) cdf_buf[(*cdf_idx)++] = (double)rel;
        sum += rel; sumsq += rel * rel;
        num_norm += (v - r) * (v - r);
        den_norm += r * r;
    }
    qsort(errors, n, sizeof(*errors), cmp_double);
    s.count = n;
    s.mean = (double)(sum / (long double)n);
    s.rms  = (double)sqrtl(sumsq / (long double)n);
    s.max  = errors[n - 1];
    s.p95  = errors[(size_t)ceill(0.95L * n) - 1];
    s.p99  = errors[(size_t)ceill(0.99L * n) - 1];
    s.normwise = (den_norm > 0) ? (double)sqrtl(num_norm / den_norm) : 0.0;
    free(errors);
    return s;
}

/* ============================================================
   Modelo de instruções RVV (unificado)
   ============================================================
   Por bloco (chunk) de no máximo VL elementos:
     IEEE narrow/wide : vsetvli + vle(vals) + vluxei(gather x) + vfmul + vfredusum = 5
     Takum  vtdot     : vsetvli + vle(vals) + vluxei(gather x) + vtdot            = 4
     Takum  twidth    : vsetvli + vle(vals) + vluxei(gather x) + vtmul + vtadd    = 5
*/
static uint64_t rvv_chunks(size_t nnz, size_t vl) {
    if (vl == 0) return 0;
    return (nnz + vl - 1) / vl;
}
static uint64_t rvv_inst_est(size_t nnz, size_t vl, int is_takum_vtdot, int is_takum_twidth) {
    uint64_t c = rvv_chunks(nnz, vl);
    if (c == 0) return 0;
    uint64_t per = 5;
    if (is_takum_vtdot)   per = 4;
    if (is_takum_twidth)  per = 5;
    return c * per;
}

/* ============================================================
   Agregação estatística por (size, format)
   ============================================================ */
typedef struct {
    double *vals;   /* amostras */
    size_t  n, cap;
} sample_vec;

static void sv_push(sample_vec *sv, double v) {
    if (sv->n == sv->cap) {
        size_t nc = sv->cap ? sv->cap * 2 : 16;
        double *nv = realloc(sv->vals, nc * sizeof(double));
        if (!nv) return;
        sv->vals = nv; sv->cap = nc;
    }
    sv->vals[sv->n++] = v;
}
static void sv_free(sample_vec *sv) { free(sv->vals); sv->vals = NULL; sv->n = sv->cap = 0; }

/* média, desvio amostral, IC95 (t≈1.96 para n≥20; suficiente aqui) */
static void sv_stats(const sample_vec *sv, double *avg, double *sd, double *ci) {
    if (sv->n == 0) { *avg = *sd = *ci = 0; return; }
    double m = 0;
    for (size_t i = 0; i < sv->n; ++i) m += sv->vals[i];
    m /= (double)sv->n;
    double var = 0;
    if (sv->n > 1) {
        for (size_t i = 0; i < sv->n; ++i) {
            double d = sv->vals[i] - m; var += d * d;
        }
        var /= (double)(sv->n - 1);
    }
    double s = sqrt(var);
    *avg = m; *sd = s; *ci = 1.96 * s / sqrt((double)sv->n);
}

/* ============================================================
   Utilitário de tempo
   ============================================================ */
static double now_sec(void) {
#ifdef _WIN32
    LARGE_INTEGER c, f;
    QueryPerformanceFrequency(&f); QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#endif
}

/* ============================================================
   main
   ============================================================ */
int main(void) {
#if !HAVE_FLOAT128
    printf("[AVISO] Compilado sem TAKUM_USE_FLOAT128. Referencia = long double (80-bit).\n");
    printf("[AVISO] Para o artigo, recompile com: -DTAKUM_USE_FLOAT128 -lquadmath\n\n");
#endif
    printf("[AVISO] Tempo medido = custo de SIMULACAO SOFTWARE (expl/logl no Takum).\n");
    printf("[AVISO] NAO representa latencia de hardware.\n\n");

    const size_t n_sizes = sizeof(MATRIX_SIZES) / sizeof(MATRIX_SIZES[0]);

    FILE *csv  = fopen("benchmark_results.csv", "w");
    FILE *stt  = fopen("benchmark_statistics.csv", "w");
    FILE *cdf  = fopen("benchmark_cdf.csv", "w");
    FILE *sum  = fopen("benchmark_summary.txt", "w");
    if (!csv || !stt || !cdf || !sum) { fprintf(stderr, "OOM\n"); return 1; }

    /* Cabeçalhos */
    fprintf(csv, "size,nnz,format,mat_id,seed,time_s,rvv_chunks,rvv_inst_est,"
                 "err_mean,err_rms,err_max,err_p95,err_p99,err_normwise\n");
    fprintf(stt, "size,format,count,"
                 "err_mean_avg,err_mean_sd,err_mean_ci95,"
                 "normwise_avg,normwise_sd,normwise_ci95,"
                 "time_ms_avg,time_ms_sd,time_ms_ci95\n");
    fprintf(cdf, "size,format,error_relative,cumulative_probability\n");

    fprintf(sum, "=== Benchmark SpMV (y = A*x) — Takum vs IEEE 754 ===\n");
    fprintf(sum, "NUM_MATRICES=%u (matrizes diferentes por tamanho)\n", NUM_MATRICES);
    fprintf(sum, "SPARSITY=%.2f  RNG_SEED=0x%X  RVV_VLEN=%u bits\n\n",
            SPARSITY, RNG_SEED, RVV_VLEN_BITS);

    fprintf(sum, "Formatos avaliados:\n");
    for (size_t i = 0; i < NUM_FORMATS; ++i) fprintf(sum, "  - %s\n", FORMAT_NAMES[i]);
    fprintf(sum, "\n");

    fprintf(sum, "=== Modelo de instrucoes RVV (por chunk) ===\n");
    fprintf(sum, "  IEEE (todos):  vsetvli + vle + vluxei + vfmul + vfredusum = 5\n");
    fprintf(sum, "  Takum vtdot:   vsetvli + vle + vluxei + vtdot              = 4\n");
    fprintf(sum, "  Takum twidth:  vsetvli + vle + vluxei + vtmul + vtadd      = 5\n");
    fprintf(sum, "  VL por formato (registrador de %u bits): ", RVV_VLEN_BITS);
    fprintf(sum, "E4M3/E5M2=%zu, F16/T16=%zu, F32/T32=%zu\n\n",
            vl_for_bits(8), vl_for_bits(16), vl_for_bits(32));

    fprintf(sum, "=== CSRs ===\n");
    fprintf(sum, "  IEEE 754 exige frm + fflags (2 CSRs) para arredondamento e excecoes.\n");
    fprintf(sum, "  Takum usa RNE fixo e nao gera excecoes (apenas NaR) => 0 CSRs.\n\n");

    for (size_t si = 0; si < n_sizes; ++si) {
        unsigned n = MATRIX_SIZES[si];
        printf("[benchmark] size=%u, %u matrizes...\n", n, NUM_MATRICES);

        /* Buffers de agregação por formato (indexado por FORMAT_NAMES) */
        sample_vec sv_err_mean[NUM_FORMATS];
        sample_vec sv_normwise[NUM_FORMATS];
        sample_vec sv_time[NUM_FORMATS];
        memset(sv_err_mean, 0, sizeof(sv_err_mean));
        memset(sv_normwise, 0, sizeof(sv_normwise));
        memset(sv_time,     0, sizeof(sv_time));

        /* CDF por formato para este size */
        double *cdf_buf[NUM_FORMATS];
        size_t  cdf_idx[NUM_FORMATS];
        for (size_t f = 0; f < NUM_FORMATS; ++f) {
            cdf_buf[f] = malloc((size_t)n * NUM_MATRICES * sizeof(double));
            cdf_idx[f] = 0;
        }

        for (unsigned m = 0; m < NUM_MATRICES; ++m) {
            uint32_t seed = RNG_SEED + (uint32_t)(si * 1000u) + m + 1u;
            csr_matrix A;
            if (csr_generate(&A, n, seed) != 0) {
                fprintf(stderr, "OOM\n"); return 1;
            }

            long double *x    = malloc(n * sizeof(long double));
            ref_t       *ref  = malloc(n * sizeof(ref_t));
            float *y_buf      = malloc(n * sizeof(float));
            uint32_t *y_tk    = malloc(n * sizeof(uint32_t));
            float *y_conv     = malloc(n * sizeof(float));
            if (!x || !ref || !y_buf || !y_tk || !y_conv) { fprintf(stderr, "OOM\n"); return 1; }

            uint32_t xs = seed ^ 0xC0FFEEu;
            if (xs == 0) xs = 1;
            for (unsigned i = 0; i < n; ++i) x[i] = rng_ld(&xs);
            spmv_ref(&A, x, ref);

            /* --- helpers para não repetir código --- */
            #define RUN_AND_EMIT(label_idx, call_expr, vl_elems, is_tk_dot, is_tk_tw, time_to_ms) do { \
                double t0 = now_sec(); (call_expr); double t1 = now_sec(); \
                for (unsigned i = 0; i < n; ++i) y_conv[i] = (float)tk_to_ld(y_tk[i], 8); /* placeholder */ \
                (void)0; \
            } while (0)
            /* Usamos macros explícitas abaixo para clareza. */

            /* ---------------- F32 ---------------- */
            {
                double t0 = now_sec();
                spmv_f32(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[6], &cdf_idx[6]);
                uint64_t vl = vl_for_bits(32);
                uint64_t ch = rvv_chunks(A.nnz, vl);
                uint64_t ie = rvv_inst_est(A.nnz, vl, 0, 0);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "F32", m, seed, t1-t0, ch, ie,
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[6], e.mean);
                sv_push(&sv_normwise[6], e.normwise);
                sv_push(&sv_time[6], (t1-t0)*1e3);
            }

            /* ---------------- F16 narrow ---------------- */
            {
                double t0 = now_sec();
                spmv_f16_narrow(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[4], &cdf_idx[4]);
                uint64_t vl = vl_for_bits(16);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "F16_N", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[4], e.mean);
                sv_push(&sv_normwise[4], e.normwise);
                sv_push(&sv_time[4], (t1-t0)*1e3);
            }

            /* ---------------- F16 wide ---------------- */
            {
                double t0 = now_sec();
                spmv_f16_wide(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[5], &cdf_idx[5]);
                uint64_t vl = vl_for_bits(16);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "F16_W", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[5], e.mean);
                sv_push(&sv_normwise[5], e.normwise);
                sv_push(&sv_time[5], (t1-t0)*1e3);
            }

            /* ---------------- E4M3 narrow / wide ---------------- */
            {
                double t0 = now_sec();
                spmv_e4m3_narrow(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[0], &cdf_idx[0]);
                uint64_t vl = vl_for_bits(8);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "E4M3_N", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[0], e.mean);
                sv_push(&sv_normwise[0], e.normwise);
                sv_push(&sv_time[0], (t1-t0)*1e3);
            }
            {
                double t0 = now_sec();
                spmv_e4m3_wide(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[1], &cdf_idx[1]);
                uint64_t vl = vl_for_bits(8);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "E4M3_W", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[1], e.mean);
                sv_push(&sv_normwise[1], e.normwise);
                sv_push(&sv_time[1], (t1-t0)*1e3);
            }

            /* ---------------- E5M2 narrow / wide ---------------- */
            {
                double t0 = now_sec();
                spmv_e5m2_narrow(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[2], &cdf_idx[2]);
                uint64_t vl = vl_for_bits(8);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "E5M2_N", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[2], e.mean);
                sv_push(&sv_normwise[2], e.normwise);
                sv_push(&sv_time[2], (t1-t0)*1e3);
            }
            {
                double t0 = now_sec();
                spmv_e5m2_wide(&A, x, y_buf);
                double t1 = now_sec();
                err_stats e = compute_err(ref, y_buf, n, cdf_buf[3], &cdf_idx[3]);
                uint64_t vl = vl_for_bits(8);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, "E5M2_W", m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[3], e.mean);
                sv_push(&sv_normwise[3], e.normwise);
                sv_push(&sv_time[3], (t1-t0)*1e3);
            }

            /* ---------------- Takum T8/T16/T32 (vtdot) ---------------- */
            const unsigned t_sew[3] = {8, 16, 32};
            const int      t_idx[3] = {7, 8, 9};
            for (int k = 0; k < 3; ++k) {
                unsigned sew = t_sew[k];
                int fi = t_idx[k];
                double t0 = now_sec();
                spmv_takum_vtdot(&A, x, sew, y_tk);
                double t1 = now_sec();
                for (unsigned i = 0; i < n; ++i) y_conv[i] = (float)tk_to_ld(y_tk[i], sew);
                err_stats e = compute_err(ref, y_conv, n, cdf_buf[fi], &cdf_idx[fi]);
                uint64_t vl = vl_for_bits(sew);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, FORMAT_NAMES[fi], m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 1, 0),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[fi], e.mean);
                sv_push(&sv_normwise[fi], e.normwise);
                sv_push(&sv_time[fi], (t1-t0)*1e3);
            }

            /* ---------------- Takum T16 twidth (vtmul+vtadd) ---------------- */
            {
                unsigned sew = 16;
                int fi = 10;
                double t0 = now_sec();
                spmv_takum_twidth(&A, x, sew, y_tk);
                double t1 = now_sec();
                for (unsigned i = 0; i < n; ++i) y_conv[i] = (float)tk_to_ld(y_tk[i], sew);
                err_stats e = compute_err(ref, y_conv, n, cdf_buf[fi], &cdf_idx[fi]);
                uint64_t vl = vl_for_bits(sew);
                fprintf(csv, "%u,%zu,%s,%u,%u,%.6f,%" PRIu64 ",%" PRIu64 ","
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, FORMAT_NAMES[fi], m, seed, t1-t0,
                        rvv_chunks(A.nnz, vl), rvv_inst_est(A.nnz, vl, 0, 1),
                        e.mean, e.rms, e.max, e.p95, e.p99, e.normwise);
                sv_push(&sv_err_mean[fi], e.mean);
                sv_push(&sv_normwise[fi], e.normwise);
                sv_push(&sv_time[fi], (t1-t0)*1e3);
            }

            csr_free(&A);
            free(x); free(ref); free(y_buf); free(y_tk); free(y_conv);
        }

        /* ---- Emitir estatísticas agregadas por formato ---- */
        for (size_t f = 0; f < NUM_FORMATS; ++f) {
            double em, es, ec, nm, ns, nc, tm, ts, tc;
            sv_stats(&sv_err_mean[f], &em, &es, &ec);
            sv_stats(&sv_normwise[f], &nm, &ns, &nc);
            sv_stats(&sv_time[f],    &tm, &ts, &tc);
            fprintf(stt, "%u,%s,%zu,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    n, FORMAT_NAMES[f], sv_err_mean[f].n,
                    em, es, ec, nm, ns, nc, tm, ts, tc);
        }

        /* ---- Emitir CDF por (size, format) ---- */
        for (size_t f = 0; f < NUM_FORMATS; ++f) {
            if (cdf_idx[f] == 0) { free(cdf_buf[f]); continue; }
            qsort(cdf_buf[f], cdf_idx[f], sizeof(double), cmp_double);
            size_t step = cdf_idx[f] / 1000 + 1;
            for (size_t i = 0; i < cdf_idx[f]; i += step) {
                fprintf(cdf, "%u,%s,%.6e,%.6f\n",
                        n, FORMAT_NAMES[f], cdf_buf[f][i],
                        (double)(i + 1) / (double)cdf_idx[f]);
            }
            free(cdf_buf[f]);
        }

        for (size_t f = 0; f < NUM_FORMATS; ++f) {
            sv_free(&sv_err_mean[f]);
            sv_free(&sv_normwise[f]);
            sv_free(&sv_time[f]);
        }
    }

    fclose(csv); fclose(stt); fclose(cdf); fclose(sum);
    printf("[benchmark] concluido.\n");
    printf("  - benchmark_results.csv\n");
    printf("  - benchmark_statistics.csv\n");
    printf("  - benchmark_cdf.csv\n");
    printf("  - benchmark_summary.txt\n");
    return 0;
}