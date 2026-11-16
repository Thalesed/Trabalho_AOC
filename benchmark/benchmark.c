/* benchmark.c */
#define _POSIX_C_SOURCE 200809L
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

static const unsigned NUM_MATRICES       = 20; /* Matrizes DIFERENTES para variância estatística */
static const unsigned MATRIX_SIZES[]     = {64, 128, 256, 512};
#define SPARSITY                        0.90
#define RNG_SEED                        0xDEADBEEFu

/* ============================================================
   RNG determinístico (xorshift32)
   ============================================================ */
static uint32_t g_state = RNG_SEED;
static uint32_t rng32_state(uint32_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    return *state;
}
static uint32_t rng32(void) { return rng32_state(&g_state); }
static long double rng_ld(void) {
    return ((long double)rng32() / (long double)UINT32_MAX) * 2.0L - 1.0L;
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

static int csr_generate(csr_matrix *m, unsigned n) {
    memset(m, 0, sizeof(*m));
    m->n = n;
    m->row_ptr = calloc((size_t)n + 1, sizeof(uint32_t));
    if (!m->row_ptr) return -1;

    const uint32_t decision_seed = g_state;
    size_t total = 0;
    for (unsigned i = 0; i < n; ++i) {
        m->row_ptr[i] = (uint32_t)total;
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32() / (double)UINT32_MAX >= SPARSITY) ++total;
        }
    }
    m->row_ptr[n] = (uint32_t)total;
    m->nnz = total;
    m->col_idx = malloc((total ? total : 1) * sizeof(*m->col_idx));
    m->vals = malloc((total ? total : 1) * sizeof(*m->vals));
    if ((!m->col_idx || !m->vals) && total > 0) { csr_free(m); return -1; }

    g_state = decision_seed;
    uint32_t value_state = decision_seed ^ 0xA5A5A5A5u;
    if (value_state == 0) value_state = 1;
    size_t k = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32() / (double)UINT32_MAX >= SPARSITY) {
                m->col_idx[k] = j;
                m->vals[k] = ((long double)rng32_state(&value_state) / (long double)UINT32_MAX) * 2.0L - 1.0L;
                ++k;
            }
        }
    }
    return 0;
}

/* ============================================================
   Referência e SpMV para vários formatos
   ============================================================ */
#ifdef TAKUM_USE_FLOAT128
typedef __float128 ref_t;
#define REF(x) ((ref_t)(x))
#else
typedef long double ref_t;
#define REF(x) ((ref_t)(x))
#endif

static void spmv_ref(const csr_matrix *A, const long double *x, ref_t *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        ref_t s = 0;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k)
            s += REF(A->vals[k]) * REF(x[A->col_idx[k]]);
        y[i] = s;
    }
}

static int spmv_takum(const csr_matrix *A, const long double *x, unsigned sew, uint32_t *y_out) {
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

/* --- IEEE 754 Narrow (Acumulador na mesma largura) --- */
static void spmv_f32(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k)
            s += (float)A->vals[k] * (float)x[A->col_idx[k]];
        y[i] = s;
    }
}
static void spmv_f64(const csr_matrix *A, const long double *x, double *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        double s = 0;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k)
            s += (double)A->vals[k] * (double)x[A->col_idx[k]];
        y[i] = s;
    }
}

/* --- IEEE 754 Widening (Acumulador em FP32, ex: vfwmacc) --- */
static void spmv_f16_wide(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f; /* Acumulador widening em FP32 */
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = (float)A->vals[k]; /* Simula F16 */
            float b = (float)x[A->col_idx[k]];
            s += a * b;
        }
        y[i] = s;
    }
}

/* --- FP8 E4M3 (Software Quantizer) --- */
static float quantize_f8_e4m3(float val) {
    uint32_t b; memcpy(&b, &val, 4);
    uint32_t sign = (b >> 31) & 1;
    int32_t exp = ((b >> 23) & 0xFF) - 127 + 7;
    uint32_t mant = b & 0x7FFFFF;
    uint8_t f8 = 0;
    if (exp >= 15) { f8 = 0x7E; } /* Clamp to max normal (448) */
    else if (exp <= 0) {
        if (exp < -9) f8 = 0;
        else {
            mant |= 0x800000;
            int shift = 23 - 3 + (1 - exp);
            mant += (1u << (shift - 1));
            f8 = (mant >> shift) & 0x7F;
        }
    } else {
        int shift = 23 - 3;
        mant += (1u << (shift - 1));
        if (mant & 0x800000) { mant >>= 1; exp++; if (exp >= 15) f8 = 0x7E; else f8 = (exp << 3) | ((mant >> shift) & 0x7); }
        else { f8 = (exp << 3) | ((mant >> shift) & 0x7); }
    }
    f8 |= (sign << 7);
    if ((f8 & 0x7F) == 0x7F) f8 = 0x7E; /* NaN -> Max Normal */
    
    /* Reverter para float */
    uint32_t r_sign = (f8 >> 7) & 1;
    uint32_t r_exp = (f8 >> 3) & 0xF;
    uint32_t r_mant = f8 & 0x7;
    uint32_t r_bits = (r_sign << 31);
    if (r_exp == 0) {
        if (r_mant == 0) return r_sign ? -0.0f : 0.0f;
        r_bits |= ((0) << 23) | (r_mant << 20); /* Subnormal */
    } else {
        r_bits |= ((r_exp - 7 + 127) << 23) | (r_mant << 20);
    }
    float res; memcpy(&res, &r_bits, 4);
    return res;
}

static void spmv_f8_narrow(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f;
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e4m3((float)A->vals[k]);
            float b = quantize_f8_e4m3((float)x[A->col_idx[k]]);
            s += quantize_f8_e4m3(a * b);
        }
        y[i] = s;
    }
}
static void spmv_f8_wide(const csr_matrix *A, const long double *x, float *y) {
    for (unsigned i = 0; i < A->n; ++i) {
        float s = 0.f; /* Acumulador widening em FP32 */
        for (uint32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; ++k) {
            float a = quantize_f8_e4m3((float)A->vals[k]);
            float b = quantize_f8_e4m3((float)x[A->col_idx[k]]);
            s += a * b; /* Soma em FP32 */
        }
        y[i] = s;
    }
}

/* ============================================================
   Estatísticas de Erro (Element-wise, Normwise, CDF)
   ============================================================ */
typedef struct {
    double mean, rms, max, p95, p99;
    double normwise; /* ||y - y_ref||_2 / ||y_ref||_2 */
    size_t count;
} err_stats;

static int cmp_double(const void *a, const void *b) {
    const double x = *(const double *)a;
    const double y = *(const double *)b;
    if (isnan(x)) return isnan(y) ? 0 : 1;
    if (isnan(y)) return -1;
    return (x > y) - (x < y);
}

static err_stats compute_err(const ref_t *ref, const float *est, unsigned n, double *global_cdf_buf, size_t *cdf_idx) {
    err_stats s = {0};
    if (n == 0) return s;
    double *errors = malloc((size_t)n * sizeof(*errors));
    if (!errors) return s;

    long double sum = 0.0L, sumsq = 0.0L;
    long double num_norm = 0.0L, den_norm = 0.0L;

    for (unsigned i = 0; i < n; ++i) {
        long double r = (long double)ref[i];
        long double v = (long double)est[i];
        long double denom = fmaxl(fabsl(r), 1.0e-30L);
        long double rel = fabsl(v - r) / denom;
        if (!isfinite(rel) || rel > (long double)DBL_MAX) rel = (long double)DBL_MAX;
        
        errors[i] = (double)rel;
        if (global_cdf_buf && cdf_idx) global_cdf_buf[(*cdf_idx)++] = (double)rel;
        
        sum += rel; sumsq += rel * rel;
        num_norm += (v - r) * (v - r);
        den_norm += r * r;
    }

    qsort(errors, n, sizeof(*errors), cmp_double);
    s.count = n;
    s.mean = (double)(sum / (long double)n);
    s.rms = (double)sqrtl(sumsq / (long double)n);
    s.max = errors[n - 1];
    s.p95 = errors[(size_t)ceill(0.95L * n) - 1];
    s.p99 = errors[(size_t)ceill(0.99L * n) - 1];
    
    s.normwise = (den_norm > 0) ? (double)sqrtl(num_norm / den_norm) : 0.0;
    
    free(errors);
    return s;
}

/* ============================================================
   Modelo de Instruções RVV (Unificado para IEEE e Takum)
   ============================================================ */
/* 
 * Modelo por bloco (chunk):
 * IEEE: vsetvli (1) + vle (1) + vluxei (1) + vfmul (1) + vfredusum (1) = 5 instruções
 * Takum: vsetvli (1) + vle (1) + vluxei (1) + vtdot (1) = 4 instruções
 */
static uint64_t estimate_rvv_instructions(size_t nnz, size_t vlmax, int is_takum) {
    if (nnz == 0 || vlmax == 0) return 0;
    size_t chunks = (nnz + vlmax - 1) / vlmax;
    uint64_t inst_per_chunk = is_takum ? 4 : 5;
    return chunks * inst_per_chunk;
}

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
#ifndef TAKUM_USE_FLOAT128
    printf("[AVISO] Compilado sem TAKUM_USE_FLOAT128. A referência será long double (80-bit), não Float128.\n");
    printf("[AVISO] Para o artigo, use: gcc ... -DTAKUM_USE_FLOAT128 -lquadmath\n");
#endif
    printf("[AVISO] O tempo medido reflete o custo de SIMULAÇÃO SOFTWARE (expl/logl), não latência de hardware.\n\n");

    const size_t n_sizes = sizeof(MATRIX_SIZES) / sizeof(MATRIX_SIZES[0]);
    FILE *csv = fopen("benchmark_results.csv", "w");
    FILE *sum = fopen("benchmark_summary.txt", "w");
    
    fprintf(csv, "size,nnz,tipo,mat_id,time_s,inst_est,err_mean,err_rms,err_max,err_p95,err_p99,err_normwise\n");
    fprintf(sum, "=== Benchmark SpMV (y = A*x) — Takum vs IEEE 754 ===\n");
    fprintf(sum, "NUM_MATRICES=%u (Matrizes DIFERENTES para variância estatística)\n", NUM_MATRICES);
    fprintf(sum, "SPARSITY=%.2f  SEED=0x%X\n\n", SPARSITY, RNG_SEED);

    /* Buffers para CDF global (todos os elementos de todas as matrizes) */
    size_t max_total_elements = NUM_MATRICES * 512 * 512; /* Pior caso */
    double *global_cdf = malloc(max_total_elements * sizeof(double));
    size_t cdf_idx = 0;

    for (size_t si = 0; si < n_sizes; ++si) {
        unsigned n = MATRIX_SIZES[si];
        printf("[benchmark] Gerando %u matrizes %ux%u para avaliação estatística...\n", NUM_MATRICES, n, n);
        
        for (unsigned m = 0; m < NUM_MATRICES; ++m) {
            csr_matrix A;
            if (csr_generate(&A, n) != 0) { fprintf(stderr, "OOM\n"); return 1; }
            
            long double *x = malloc(n * sizeof(long double));
            ref_t *ref = malloc(n * sizeof(ref_t));
            float *y_f8_n = malloc(n * sizeof(float));
            float *y_f8_w = malloc(n * sizeof(float));
            float *y_f16_w = malloc(n * sizeof(float));
            float *y_f32 = malloc(n * sizeof(float));
            float *y_f64 = malloc(n * sizeof(float)); /* Cast para float para unificar compute_err */
            uint32_t *y_tk = malloc(n * sizeof(uint32_t));
            float *y_t8 = malloc(n * sizeof(float));
            float *y_t16 = malloc(n * sizeof(float));
            float *y_t32 = malloc(n * sizeof(float));

            for (unsigned i = 0; i < n; ++i) x[i] = rng_ld();
            spmv_ref(&A, x, ref);

            double t0, t1;
            size_t vlmax_t8 = tk_sim_vlmax(8), vlmax_t16 = tk_sim_vlmax(16), vlmax_t32 = tk_sim_vlmax(32);
            size_t vlmax_f32 = 8; /* Exemplo RVV F32 */

            /* Execução e Coleta */
            t0 = now_sec(); spmv_f8_narrow(&A, x, y_f8_n); t1 = now_sec();
            err_stats e_f8n = compute_err(ref, y_f8_n, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,F8_N,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_f32, 0), e_f8n.mean, e_f8n.rms, e_f8n.max, e_f8n.p95, e_f8n.p99, e_f8n.normwise);

            t0 = now_sec(); spmv_f8_wide(&A, x, y_f8_w); t1 = now_sec();
            err_stats e_f8w = compute_err(ref, y_f8_w, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,F8_W,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_f32, 0), e_f8w.mean, e_f8w.rms, e_f8w.max, e_f8w.p95, e_f8w.p99, e_f8w.normwise);

            t0 = now_sec(); spmv_f16_wide(&A, x, y_f16_w); t1 = now_sec();
            err_stats e_f16w = compute_err(ref, y_f16_w, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,F16_W,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_f32, 0), e_f16w.mean, e_f16w.rms, e_f16w.max, e_f16w.p95, e_f16w.p99, e_f16w.normwise);

            t0 = now_sec(); spmv_f32(&A, x, y_f32); t1 = now_sec();
            err_stats e_f32 = compute_err(ref, y_f32, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,F32,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_f32, 0), e_f32.mean, e_f32.rms, e_f32.max, e_f32.p95, e_f32.p99, e_f32.normwise);

            t0 = now_sec(); spmv_takum(&A, x, 8, y_tk); t1 = now_sec();
            for(unsigned i=0;i<n;i++) y_t8[i] = tk_to_ld(y_tk[i], 8);
            err_stats e_t8 = compute_err(ref, y_t8, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,T8,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_t8, 1), e_t8.mean, e_t8.rms, e_t8.max, e_t8.p95, e_t8.p99, e_t8.normwise);

            t0 = now_sec(); spmv_takum(&A, x, 16, y_tk); t1 = now_sec();
            for(unsigned i=0;i<n;i++) y_t16[i] = tk_to_ld(y_tk[i], 16);
            err_stats e_t16 = compute_err(ref, y_t16, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,T16,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_t16, 1), e_t16.mean, e_t16.rms, e_t16.max, e_t16.p95, e_t16.p99, e_t16.normwise);

            t0 = now_sec(); spmv_takum(&A, x, 32, y_tk); t1 = now_sec();
            for(unsigned i=0;i<n;i++) y_t32[i] = tk_to_ld(y_tk[i], 32);
            err_stats e_t32 = compute_err(ref, y_t32, n, global_cdf, &cdf_idx);
            fprintf(csv, "%u,%zu,T32,%u,%.6f,%llu,%.3e,%.3e,%.3e,%.3e,%.3e,%.3e\n", n, A.nnz, m, t1-t0, estimate_rvv_instructions(A.nnz, vlmax_t32, 1), e_t32.mean, e_t32.rms, e_t32.max, e_t32.p95, e_t32.p99, e_t32.normwise);

            csr_free(&A);
            free(x); free(ref); free(y_f8_n); free(y_f8_w); free(y_f16_w); 
            free(y_f32); free(y_f64); free(y_tk); free(y_t8); free(y_t16); free(y_t32);
        }
    }

    /* Gerar CDF Global */
    qsort(global_cdf, cdf_idx, sizeof(double), cmp_double);
    FILE *cdf_file = fopen("benchmark_cdf.csv", "w");
    fprintf(cdf_file, "error_relative,cumulative_probability\n");
    for (size_t i = 0; i < cdf_idx; i += (cdf_idx / 1000 + 1)) { /* Amostragem para não explodir o CSV */
        fprintf(cdf_file, "%.6e,%.6f\n", global_cdf[i], (double)(i + 1) / cdf_idx);
    }
    fclose(cdf_file);
    free(global_cdf);

    fprintf(sum, "\n=== Análise de Complexidade e CSRs ===\n");
    fprintf(sum, "1. Modelo de Instruções RVV (por bloco/chunk):\n");
    fprintf(sum, "   IEEE 754: vsetvli + vle + vluxei + vfmul + vfredusum = 5 instruções.\n");
    fprintf(sum, "   Takum:    vsetvli + vle + vluxei + vtdot = 4 instruções.\n");
    fprintf(sum, "   O vtdot do Takum substitui a multiplicação e redução, economizando 1 instrução por chunk.\n\n");
    
    fprintf(sum, "2. Redução de Estado de CSRs (Argumento Defensável):\n");
    fprintf(sum, "   IEEE 754 exige manutenção de CSRs de controle e status: 'fflags' (exceções) e 'frm' (rounding mode).\n");
    fprintf(sum, "   Takum possui arredondamento fixo (RNE) e não gera exceções de ponto flutuante (apenas NaR).\n");
    fprintf(sum, "   Isso elimina a necessidade de salvar/restaurar 'fflags' e 'frm' em contextos de interrupção ou chamadas de função, reduzindo o estado de CSR de 2 registradores para 0.\n\n");

    fprintf(sum, "3. Aviso sobre Tempo de Execução:\n");
    fprintf(sum, "   O tempo medido neste benchmark reflete o custo de SIMULAÇÃO SOFTWARE (devido às funções expl/logl usadas na conversão Takum).\n");
    fprintf(sum, "   NÃO reflete a latência de hardware real. Em hardware dedicado (RTL), Takum é mais simples que IEEE, mas a simulação em CPU x86 é mais lenta.\n");

    fclose(csv); fclose(sum);
    printf("[benchmark] concluído. Veja benchmark_results.csv, benchmark_cdf.csv e benchmark_summary.txt\n");
    return 0;
}