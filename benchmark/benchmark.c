/* benchmark.c — Avaliação experimental SpMV (y = A * x) em Takum vs IEEE 754
 *   - Multiplicação matriz-vetor com matrizes esparsas (CSR)
 *   - Distribuição do erro relativo acumulado vs Float128
 *   - Redução de instruções e CSRs (contadores por tipo)
 *   - Tempo de execução cronometrado
 *
 * Compilar:  make benchmark
 * Executar:                        ./build/benchmark
 * Saídas:                          benchmark_results.csv
 *                                  benchmark_summary.txt
 *                                  benchmark_hist_T8.csv / T16 / T32
 */

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

static const unsigned REPEATS_PER_SIZE = 20;                  /* N */
static const unsigned MATRIX_SIZES[]   = {64, 128, 256, 512}; /* Tamanho das matrizes e vatores */
#define SPARSITY                      0.90                    /* 90% zeros */
#define RNG_SEED                      0xDEADBEEFu

/* ============================================================
 * RNG determinístico (xorshift32)
 * ============================================================ */
static uint32_t g_state = RNG_SEED;
static uint32_t rng32_state(uint32_t *state) {
    *state ^= *state << 13;
    *state ^= *state >> 17;
    *state ^= *state << 5;
    return *state;
}
static uint32_t rng32(void) {
    return rng32_state(&g_state);
}
static long double rng_ld(void) { /* U(-1,1) */
    return ((long double)rng32() / (long double)UINT32_MAX) * 2.0L - 1.0L;
}

/* ============================================================
 * Matriz esparsa em CSR
 * ============================================================ */
typedef struct {
    unsigned n;          /* dimensões n x n */
    size_t nnz;          /* número de não-zeros */
    uint32_t *row_ptr;   /* [n+1] */
    uint32_t *col_idx;   /* [nnz] */
    long double *vals;   /* [nnz] */
} csr_matrix;

static void csr_free(csr_matrix *m) {
    free(m->row_ptr); free(m->col_idx); free(m->vals);
    memset(m, 0, sizeof(*m));
}

/* Gera CSR aleatório com esparsidade SPARSITY e valores em U(-1,1). */
static int csr_generate(csr_matrix *m, unsigned n) {
    memset(m, 0, sizeof(*m));
    m->n = n;
    m->row_ptr = calloc((size_t)n + 1, sizeof(uint32_t));
    if (!m->row_ptr) return -1;

    /* Preserve the decision RNG state so both passes use identical draws. */
    const uint32_t decision_seed = g_state;
    size_t total = 0;
    for (unsigned i = 0; i < n; ++i) {
        if (total > UINT32_MAX) {
            csr_free(m);
            return -1;
        }
        m->row_ptr[i] = (uint32_t)total;
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32() / (double)UINT32_MAX >= SPARSITY)
                ++total;
        }
    }
    if (total > UINT32_MAX) {
        csr_free(m);
        return -1;
    }
    m->row_ptr[n] = (uint32_t)total;
    m->nnz = total;

    m->col_idx = malloc((total ? total : 1) * sizeof(*m->col_idx));
    m->vals = malloc((total ? total : 1) * sizeof(*m->vals));
    if (!m->col_idx || !m->vals) {
        csr_free(m);
        return -1;
    }

    /* Replay only the sparsity decisions; values use an independent RNG. */
    g_state = decision_seed;
    uint32_t value_state = decision_seed ^ 0xA5A5A5A5u;
    if (value_state == 0) value_state = 1;
    size_t k = 0;
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
            if ((double)rng32() / (double)UINT32_MAX >= SPARSITY) {
                if (k >= total) {
                    csr_free(m);
                    return -1;
                }
                m->col_idx[k] = j;
                m->vals[k] = ((long double)rng32_state(&value_state) /
                              (long double)UINT32_MAX) * 2.0L - 1.0L;
                ++k;
            }
        }
    }
    if (k != total) {
        csr_free(m);
        return -1;
    }
    return 0;
}

/*  Tipo: referência Float128 (precisão quádrupla) */
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

static int spmv_takum(const csr_matrix *A, const long double *x,
                      unsigned sew, uint32_t *y_out,
                      uint64_t *retired, uint64_t *elements,
                      uint64_t *configs) {
    tk_sim sim; tk_sim_init(&sim);
    size_t vlmax = tk_sim_vlmax(sew);
    if (vlmax == 0) return -1;

    const uint64_t r0 = sim.retired, e0 = sim.elements, c0 = sim.configurations;

    for (unsigned i = 0; i < A->n; ++i) {
        uint32_t start = A->row_ptr[i], end = A->row_ptr[i + 1];
        uint32_t nnz = end - start;
        if (nnz == 0) {
            if (tk_sim_setvl(&sim, 1, sew) != TK_SIM_OK) return -1;
            long double zero = 0;
            if (tk_sim_load(&sim, 1, &zero, 1) != TK_SIM_OK) return -1;
            uint32_t insn = tk_instruction(TK_VTDOT, 3, 1, 1);
            if (tk_sim_execute(&sim, insn) != TK_SIM_OK) return -1;
            if (tk_sim_read(&sim, 3, sew, 0, &y_out[i]) != TK_SIM_OK) return -1;
            continue;
        }
        /* processa linha em blocos de vlmax */
        long double acc = 0;
        uint32_t remain = nnz, off = start;
        while (remain) {
            uint32_t chunk = remain < vlmax ? remain : (uint32_t)vlmax;
            long double a_buf[TK_VLEN_BYTES], b_buf[TK_VLEN_BYTES];
            for (uint32_t t = 0; t < chunk; ++t) {
                a_buf[t] = A->vals[off + t];
                b_buf[t] = x[A->col_idx[off + t]];
            }
            if (tk_sim_setvl(&sim, chunk, sew) != TK_SIM_OK) return -1;
            if (tk_sim_load(&sim, 1, a_buf, chunk) != TK_SIM_OK) return -1;
            if (tk_sim_load(&sim, 2, b_buf, chunk) != TK_SIM_OK) return -1;
            uint32_t insn = tk_instruction(TK_VTDOT, 3, 1, 2);
            if (tk_sim_execute(&sim, insn) != TK_SIM_OK) return -1;
            acc += sim.dot_accumulator;
            remain -= chunk; off += chunk;
        }
        y_out[i] = tk_from_ld(acc, sew);
    }

    if (retired)  *retired  = sim.retired  - r0;
    if (elements) *elements = sim.elements - e0;
    if (configs)  *configs  = sim.configurations - c0;
    return 0;
}

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

/* Estatísticas de erro relativo  */
typedef struct {
    double mean, rms, max, p95, p99;
    size_t count;
} err_stats;

static int cmp_double(const void *a, const void *b) {
    const double x = *(const double *)a;
    const double y = *(const double *)b;

    /* A função de comparação também define uma ordem para valores não finitos. */
    if (isnan(x)) return isnan(y) ? 0 : 1;
    if (isnan(y)) return -1;
    return (x > y) - (x < y);
}

static err_stats compute_err(const ref_t *ref, const long double *est, unsigned n) {
    err_stats s = {0};
    if (n == 0) return s;

    double *errors = malloc((size_t)n * sizeof(*errors));
    if (!errors) {
        fprintf(stderr, "erro: falha ao alocar vetor de erros\n");
        return s;
    }

    /*
     * Calcula e armazena exatamente os mesmos erros usados nas estatísticas.
     * O acumulador long double reduz o risco de overflow na soma.
     * NaN/Inf são tratados como erro máximo, nunca como erro zero.
     */
    long double sum = 0.0L;
    long double sumsq = 0.0L;
    for (unsigned i = 0; i < n; ++i) {
        const long double r = (long double)ref[i];
        const long double v = est[i];
        const long double denom = fmaxl(fabsl(r), 1.0e-30L);
        long double rel = fabsl(v - r) / denom;

        if (!isfinite(rel) || rel > (long double)DBL_MAX)
            rel = (long double)DBL_MAX;
        if (rel < 0.0L) rel = (long double)DBL_MAX;

        errors[i] = (double)rel;
        sum += rel;
        sumsq += rel * rel;
    }

    qsort(errors, n, sizeof(*errors), cmp_double);

    s.count = n;
    s.mean = (double)(sum / (long double)n);
    s.rms = (double)sqrtl(sumsq / (long double)n);
    s.max = errors[n - 1];

    /* Quantis por nearest-rank: índices crescentes e limitados a [0, n-1]. */
    size_t i95 = (size_t)ceill(0.95L * (long double)n) - 1;
    size_t i99 = (size_t)ceill(0.99L * (long double)n) - 1;
    if (i95 >= n) i95 = n - 1;
    if (i99 >= n) i99 = n - 1;
    s.p95 = errors[i95];
    s.p99 = errors[i99];

    /* Essas relações devem sempre valer para uma distribuição ordenada. */
    if (s.mean > s.max || s.p95 > s.p99 || s.p99 > s.max) {
        fprintf(stderr,
                "erro interno nas estatísticas: mean=%g max=%g p95=%g p99=%g\n",
                s.mean, s.max, s.p95, s.p99);
    }

    free(errors);
    return s;
}

/* Histograma log2 dos erros (para distribuição) */
static void hist_log2(const ref_t *ref, const long double *est, unsigned n,
                      unsigned *bins, unsigned nbins, double lo, double hi) {
    memset(bins, 0, nbins * sizeof(unsigned));
    double step = (hi - lo) / nbins;
    for (unsigned i = 0; i < n; ++i) {
        double r = (double)ref[i], v = (double)est[i];
        double denom = fabs(r) > 1e-300 ? fabs(r) : 1e-300;
        double rel = fabs(v - r) / denom;
        if (rel <= 0) rel = 1e-300;
        double x = log2(rel);
        int b = (int)((x - lo) / step);
        if (b < 0) b = 0;
        if (b >= (int)nbins) b = (int)nbins - 1;
        bins[b]++;
    }
}

/* Cronômetro */
static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}


int main(void) {
    const size_t n_sizes = sizeof(MATRIX_SIZES) / sizeof(MATRIX_SIZES[0]);
    const unsigned nbins = 40;
    const double hist_lo = -30.0, hist_hi = 10.0; /* log2(erro) */

    FILE *csv = fopen("benchmark_results.csv", "w");
    FILE *sum = fopen("benchmark_summary.txt", "w");
    if (!csv || !sum) { fprintf(stderr, "erro abrindo arquivos\n"); return 1; }

    fprintf(csv, "size,nnz,tipo,repeat,time_s,csr_writes_est,"
                 "inst_retired,elements,configs,err_mean,err_rms,err_max,"
                 "err_p95,err_p99,out0\n");
    fprintf(sum, "=== Benchmark SpMV (y = A*x) — Takum vs IEEE 754 ===\n");
    fprintf(sum, "REPEATS_PER_SIZE=%u  SPARSITY=%.2f  SEED=0x%X\n\n",
            REPEATS_PER_SIZE, SPARSITY, RNG_SEED);
    fprintf(sum, "%-6s %-8s %-8s %-12s %-12s %-12s %-12s %-12s\n",
            "size", "nnz", "tipo", "time_ms(mean)", "err_mean",
            "err_max", "err_p95", "err_p99");

    /* acumuladores por tipo para resumo final */
    struct {
        const char *name;
        double t_sum; size_t t_cnt;
        double err_mean_acc, err_max_acc, err_p95_acc, err_p99_acc;
        uint64_t inst_acc, elem_acc, cfg_acc, csr_acc;
        size_t runs;
    } acc[] = {
        {"T8",  0,0, 0,0,0,0, 0,0,0,0, 0},
        {"T16", 0,0, 0,0,0,0, 0,0,0,0, 0},
        {"T32", 0,0, 0,0,0,0, 0,0,0,0, 0},
        {"F32", 0,0, 0,0,0,0, 0,0,0,0, 0},
        {"F64", 0,0, 0,0,0,0, 0,0,0,0, 0},
    };
    const size_t NB = sizeof(acc) / sizeof(acc[0]);

    /* buffers temporários (malloc para tamanho grande) */
    unsigned *hist_buf = calloc(nbins, sizeof(unsigned));
    if (!hist_buf) return 1;

    for (size_t si = 0; si < n_sizes; ++si) {
        unsigned n = MATRIX_SIZES[si];
        printf("[benchmark] gerando matriz %ux%u (sparsity %.2f)\n",
               n, n, SPARSITY);
        csr_matrix A;
        if (csr_generate(&A, n) != 0) { fprintf(stderr, "OOM\n"); return 1; }

        long double *x = malloc(n * sizeof(long double));
        long double *y_ref = malloc(n * sizeof(long double));
        long double *y_t8  = malloc(n * sizeof(long double));
        long double *y_t16 = malloc(n * sizeof(long double));
        long double *y_t32 = malloc(n * sizeof(long double));
        long double *y_f32 = malloc(n * sizeof(long double));
        long double *y_f64 = malloc(n * sizeof(long double));
        ref_t *ref = malloc(n * sizeof(ref_t));
        uint32_t *y_tk = malloc(n * sizeof(uint32_t));
        if (!x || !y_ref || !y_t8 || !y_t16 || !y_t32 || !y_f32 || !y_f64 || !ref || !y_tk) {
            fprintf(stderr, "OOM\n"); return 1;
        }

        /* gera vetor x fora do cronômetro */
        for (unsigned i = 0; i < n; ++i) x[i] = rng_ld();

        /* referência Float128 */
        spmv_ref(&A, x, ref);
        for (unsigned i = 0; i < n; ++i) y_ref[i] = (long double)ref[i];

        /* -------- Float32 -------- */
        {
            double t_sum = 0; err_stats es = {0};
            float *yf = malloc(n * sizeof(float));
            for (unsigned r = 0; r < REPEATS_PER_SIZE; ++r) {
                double t0 = now_sec();
                spmv_f32(&A, x, yf);
                double t1 = now_sec();
                t_sum += (t1 - t0);
                for (unsigned i = 0; i < n; ++i) y_f32[i] = (long double)yf[i];
                err_stats e = compute_err(ref, y_f32, n);
                es = e; /* último run é representativo */
                fprintf(csv, "%u,%zu,F32,%u,%.9f,%u,0,0,0,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, r, (t1 - t0), (unsigned)A.nnz,
                        e.mean, e.rms, e.max, e.p95, e.p99,
                        (double)(n ? y_f32[0] : 0));
            }
            free(yf);
            double t_mean = t_sum / REPEATS_PER_SIZE;
            acc[3].t_sum += t_mean; acc[3].t_cnt++;
            acc[3].err_mean_acc += es.mean; acc[3].err_max_acc += es.max;
            acc[3].err_p95_acc += es.p95;  acc[3].err_p99_acc += es.p99;
            acc[3].csr_acc += (uint64_t)A.nnz; /* IEEE trad.: 1 CSR write por elem (upper bound) */
            acc[3].runs++;
            fprintf(sum, "%-6u %-8zu %-8s %-12.4f %-12.3e %-12.3e %-12.3e %-12.3e\n",
                    n, A.nnz, "F32", t_mean * 1e3, es.mean, es.max, es.p95, es.p99);
        }

        /* -------- Float64 -------- */
        {
            double t_sum = 0; err_stats es = {0};
            double *yd = malloc(n * sizeof(double));
            for (unsigned r = 0; r < REPEATS_PER_SIZE; ++r) {
                double t0 = now_sec();
                spmv_f64(&A, x, yd);
                double t1 = now_sec();
                t_sum += (t1 - t0);
                for (unsigned i = 0; i < n; ++i) y_f64[i] = (long double)yd[i];
                err_stats e = compute_err(ref, y_f64, n);
                es = e;
                fprintf(csv, "%u,%zu,F64,%u,%.9f,%u,0,0,0,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, r, (t1 - t0), (unsigned)A.nnz,
                        e.mean, e.rms, e.max, e.p95, e.p99,
                        (double)(n ? y_f64[0] : 0));
            }
            free(yd);
            double t_mean = t_sum / REPEATS_PER_SIZE;
            acc[4].t_sum += t_mean; acc[4].t_cnt++;
            acc[4].err_mean_acc += es.mean; acc[4].err_max_acc += es.max;
            acc[4].err_p95_acc += es.p95;  acc[4].err_p99_acc += es.p99;
            acc[4].csr_acc += (uint64_t)A.nnz;
            acc[4].runs++;
            fprintf(sum, "%-6u %-8zu %-8s %-12.4f %-12.3e %-12.3e %-12.3e %-12.3e\n",
                    n, A.nnz, "F64", t_mean * 1e3, es.mean, es.max, es.p95, es.p99);
        }

        /* -------- Takum T8/T16/T32 -------- */
        for (unsigned k = 0; k < 3; ++k) {
            unsigned sew = 8u << k;
            const char *label = acc[k].name;
            double t_sum = 0;
            err_stats es = {0};
            uint64_t r_sum = 0, e_sum = 0, c_sum = 0;
            for (unsigned r = 0; r < REPEATS_PER_SIZE; ++r) {
                uint64_t rr = 0, ee = 0, cc = 0;
                double t0 = now_sec();
                int rc = spmv_takum(&A, x, sew, y_tk, &rr, &ee, &cc);
                double t1 = now_sec();
                if (rc != 0) { fprintf(stderr, "falha spmv_takum\n"); return 1; }
                t_sum += (t1 - t0);
                r_sum += rr; e_sum += ee; c_sum += cc;
                for (unsigned i = 0; i < n; ++i)
                    (k == 0 ? (y_t8[i]  = tk_to_ld(y_tk[i], 8))  :
                     k == 1 ? (y_t16[i] = tk_to_ld(y_tk[i], 16)) :
                              (y_t32[i] = tk_to_ld(y_tk[i], 32)));
                long double *yv = k == 0 ? y_t8 : k == 1 ? y_t16 : y_t32;
                err_stats e = compute_err(ref, yv, n);
                es = e;
                fprintf(csv, "%u,%zu,%s,%u,%.9f,0,%llu,%llu,%llu,"
                             "%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                        n, A.nnz, label, r, (t1 - t0),
                        (unsigned long long)rr, (unsigned long long)ee,
                        (unsigned long long)cc,
                        e.mean, e.rms, e.max, e.p95, e.p99,
                        (double)(n ? yv[0] : 0));
            }
            double t_mean = t_sum / REPEATS_PER_SIZE;
            acc[k].t_sum += t_mean; acc[k].t_cnt++;
            acc[k].err_mean_acc += es.mean; acc[k].err_max_acc += es.max;
            acc[k].err_p95_acc += es.p95;  acc[k].err_p99_acc += es.p99;
            acc[k].inst_acc += r_sum; acc[k].elem_acc += e_sum;
            acc[k].cfg_acc  += c_sum; acc[k].csr_acc  += 0; /* sem CSR */
            acc[k].runs++;
            fprintf(sum, "%-6u %-8zu %-8s %-12.4f %-12.3e %-12.3e %-12.3e %-12.3e\n",
                    n, A.nnz, label, t_mean * 1e3, es.mean, es.max, es.p95, es.p99);
        }

        /* histogramas de erro (última repetição de cada tipo Takum) */
        for (unsigned k = 0; k < 3; ++k) {
            long double *yv = k == 0 ? y_t8 : k == 1 ? y_t16 : y_t32;
            hist_log2(ref, yv, n, hist_buf, nbins, hist_lo, hist_hi);
            char fname[64];
            snprintf(fname, sizeof(fname), "benchmark_hist_%s.csv", acc[k].name);
            FILE *hf = fopen(fname, "w");
            if (hf) {
                fprintf(hf, "log2_err_bucket,count\n");
                double step = (hist_hi - hist_lo) / nbins;
                for (unsigned b = 0; b < nbins; ++b)
                    fprintf(hf, "%.2f,%u\n", hist_lo + (b + 0.5) * step, hist_buf[b]);
                fclose(hf);
            }
        }

        csr_free(&A);
        free(x); free(y_ref); free(y_t8); free(y_t16); free(y_t32);
        free(y_f32); free(y_f64); free(ref); free(y_tk);
    }

    /* ---------- Resumo agregado + CSRs ---------- */
    fprintf(sum, "\n=== Médias por tipo (todos os tamanhos) ===\n");
    fprintf(sum, "%-6s %-14s %-12s %-12s %-12s %-12s %-14s %-14s %-10s\n",
            "fmt", "time_ms(mean)", "err_mean", "err_max", "err_p95",
            "err_p99", "inst_retired", "elements", "CSR_writes");
    for (size_t b = 0; b < NB; ++b) {
        if (!acc[b].runs) continue;
        fprintf(sum, "%-6s %-14.4f %-12.3e %-12.3e %-12.3e %-12.3e "
                     "%-14llu %-14llu %-10llu\n",
                acc[b].name,
                acc[b].t_sum / acc[b].t_cnt * 1e3,
                acc[b].err_mean_acc / acc[b].runs,
                acc[b].err_max_acc  / acc[b].runs,
                acc[b].err_p95_acc  / acc[b].runs,
                acc[b].err_p99_acc  / acc[b].runs,
                (unsigned long long)acc[b].inst_acc,
                (unsigned long long)acc[b].elem_acc,
                (unsigned long long)acc[b].csr_acc);
    }

    fprintf(sum, "\n=== Complexidade de instruções / CSRs ===\n");
    fprintf(sum, "Takum (Zvtakum): nenhum CSR dedicado (sem frm, fflags, "
                 "NaN-boxing). fflags não é atualizado por operações Takum.\n");
    fprintf(sum, "IEEE 754 (RVV tradicional): estimativa de 1 escrita de fflags "
                 "por elemento multiplicado (upper bound; na prática só em "
                 "exceções, mas o CSR existe e deve ser preservado).\n");
    fprintf(sum, "Redução de CSRs: T8/T16/T32 = 0 escritas vs ~nnz para IEEE.\n");
    fprintf(sum, "\nArquivos gerados: benchmark_results.csv, "
                 "benchmark_hist_T8.csv, benchmark_hist_T16.csv, "
                 "benchmark_hist_T32.csv\n");

    free(hist_buf);
    fclose(csv); fclose(sum);
    printf("[benchmark] concluído. Veja benchmark_results.csv e benchmark_summary.txt\n");
    return 0;
}