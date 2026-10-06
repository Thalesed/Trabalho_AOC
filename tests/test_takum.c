#include "../src/rvv_takum.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long checks;
#define CHECK(expr) do { ++checks; if (!(expr)) { \
    fprintf(stderr,"Falha %s:%d: %s\n",__FILE__,__LINE__,#expr); exit(EXIT_FAILURE); \
} } while (0)
static uint32_t seed = UINT32_C(0x31415926);
static uint32_t random_word(void) { seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5; return seed; }
static uint32_t lane(const tk_sim *s, unsigned r, unsigned n, size_t i) {
    uint32_t result;
    CHECK(tk_sim_read(s,r,n,i,&result) == TK_SIM_OK);
    return result;
}
static void scalar(void) {
    CHECK(!tk_width_valid(0)); CHECK(!tk_width_valid(64));
    CHECK(tk_mask(32) == UINT32_MAX);
    CHECK(tk_nar(32) == UINT32_C(0x80000000));
    takum_fields f;
    CHECK(!tk_extract(0,7,&f)); CHECK(!tk_extract(0,8,NULL));
    CHECK(isnan(tk_to_ld(0,7)));
    CHECK(tk_extract(1,8,&f));
    CHECK(f.S == 0 && f.D == 0 && f.R == 0 && f.r == 7);
    CHECK(f.C == 1 && f.characteristic_bits == 3 && f.fraction_bits == 0);
    CHECK(f.c == -239 && f.l == -239);
    CHECK(tk_extract(0x3800,16,&f));
    CHECK(f.R == 7 && f.r == 0 && f.c == -1 && f.fraction_bits == 11);
    CHECK(tk_to_ld(0x4000,16) == 1 && tk_to_ld(0xc000,16) == -1);
    CHECK(tk_to_log(0x7fff,16) == 254.9375L);
    CHECK(tk_to_log(0x7f,8) == 239);
    CHECK(tk_to_log(UINT32_C(0x7fffffff),32) == 254.99999904632568359375L);
    for (unsigned n = 8; n <= 32; n *= 2) {
        uint32_t one = UINT32_C(1) << (n-2), nar = tk_nar(n), mask = tk_mask(n);
        CHECK(tk_from_ld(0,n) == 0); CHECK(tk_from_ld(-0.0L,n) == 0);
        CHECK(tk_to_ld(0,n) == 0); CHECK(isnan(tk_to_ld(nar,n)));
        CHECK(tk_from_ld(NAN,n) == nar); CHECK(tk_from_ld(INFINITY,n) == nar);
        CHECK(tk_from_ld(-INFINITY,n) == nar);
        CHECK(tk_from_ld(1,n) == one); CHECK(tk_from_ld(-1,n) == ((0-one) & mask));
        CHECK(tk_from_ld(1e300L,n) == nar-1); CHECK(tk_from_ld(-1e300L,n) == nar+1);
        CHECK(tk_from_ld(1e-300L,n) == 1); CHECK(tk_from_ld(-1e-300L,n) == mask);
        CHECK(tk_mul(0,one,n) == 0); CHECK(tk_mul(0,nar,n) == nar);
        CHECK(tk_add(one,nar,n) == nar); CHECK(tk_sub(one,one,n) == 0);
        CHECK(tk_convert(0,n,8) == 0); CHECK(tk_convert(nar,n,8) == 0x80);
        uint32_t total = n < 32 ? UINT32_C(1) << n : 100000;
        for (uint32_t i = 0; i < total; ++i) {
            uint32_t w = n < 32 ? i : random_word();
            long double x = tk_to_ld(w,n);
            CHECK(tk_from_ld(x,n) == w);
            CHECK(tk_mul(w,one,n) == w);
            CHECK(tk_add(w,0,n) == w);
            if (w != nar) {
                uint32_t neg = (0-w) & mask;
                CHECK(tk_to_ld(neg,n) == -x);
                CHECK(tk_add(w,neg,n) == 0);
            }
        }
        uint32_t a = tk_from_ld(2.5L,n), b = tk_from_ld(3,n);
        long double tolerance = n == 8 ? 0.1L : n == 16 ? 0.002L : 1e-7L;
        CHECK(fabsl(tk_to_ld(tk_add(a,b,n),n)-5.5L)/5.5L < tolerance);
        CHECK(fabsl(tk_to_ld(tk_mul(a,b,n),n)-7.5L)/7.5L < tolerance);
    }
    CHECK(tk_convert(0x4080,16,8) == 0x40);
    CHECK(tk_mul(0x01,0x59,8) == 0x02); /* Empate cujo prefixo e impar. */
    CHECK(tk_convert(0x4180,16,8) == 0x42);
    CHECK(tk_convert(1,16,8) == 1); CHECK(tk_convert(0xffff,16,8) == 0xff);
    CHECK(tk_convert(0x7fff,16,8) == 0x7f); CHECK(tk_convert(0x8001,16,8) == 0x81);
    for (uint32_t i = 0; i < 256; ++i) {
        CHECK(tk_convert(i,8,32) == i << 24);
        CHECK(tk_convert(tk_convert(i,8,32),32,8) == i);
    }
    CHECK(takum_to_float(float_to_takum(1)) == 1);
    CHECK(takum_sub(float_to_takum(1),float_to_takum(1)) == 0);
    CHECK(takum_multiply(float_to_takum(1),float_to_takum(1)) == float_to_takum(1));
    CHECK(takum_add(float_to_takum(1),0) == float_to_takum(1));
}
static void vectors(void) {
    for (unsigned n = 8; n <= 32; n *= 2) {
        tk_sim s, before;
        tk_sim_init(&s);
        CHECK(tk_sim_execute(&s,0) == TK_SIM_BAD_CONFIG);
        CHECK(tk_sim_setvl(&s,999,n) == TK_SIM_OK);
        CHECK(s.vl == tk_sim_vlmax(n));
        long double full[TK_VLEN_BYTES];
        for (size_t i = 0; i < s.vl; ++i) full[i] = (long double)i - 5;
        CHECK(tk_sim_load(&s,1,full,s.vl) == TK_SIM_OK);
        CHECK(tk_sim_load(&s,2,full,s.vl) == TK_SIM_OK);
        before = s;
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTMUL,2,1,2)) == TK_SIM_OK);
        for (size_t i = 0; i < s.vl; ++i)
            CHECK(lane(&s,2,n,i) == tk_mul(lane(&before,1,n,i),lane(&before,2,n,i),n));
        tk_sim_init(&s);
        CHECK(tk_sim_setvl(&s,4,n) == TK_SIM_OK);
        const long double a[] = {1,2,-3,0.5L}, b[] = {4,-1,2,0};
        CHECK(tk_sim_load(&s,1,a,4) == TK_SIM_OK);
        CHECK(tk_sim_load(&s,2,b,4) == TK_SIM_OK);
        CHECK(tk_sim_load(&s,32,a,4) == TK_SIM_BAD_ARGUMENT);
        CHECK(tk_sim_load(&s,1,a,5) == TK_SIM_BAD_ARGUMENT);
        for (unsigned op = 0; op < 3; ++op) {
            memset(s.v[3],0x5a,TK_VLEN_BYTES);
            CHECK(tk_sim_execute(&s,tk_instruction((tk_opcode)op,3,1,2)) == TK_SIM_OK);
            for (size_t i = 0; i < 4; ++i) {
                uint32_t x = lane(&s,1,n,i), y = lane(&s,2,n,i);
                uint32_t expected = op == 0 ? tk_add(x,y,n) : op == 1 ? tk_sub(x,y,n) : tk_mul(x,y,n);
                CHECK(lane(&s,3,n,i) == expected);
            }
            for (size_t i = 4*n/8; i < TK_VLEN_BYTES; ++i) CHECK(s.v[3][i] == 0x5a);
        }
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTDOT,6,1,2)) == TK_SIM_OK);
        long double expected = 0;
        for (size_t i = 0; i < 4; ++i)
            expected += tk_to_ld(lane(&s,1,n,i),n) * tk_to_ld(lane(&s,2,n,i),n);
        CHECK(s.dot_accumulator == expected);
        CHECK(lane(&s,6,32,0) == tk_from_ld(expected,32));
        before = s;
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VNCVT,1,1,0)) ==
              (n == 8 ? TK_SIM_ILLEGAL : TK_SIM_OK));
        if (n != 8) {
            for (size_t i = 0; i < 4; ++i)
                CHECK(lane(&s,1,n/2,i) == tk_convert(lane(&before,1,n,i),n,n/2));
            s = before;
        } else CHECK(memcmp(&s,&before,sizeof(s)) == 0);
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTADD,1,1,2)) == TK_SIM_OK);
        for (size_t i = 0; i < 4; ++i)
            CHECK(lane(&s,1,n,i) == tk_add(lane(&before,1,n,i),lane(&before,2,n,i),n));
        before = s;
        CHECK(tk_sim_execute(&s,0xffffffff) == TK_SIM_ILLEGAL);
        CHECK(memcmp(&s,&before,sizeof(s)) == 0);
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTADD,3,1,2) | (1u << 12)) == TK_SIM_ILLEGAL);
        CHECK(memcmp(&s,&before,sizeof(s)) == 0);
        CHECK(tk_sim_setvl(&s,0,n) == TK_SIM_OK);
        before = s;
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTADD,3,1,2)) == TK_SIM_OK);
        CHECK(memcmp(s.v,before.v,sizeof(s.v)) == 0);
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTDOT,6,1,2)) == TK_SIM_OK);
        CHECK(lane(&s,6,32,0) == 0);
        CHECK(s.retired == before.retired + 2 && s.elements == before.elements);
        CHECK(s.pc == before.pc + 8);
        size_t done;
        uint32_t program[] = {tk_instruction(TK_VTADD,3,1,2),0xffffffff};
        CHECK(tk_sim_run(&s,program,2,&done) == TK_SIM_ILLEGAL && done == 1);
        CHECK(tk_sim_setvl(&s,1,n) == TK_SIM_OK);
        const long double invalid[] = {NAN};
        CHECK(tk_sim_load(&s,2,invalid,1) == TK_SIM_OK);
        CHECK(tk_sim_execute(&s,tk_instruction(TK_VTDOT,6,1,2)) == TK_SIM_OK);
        CHECK(lane(&s,6,32,0) == tk_nar(32));
        s.vl = 1000;
        CHECK(tk_sim_execute(&s,program[0]) == TK_SIM_BAD_CONFIG);
    }
    CHECK(tk_instruction(TK_VTADD,32,1,2) == UINT32_MAX);
    CHECK(tk_instruction(TK_VNCVT,3,1,2) == UINT32_MAX);
    CHECK(tk_sim_setvl(NULL,4,16) == TK_SIM_BAD_ARGUMENT);
}
static void randomized_vectors(void) {
    for (unsigned trial = 0; trial < 3000; ++trial) {
        unsigned n = 8u << (random_word()%3), op = random_word()%5;
        if (n == 8 && op == TK_VNCVT) op = TK_VTMUL;
        tk_sim s, expected;
        tk_sim_init(&s);
        CHECK(tk_sim_setvl(&s,random_word()%(tk_sim_vlmax(n)+1),n) == TK_SIM_OK);
        for (unsigned r = 0; r < TK_VREG_COUNT; ++r)
            for (unsigned i = 0; i < TK_VLEN_BYTES; ++i)
                s.v[r][i] = (uint8_t)random_word();
        unsigned a = random_word()%32, b = random_word()%32;
        unsigned d = trial%3 == 0 ? a : trial%3 == 1 ? b : random_word()%32;
        if (op == TK_VNCVT) b = 0;
        expected = s;
        long double sum = 0;
        for (size_t i = 0; i < s.vl; ++i) {
            uint32_t x = lane(&s,a,n,i), y = lane(&s,b,n,i), z = 0;
            unsigned width = op == TK_VNCVT ? n/2 : n;
            if (op == TK_VTDOT) {
                sum += tk_to_ld(x,n)*tk_to_ld(y,n);
                continue;
            }
            if (op == TK_VTADD) z = tk_add(x,y,n);
            else if (op == TK_VTSUB) z = tk_sub(x,y,n);
            else if (op == TK_VTMUL) z = tk_mul(x,y,n);
            else z = tk_convert(x,n,n/2);
            for (unsigned k = 0; k < width/8; ++k)
                expected.v[d][i*(width/8)+k] = (uint8_t)(z>>(8*k));
        }
        if (op == TK_VTDOT) {
            uint32_t z = tk_from_ld(sum,32);
            for (unsigned k = 0; k < 4; ++k) expected.v[d][k] = (uint8_t)(z>>(8*k));
        }
        CHECK(tk_sim_execute(&s,tk_instruction((tk_opcode)op,d,a,b)) == TK_SIM_OK);
        CHECK(memcmp(s.v,expected.v,sizeof(s.v)) == 0);
        CHECK(s.pc == 4 && s.retired == 1 && s.elements == s.vl);
        expected = s;
        /* Reserved funct3 values must leave the entire state unchanged. */
        CHECK(tk_sim_execute(&s,tk_instruction((tk_opcode)op,d,a,b) | (1u<<12)) == TK_SIM_ILLEGAL);
        CHECK(memcmp(&s,&expected,sizeof(s)) == 0);
    }
}
int main(void) {
    scalar(); vectors(); randomized_vectors();
    printf("OK: %lu verificacoes (T8/T16 exaustivos, 100000 palavras T32, ISA vetorial).\n",checks);
    return EXIT_SUCCESS;
}
