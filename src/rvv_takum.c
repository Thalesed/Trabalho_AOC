#include "rvv_takum.h"
#include <math.h>
#include <string.h>

static uint32_t read_lane(const uint8_t *v, unsigned width, size_t i) {
    unsigned bytes = width / 8;
    uint32_t result = 0;
    for (unsigned b = 0; b < bytes; ++b)
        result |= (uint32_t)v[i * bytes + b] << (8 * b);
    return result;
}
static void write_lane(uint8_t *v, unsigned width, size_t i, uint32_t word) {
    unsigned bytes = width / 8;
    for (unsigned b = 0; b < bytes; ++b)
        v[i * bytes + b] = (uint8_t)(word >> (8 * b));
}
void tk_sim_init(tk_sim *s) {
    if (s) memset(s, 0, sizeof(*s));
}
size_t tk_sim_vlmax(unsigned n) { return tk_width_valid(n) ? TK_VLEN_BYTES * 8 / n : 0; }
tk_sim_status tk_sim_setvl(tk_sim *s, size_t avl, unsigned n) {
    if (!s || !tk_width_valid(n)) return TK_SIM_BAD_ARGUMENT;
    s->sew = n;
    size_t max = tk_sim_vlmax(n);
    s->vl = avl < max ? avl : max;
    ++s->configurations;
    return TK_SIM_OK;
}
tk_sim_status tk_sim_load(tk_sim *s, unsigned reg, const long double *x, size_t count) {
    if (!s || reg >= TK_VREG_COUNT || (!x && count)) return TK_SIM_BAD_ARGUMENT;
    if (!tk_width_valid(s->sew) || s->vl > tk_sim_vlmax(s->sew)) return TK_SIM_BAD_CONFIG;
    if (count > s->vl) return TK_SIM_BAD_ARGUMENT;
    for (size_t i = 0; i < count; ++i)
        write_lane(s->v[reg], s->sew, i, tk_from_ld(x[i],s->sew));
    return TK_SIM_OK;
}
tk_sim_status tk_sim_read(const tk_sim *s, unsigned reg, unsigned width,
                          size_t lane, uint32_t *word) {
    if (!s || !word || reg >= TK_VREG_COUNT || !tk_width_valid(width) ||
        lane >= tk_sim_vlmax(width)) return TK_SIM_BAD_ARGUMENT;
    *word = read_lane(s->v[reg], width, lane);
    return TK_SIM_OK;
}
uint32_t tk_instruction(tk_opcode op, unsigned vd, unsigned vs1, unsigned vs2) {
    if ((unsigned)op > TK_VNCVT || vd >= 32 || vs1 >= 32 || vs2 >= 32 ||
        (op == TK_VNCVT && vs2 != 0)) return UINT32_MAX;
    return ((uint32_t)op << 25) | (vs2 << 20) | (vs1 << 15) | (vd << 7) | TK_CUSTOM_OPCODE;
}
tk_sim_status tk_sim_execute(tk_sim *s, uint32_t insn) {
    if (!s) return TK_SIM_BAD_ARGUMENT;
    if (!tk_width_valid(s->sew) || s->vl > tk_sim_vlmax(s->sew)) return TK_SIM_BAD_CONFIG;
    unsigned op = insn >> 25, vd = (insn >> 7) & 31;
    unsigned vs1 = (insn >> 15) & 31, vs2 = (insn >> 20) & 31;
    if ((insn & 127) != TK_CUSTOM_OPCODE || ((insn >> 12) & 7) != 0 ||
        op > TK_VNCVT || (op == TK_VNCVT && (vs2 != 0 || s->sew == 8)))
        return TK_SIM_ILLEGAL;
    /* Snapshot garante vd == vs1/vs2 inclusive na conversao estreita. */
    uint8_t a[TK_VLEN_BYTES], b[TK_VLEN_BYTES];
    memcpy(a,s->v[vs1],sizeof(a)); memcpy(b,s->v[vs2],sizeof(b));
    if (op == TK_VTDOT) {
        long double sum = 0;
        for (size_t i = 0; i < s->vl; ++i)
            sum += tk_to_ld(read_lane(a,s->sew,i),s->sew) *
                   tk_to_ld(read_lane(b,s->sew,i),s->sew);
        s->dot_accumulator = sum;
        /* Resultado T32 na lane 0; demais bytes preservados. vl=0 gera zero. */
        write_lane(s->v[vd],32,0,tk_from_ld(sum,32));
    } else {
        unsigned out_width = op == TK_VNCVT ? s->sew / 2 : s->sew;
        for (size_t i = 0; i < s->vl; ++i) {
            uint32_t x = read_lane(a,s->sew,i), y = read_lane(b,s->sew,i), z;
            switch (op) {
            case TK_VTADD: z = tk_add(x,y,s->sew); break;
            case TK_VTSUB: z = tk_sub(x,y,s->sew); break;
            case TK_VTMUL: z = tk_mul(x,y,s->sew); break;
            default: z = tk_convert(x,s->sew,out_width); break;
            }
            write_lane(s->v[vd],out_width,i,z);
        }
    }
    ++s->retired; s->elements += s->vl; s->pc += 4;
    return TK_SIM_OK;
}
tk_sim_status tk_sim_run(tk_sim *s, const uint32_t *p, size_t count, size_t *done) {
    if (done) *done = 0;
    if (!s || (!p && count)) return TK_SIM_BAD_ARGUMENT;
    for (size_t i = 0; i < count; ++i) {
        tk_sim_status status = tk_sim_execute(s,p[i]);
        if (status != TK_SIM_OK) return status;
        if (done) *done = i + 1;
    }
    return TK_SIM_OK;
}
const char *tk_sim_error(tk_sim_status status) {
    switch (status) {
    case TK_SIM_OK: return "sucesso";
    case TK_SIM_BAD_ARGUMENT: return "argumento invalido";
    case TK_SIM_BAD_CONFIG: return "configuracao vetorial invalida";
    default: return "instrucao ilegal";
    }
}
