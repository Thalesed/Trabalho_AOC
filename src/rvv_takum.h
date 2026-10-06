#ifndef AOC_RVV_TAKUM_H
#define AOC_RVV_TAKUM_H
#include <stddef.h>
#include <stdint.h>
#include "../takum.h"
#define TK_VLEN_BYTES 32u /* VLEN = 256 bits, LMUL = 1 */
#define TK_VREG_COUNT 32u
#define TK_CUSTOM_OPCODE 0x0bu

typedef enum {
    TK_VTADD = 0, TK_VTSUB = 1, TK_VTMUL = 2, TK_VTDOT = 3, TK_VNCVT = 4
} tk_opcode;
typedef enum {
    TK_SIM_OK = 0, TK_SIM_BAD_ARGUMENT, TK_SIM_BAD_CONFIG, TK_SIM_ILLEGAL
} tk_sim_status;
typedef struct {
    uint8_t v[TK_VREG_COUNT][TK_VLEN_BYTES];
    unsigned sew;
    size_t vl;
    uint64_t pc, retired, elements, configurations;
    long double dot_accumulator;
} tk_sim;

void tk_sim_init(tk_sim *sim);
size_t tk_sim_vlmax(unsigned sew);
tk_sim_status tk_sim_setvl(tk_sim *sim, size_t avl, unsigned sew);
tk_sim_status tk_sim_load(tk_sim *sim, unsigned reg, const long double *x, size_t count);
tk_sim_status tk_sim_read(const tk_sim *sim, unsigned reg, unsigned width,
                          size_t lane, uint32_t *word);
/* Codificacao didatica local; UINT32_MAX sinaliza argumento invalido. */
uint32_t tk_instruction(tk_opcode op, unsigned vd, unsigned vs1, unsigned vs2);
tk_sim_status tk_sim_execute(tk_sim *sim, uint32_t instruction);
tk_sim_status tk_sim_run(tk_sim *sim, const uint32_t *program, size_t count,
                         size_t *completed);
const char *tk_sim_error(tk_sim_status status);
#endif
