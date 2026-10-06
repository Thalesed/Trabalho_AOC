#include "../src/rvv_takum.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

static void show(const tk_sim *s, const char *label, unsigned reg, unsigned width, size_t lanes) {
    printf("%s (T%u): [",label,width);
    for (size_t i = 0; i < lanes; ++i) {
        uint32_t word = 0;
        if (tk_sim_read(s,reg,width,i,&word) != TK_SIM_OK) exit(EXIT_FAILURE);
        printf("%s%.7Lf",i ? ", " : "",tk_to_ld(word,width));
    }
    puts("]");
}
int main(int argc, char **argv) {
    unsigned width = 16;
    if (argc > 2) return EXIT_FAILURE;
    if (argc == 2) {
        char *end;
        unsigned long arg = strtoul(argv[1],&end,10);
        if (!*argv[1] || *end || (arg != 8 && arg != 16 && arg != 32)) {
            fputs("Uso: demo_rvv [8|16|32]\n",stderr); return EXIT_FAILURE;
        }
        width = (unsigned)arg;
    }
    tk_sim s;
    tk_sim_init(&s);
    const long double a[] = {1,2,-3,0.5L}, b[] = {4,-1,2,0};
    if (tk_sim_setvl(&s,4,width) != TK_SIM_OK ||
        tk_sim_load(&s,1,a,4) != TK_SIM_OK || tk_sim_load(&s,2,b,4) != TK_SIM_OK)
        return EXIT_FAILURE;
    uint32_t program[] = {
        tk_instruction(TK_VTADD,3,1,2), tk_instruction(TK_VTSUB,4,1,2),
        tk_instruction(TK_VTMUL,5,1,2), tk_instruction(TK_VTDOT,6,1,2),
        tk_instruction(TK_VNCVT,7,1,0)
    };
    size_t done;
    tk_sim_status status = tk_sim_run(&s,program,width == 8 ? 4 : 5,&done);
    if (status != TK_SIM_OK) { fprintf(stderr,"%s\n",tk_sim_error(status)); return EXIT_FAILURE; }
    show(&s,"v1",1,width,4); show(&s,"v2",2,width,4);
    show(&s,"vtadd.vv v3,v1,v2",3,width,4);
    show(&s,"vtsub.vv v4,v1,v2",4,width,4);
    show(&s,"vtmul.vv v5,v1,v2",5,width,4);
    show(&s,"vtdot.vv v6,v1,v2",6,32,1);
    if (width > 8) show(&s,"vncvt.t.t v7,v1",7,width/2,4);
    printf("Acumulador antes da quantizacao T32: %.12Lf\n",s.dot_accumulator);
    printf("PC=%" PRIu64 "; instrucoes=%" PRIu64 "; elementos=%" PRIu64 "\n",
           s.pc,s.retired,s.elements);
    return EXIT_SUCCESS;
}
