# Arquitetura funcional

## Variante numérica

O projeto usa a representação Takum **logarítmica**.
Na libtakum essa família se chama `takum_log8/16/32`. A família chamada apenas
`takum8/16/32` na biblioteca atual é linear; comparar as duas como se fossem a
mesma codificação produziria resultados errados.

Para uma palavra finita e não nula:

```text
x = (-1)^S * exp(l / 2)
l = (-1)^S * (c + m)
m = F / 2^p
r = R, se D=1; r = 7-R, se D=0
c = bias[D*8 + R] + C_completo
```

`S` é o sinal, `D` a direção, `R` o campo de regime de três bits, `r` a largura
nominal da característica, `C` sua parte armazenada, `F` a fração e `p` a
quantidade de bits da fração. A tabela `C_BIAS` é indexada pelos bits DR brutos.

T16 e T32 têm `p = n-5-r`. Em T8 podem faltar bits da característica: os três
bits disponíveis depois de S/D/R são usados primeiro para C. Se `r>3`, os
`r-3` bits menos significativos ausentes de C são completados com zeros e F
fica vazio. Isso equivale a acrescentar zeros à direita para decodificar uma
palavra T8 como T32. Não há deslocamentos negativos.

Exemplos:

| Palavra | Valor/caso |
| --- | --- |
| T16 `0x0000` | Zero |
| T16 `0x8000` | NaR (Not a Real) |
| T16 `0x4000` | +1 |
| T16 `0xC000` | -1 |
| T16 `0x3800` | exp(-1/2) |
| T8 `0x01` | exp(-239/2), menor positivo |
| T8 `0x7F` | exp(239/2), maior positivo |

O sinal também entra no logaritmo decodificado. Uma palavra negativa não é
obtida apenas ligando o bit S de uma palavra positiva. A negação corresponde
a complemento de dois da palavra completa, salvo zero e NaR.

## Codificação e operações

`tk_from_ld` trata zero, NaN e infinitos antes de calcular `2*log(abs(x))`.
Números finitos fora do alcance saturam para o maior ou menor número não nulo
representável, preservando o sinal. NaN e infinitos viram NaR. Existe um único
zero. NaR se propaga inclusive em `0*NaR`.

`tk_from_log` calcula c/m, escolhe DR e compõe a palavra sem arredondar
separadamente o campo C. Um carry de arredondamento pode atravessar F/C/R/D;
a palavra inteira é arredondada ao mais próximo, com empate para palavra par.
O arredondamento é feito no espaço da codificação/logaritmo, não pela distância
linear entre números reais.

Na compilação padrão, soma e subtração decodificam, calculam e recodificam
com intermediários IEEE binary128 (`__float128`, libquadmath).
A API pública de entrada/saída usa `long double`. O backend alternativo
`USE_FLOAT128=0` usa `long double` e pode errar perto de fronteiras.
Multiplicação soma os logaritmos decodificados e combina os sinais, dispensando
as conversões exponenciais intermediárias. Este é um modelo funcional, não
uma ALU de portas lógicas. Funções matemáticas do hospedeiro podem influenciar
casos próximos de um limite de arredondamento.

`tk_convert` amplia por deslocamento à esquerda e reduz por arredondamento dos
bits descartados, com saturação que preserva zero/NaR. Assim, conversões entre
larguras não exigem passar por um float intermediário.

## Estado vetorial

- 32 registradores, cada um com VLEN=256 bits, armazenados em bytes little-endian.
- LMUL fixo em 1. SEW pode ser 8, 16 ou 32 bits.
- `tk_sim_setvl(avl,sew)` fixa `VL=min(AVL,VLEN/SEW)` e conta uma configuração.
- `tk_sim_load` inicializa lanes a partir de valores reais e quantiza para SEW.
- Cargas/configuração são APIs do hospedeiro, separadas das instruções customizadas.
- `pc` começa em zero e avança quatro bytes por instrução válida.
- `retired` conta instruções customizadas; `elements` soma VL por execução.
- Instruções ilegais não alteram registradores, PC ou contadores.
- Não há pipeline, latência, máscaras, memória de programa RISC-V ou execução da ISA base.

A rotina `tk_sim_run` recebe uma sequência de palavras de instrução, executa o
decodificador uma palavra por vez e interrompe na primeira instrução ilegal.
`completed` informa quantas instruções foram concluídas. Isso permite escrever
programas da extensão em vez de chamar somente operações em arrays.

## Codificação das instruções

Palavra de 32 bits com campos no estilo R e opcode local `0x0B` (custom-0).
Esta é uma escolha do grupo para o simulador, não uma codificação oficial RVV.

| Campo | Bits | Valor |
| --- | --- | --- |
| funct7 | 31:25 | Código da operação (0–4) |
| vs2 | 24:20 | Segundo registrador fonte; zero na conversão |
| vs1 | 19:15 | Primeiro registrador fonte |
| funct3 | 14:12 | Zero; demais valores são ilegais |
| vd | 11:7 | Registrador destino |
| opcode | 6:0 | `0x0B` |

| funct7 | Mnemônico | Semântica |
| --- | --- | --- |
| 0 | `vtadd.vv vd,vs1,vs2` | `vd[i] = TkSEW(vs1[i] + vs2[i])`, i < VL |
| 1 | `vtsub.vv vd,vs1,vs2` | `vd[i] = TkSEW(vs1[i] - vs2[i])`, i < VL |
| 2 | `vtmul.vv vd,vs1,vs2` | `vd[i] = TkSEW(vs1[i] * vs2[i])`, i < VL |
| 3 | `vtdot.vv vd,vs1,vs2` | Soma de produtos em long double; resultado T32 em vd[0] |
| 4 | `vncvt.t.t vd,vs1` | Converte VL elementos de SEW para SEW/2, SEW=16 ou 32 |

Os fontes são copiados antes da execução para permitir destinos que coincidem
com as fontes. Lanes de cauda e bytes não escritos são preservados.
`vncvt.t.t` não altera SEW/VL: para usar o destino em outra instrução aritmética,
configurar a nova SEW explicitamente. Seu campo vs2 é reservado e precisa ser
zero; T8→T4 é ilegal. Não é implementada conversão direta T32→T8 nesta instrução.

`vtdot.vv` multiplica os operandos quantizados em precisão do hospedeiro e não
quantiza cada produto para T8/T16/T32. A soma é sequencial em `long double` e só
a saída é convertida para T32. T32 não fica mais largo na saída, mas seu
acumulador possui mais precisão. Não é um acumulador exato nem Float128.
VL=0 produz zero no produto escalar; nas operações elemento a elemento nenhuma
lane é escrita. VL=0 ainda aposenta uma instrução válida.

## Exemplo de programa

```c
tk_sim s;
tk_sim_init(&s);
tk_sim_setvl(&s,4,16);
long double x[] = {1,2,-3,0.5L}, y[] = {4,-1,2,0};
tk_sim_load(&s,1,x,4);
tk_sim_load(&s,2,y,4);
uint32_t program[] = {
    tk_instruction(TK_VTADD,3,1,2),
    tk_instruction(TK_VTDOT,6,1,2)
};
size_t completed;
tk_sim_status status = tk_sim_run(&s,program,2,&completed);
/* Conferir status antes de consumir resultados. */
```

O produto escalar ideal desses vetores é -4. Com operandos T16, a saída T32 da demonstração é
aproximadamente -4,0020765 porque 2, 3 e 4 já foram quantizados para Takum antes
da execução. Nos benchmarks, o erro de quantização dos operandos deve ser
separado do erro de acumulação.
