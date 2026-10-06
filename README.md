# Trabalho AOC — Takum e extensão vetorial em C

Continuação do protótipo `Takum.c` do Thales, baseada no commit
`cefd935cb2c79c17a024ed53b03308f87336789b`.

Implementação didática em **C11**, com aritmética **Takum logarítmica** T8/T16/T32
(fase 1) e simulador funcional de uma extensão vetorial proposta (fase 2).
O artigo citado nas diretrizes descreve Takum logarítmico. A libtakum atual
separa essa variante (`takum_log*`) da variante linear (`takum*`); os testes
usam explicitamente a primeira.

## Executar

A compilação principal precisa de GCC, Make, libm e libquadmath (Float128).
Linux, WSL ou ambiente equivalente com GCC são suficientes.

```sh
make
make test
./build/demo_scalar
./build/demo_rvv 16
```

Também funcionam `./build/demo_rvv 8` e `./build/demo_rvv 32`.
`make demo` executa todas as demonstrações.

O caminho simples abaixo usa intermediários long double e tem precisão menor
perto das fronteiras de arredondamento; para validação use o Makefile padrão:

```sh
gcc -std=c11 -Wall -Wextra Takum.c -lm -o takum_demo
./takum_demo
```

Também existe `make clean && make USE_FLOAT128=0` para compilar sem quadmath.
Execute `make clean` ao trocar o backend numérico.

Para usar o codec como biblioteca, compilar com `-DTAKUM_NO_MAIN` e incluir
`takum.h`. O Makefile já faz isso no simulador e nos testes.

## Implementado

| Etapa | Funcionalidade |
| --- | --- |
| Fase 1 | Extração S/D/R/C/F; conversões; soma, subtração e multiplicação T8/T16/T32 |
| Fase 1 | Zero, NaR, saturação e arredondamento ao mais próximo com empate para par |
| Fase 2 | 32 registradores de 256 bits, SEW, VL, PC e decodificador de palavras de 32 bits |
| Fase 2 | `vtadd.vv`, `vtsub.vv`, `vtmul.vv`, `vtdot.vv`, `vncvt.t.t` |
| Validação | Testes próprios e comparação opcional com a libtakum |

As funções T16 originais (`float_to_takum`, `takum_to_float`, `takum_add`,
`takum_sub`, `takum_multiply`) permanecem disponíveis. As palavras binárias
mudam porque foram corrigidas para o formato logarítmico de referência.

## Testes

```sh
make test
make precision
make sanitize
```

`make test` inclui todos os 256 padrões T8, todos os 65.536 padrões T16 e
100.000 padrões T32, além de testes do simulador. `make sanitize` verifica
acessos à memória e comportamento indefinido. O alvo desativa apenas a busca
de vazamentos do LeakSanitizer, incompatível com o ambiente de execução usado
na validação; AddressSanitizer e UndefinedBehaviorSanitizer permanecem ativos.
O código não usa alocação dinâmica.

A libtakum **não é uma dependência de execução**. Para repetir a validação
independente, baixe a referência fora deste repositório:

```sh
git clone https://github.com/takum-arithmetic/libtakum.git ../reference-libtakum
git -C ../reference-libtakum checkout 776ac0ca7855f640100731ed8bbb0c7cf83b2d94
make oracle LIBTAKUM=../reference-libtakum
```

Resultados obtidos em 06/10/2026 com GCC/Linux:

- 903.584 verificações internas aprovadas, incluindo 3.000 cenários vetoriais aleatórios.
- 1.328.424 comparações com decodificador/quantizador independente em Float128.
- 1.427.969 comparações com `takum_log8/16/32` aprovadas.
- Compilação com `-Wall -Wextra -Wpedantic -Werror` aprovada.
- AddressSanitizer/UndefinedBehaviorSanitizer aprovados.

Os resultados são de correção funcional, não de desempenho ou redução de área.
Escopo exato e tolerâncias: [docs/validacao.md](docs/validacao.md).

## Organização

| Arquivo | Papel |
| --- | --- |
| `Takum.c` | Codec, aritmética escalar, API T16 original e demonstração escalar |
| `takum.h` | API T8/T16/T32 e campos decodificados |
| `src/rvv_takum.c` e `.h` | Estado, codificação e execução de instruções vetoriais |
| `examples/demo_rvv.c` | Programa de demonstração das cinco operações |
| `tests/test_takum.c` | Testes internos |
| `tests/test_oracle.c` | Comparação independente com libtakum |
| `tests/test_float128.c` | Fronteiras de arredondamento e soma/subtração em Float128 |
| `docs/arquitetura.md` | Equações, semântica e codificação da proposta |
| `docs/continuidade.md` | Alterações em relação ao protótipo e próxima entrega |

## Limites e próximas etapas

É um simulador funcional em nível de instrução, com codificação **local da
proposta Zvtakum**. Não é um processador RISC-V completo nem uma integração
com Spike/QEMU; não estima ciclos, área ou potência. VLEN=256, LMUL=1,
sem máscaras e sem `vstart`. A configuração e as cargas usam APIs do hospedeiro.

A acumulação de `vtdot.vv` usa `long double` e produz T32. `long double` depende
da plataforma e **não é necessariamente Float128**. Essa escolha precisa ser
mantida explícita nas medições e no artigo.

A fase 3 ainda requer benchmark de matriz esparsa, baseline FP8/FP16/FP32,
referência Float128 e análise de erros e contagem de instruções/CSRs. Relatório
IEEE em inglês e slides também continuam pendentes. O PDF de HUB é uma
referência de outro formato, e não foi usado para definir a codificação Takum.

Referências: [artigo Takum](https://arxiv.org/abs/2404.18603),
[libtakum](https://github.com/takum-arithmetic/libtakum).
