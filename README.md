# Takum-RVV — AOC 2026/2

Implementação em C de Takum logarítmico T8, T16 e T32 e um simulador funcional de instruções vetoriais.

## Implementado

- Extração dos campos S, D, R, C e F.
- Conversão, soma, subtração e multiplicação em T8/T16/T32.
- Tratamento de zero, NaR, saturação e arredondamento.
- Instruções `vtadd.vv`, `vtsub.vv`, `vtmul.vv`, `vtdot.vv` e `vncvt.t.t`.
- 32 registradores vetoriais de 256 bits, com SEW, VL e contador de programa.

## Executar

Requisitos: GCC, Make, libm e libquadmath. No Linux ou WSL:

```sh
make
./build/demo_scalar
./build/demo_rvv 16
```

A demonstração vetorial também aceita `8` e `32`.

## Testes

```sh
make test precision sanitize
```

Resultados:

- 903.584 verificações internas, incluindo 3.000 cenários vetoriais aleatórios.
- 1.328.424 comparações com uma referência em Float128.
- 1.427.969 comparações com a libtakum.

Para repetir a comparação com a libtakum:

```sh
git clone https://github.com/takum-arithmetic/libtakum.git ../reference-libtakum
git -C ../reference-libtakum checkout 776ac0ca7855f640100731ed8bbb0c7cf83b2d94
make oracle
```

O build padrão usa Float128 na codificação, soma e subtração. Para compilar sem quadmath, use `make clean` seguido de `make USE_FLOAT128=0`. Essa opção usa `long double` e tem precisão menor perto das fronteiras de arredondamento.

## Escopo e pendências

O simulador executa as cinco instruções propostas com codificação própria. Não executa a ISA RISC-V completa nem mede ciclos. O produto escalar acumula em `long double` e retorna T32.

Faltam os benchmarks de matrizes esparsas, a comparação com FP8/FP16/FP32, a análise de erros e de instruções/CSRs, o artigo e os slides.

Detalhes da arquitetura e dos testes em [docs/arquitetura.md](docs/arquitetura.md) e [docs/validacao.md](docs/validacao.md).
