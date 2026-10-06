# Testes

Os testes foram executados com GCC, C11, `-Wall -Wextra -Wpedantic -Werror`, AddressSanitizer e UndefinedBehaviorSanitizer.

`make test` executa 903.584 verificações. O conjunto cobre todos os padrões T8 e T16, 100.000 padrões T32, conversões, operações aritméticas, redução de largura, zero, NaR, saturação, arredondamento, sobreposição de registradores, `VL`, `SEW`, cauda dos vetores, instruções inválidas e sequências de instruções.

`make precision` executa 1.328.424 comparações com um decodificador e quantizador independentes em Float128, incluindo valores próximos às fronteiras de arredondamento.

`make oracle` executa 1.427.969 comparações com `takum_log8`, `takum_log16` e `takum_log32` da libtakum:

```sh
git clone https://github.com/takum-arithmetic/libtakum.git ../reference-libtakum
git -C ../reference-libtakum checkout 776ac0ca7855f640100731ed8bbb0c7cf83b2d94
make oracle LIBTAKUM=../reference-libtakum
```

As palavras produzidas pela codificação, pelas operações e pelas conversões coincidem com a referência nos casos comparados. A decodificação usa tolerância numérica porque as bibliotecas podem usar constantes internas diferentes.
