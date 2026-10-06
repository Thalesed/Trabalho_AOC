# Continuidade do trabalho do Thales

## Base e correções

Partimos de `Thalesed/Trabalho_AOC`, branch main, commit
`cefd935cb2c79c17a024ed53b03308f87336789b`. O protótipo tinha a tabela de bias,
conversões, soma/subtração/multiplicação e exemplos T16.

As funções originais foram preservadas como wrappers de uma implementação
comum de largura parametrizada. A tabela continua em `Takum.c`, agora com
indexação direta pelos bits DR.

Correções relevantes:

1. Declarações de funções no header para evitar uso antes da declaração.
2. Separação entre R armazenado e r efetivo, particularmente quando D=0.
3. Base logarítmica consistente: `exp(l/2)` e `2*log(abs(x))`.
4. Aplicação do sinal também ao logaritmo; codificação negativa correta.
5. Decodificação explícita de zero e NaR antes da equação normal.
6. Saturação correta, sem usar `0xFFFF` como marcador genérico de overflow.
7. Arredondamento com carry pela palavra completa e empate para par.
8. Decodificação T8 com característica truncada e completação por zeros.
9. Intermediários Float128 na codificação, soma e subtração do build padrão;
   multiplicação no domínio logarítmico. A acumulação vetorial permanece long double.

Essas correções alteram padrões binários e resultados do protótipo. Manter o
nome de uma função não torna as antigas palavras Takum compatíveis com a
codificação corrigida.

## Primeira entrega — 11/10/2026

O cronograma do roteiro marca a ALU básica T8/T16/T32 para essa data. O código
agora inclui essas larguras e as cinco operações vetoriais propostas, com
validação reproduzível. Para demonstrar:

```sh
make test
./build/demo_scalar
./build/demo_rvv 16
```

O relatório do simulador distingue instruções customizadas e configurações;
o contador de elementos não representa instruções escalares retiradas nem
ciclos. Ainda não há experimento que demonstre redução de instruções ou CSRs.

O roteiro dá C++/Python como exemplos de simulador; esta continuação adota C11
conforme a abordagem em C apresentada pelo grupo. É uma implementação própria
em nível de instrução da extensão, sem dependência de um simulador externo.

## Integrar as alterações

A continuação fica na branch `feat/takum-rvv-c` do repositório do Thales.
O envio utiliza a conta lucaspimentab pelo navegador; a integração API retornou
403 mesmo com permissão de colaborador. Para usar a versão revisada:

```sh
git clone https://github.com/Thalesed/Trabalho_AOC.git
cd Trabalho_AOC
git switch feat/takum-rvv-c
make
make test precision
./build/demo_rvv 16
```

A main permanece como base original até a integração da branch pelo grupo.

## Trabalho seguinte

1. Definir datasets/matrizes esparsas e kernel CSR de matriz-vetor.
2. Fixar formatos baseline: identificar variante FP8 explicitamente, além de FP16/FP32.
3. Implementar a referência Float128 para o benchmark (a referência atual valida
   o codec e a aritmética, não executa matrizes esparsas).
4. Medir erros relativos; tratar referência zero com erro absoluto separado.
5. Contar operações e instruções com a mesma política de acumulação em cada baseline.
6. Estudar quais CSRs são necessários na proposta e comparar com baseline real.
7. Escrever artigo IEEE em inglês, até seis páginas, e apresentação de 15 minutos.

Converter este C automaticamente para VHDL não basta para obter RTL
sintetizável. `logl`/`expl`, acumulador e controle vetorial precisariam de uma
arquitetura de hardware e validação específica. Para a continuação atual,
o entregável é o simulador funcional em C.
