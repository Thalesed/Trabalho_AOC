# Validação em 06/10/2026

Ambiente: Linux, GCC, C11; versão da libtakum
`776ac0ca7855f640100731ed8bbb0c7cf83b2d94`.

## Testes internos

Comando: `make test`. Resultado: **903.584 verificações aprovadas**.

- Todos os padrões de T8 e T16: ida/volta, identidade multiplicativa,
  identidade aditiva, negação e cancelamento (exceto NaR).
- 100.000 padrões T32 com gerador xorshift e seed fixo.
- Campos S/D/R/r/C/F, incluindo característica incompleta de T8.
- Zero, -0, NaR, NaN, infinitos, saturação positiva/negativa.
- Conversões, empates pares/ímpares e preservação de zero/NaR.
- Execução vetorial em três SEWs, uso de todas as lanes até VLMAX, redução de largura,
  destino sobreposto à fonte, preservação da cauda e VL=0.
- Produto escalar com soma ampliada e propagação de NaR.
- Rejeição de instruções/configurações inválidas sem alteração de estado.
- Execução de programa interrompida na primeira instrução ilegal.
- 3.000 cenários vetoriais aleatórios, comparando todos os bytes dos 32
  registradores; destinos sobrepostos e instruções ilegais sem efeitos.

Comando: `make sanitize`. Mesmos testes aprovados com AddressSanitizer e
UndefinedBehaviorSanitizer. LeakSanitizer está desativado no alvo porque não
consegue inspecionar processos neste ambiente; não existe malloc/free no projeto.

## Comparação independente

Comando: `make oracle LIBTAKUM=../reference-libtakum`.
Resultado: **1.427.969 comparações aprovadas**.

| Comparação | Cobertura | Critério |
| --- | --- | --- |
| Decodificação T8/T16 | Todas as palavras | Zero/NaR exatos; valores finitos com erro relativo < 3e-14 |
| Decodificação T32 | 100.000 palavras | Mesmo critério |
| Codificação a partir de double | 100.000 valores por largura, mais especiais | Palavra binária idêntica |
| Soma/subtração/multiplicação T8 | Todos os 65.536 pares em cada operação | Palavra binária idêntica |
| Soma/subtração/multiplicação T16/T32 | 100.000 pares por largura em cada operação | Palavra binária idêntica |
| Conversão T16→T8 | Todas as palavras | Palavra binária idêntica |
| Conversão T32→T16 | 100.000 palavras | Palavra binária idêntica |

Valores de codificação são amostrados como double de sinais variados e
logaritmo entre -300 e +300; incluem regiões de saturação. Seeds fixos no código.

A tolerância da decodificação considera a diferença entre `exp(l/2)` do projeto
e a potência da constante √e usada pela biblioteca. Nos extremos, essa
constante de precisão finita amplifica a diferença. Não há tolerância para
comparação de palavras produzidas pela codificação/operações/conversões.

Não se trata de prova formal de arredondamento correto para todo valor real,
de validação exaustiva dos 2^32 padrões T32, nem de benchmark do projeto.
Testes por amostragem não excluem diferenças em casos não cobertos, sobretudo
perto de limites de arredondamento e em outra implementação de libm/long double.

## Revisão independente de 06/10/2026

O ZIP recebido continha a versão long double, sem os intermediários Float128
mencionados na conversa. Esta revisão adicionou libquadmath ao build padrão.
A entrada pública long double é promovida antes de logq; soma/subtração usam
expq/logq sem reduzir o intermediário para long double. Multiplicação continua
somando logaritmos exatos no domínio binário do formato.

`make precision`: **1.328.424 comparações aprovadas**. O teste tem decoder próprio
baseado na palavra expandida T32 e codificador por busca binária das fronteiras
logarítmicas. Não chama tk_extract ou tk_from_log para gerar a referência.
Cobre todas as fronteiras positivas internas T8/T16 e 100.000 fronteiras T32,
com o long double central e seus vizinhos nextafterl, nos dois sinais. Também
cobre todos os pares T8 em soma/subtração e 100.000 pares T16/T32.

O mesmo teste aplicado ao código antigo falha já na fronteira T8 entre 0x01 e
0x02: um valor imediatamente abaixo da fronteira foi arredondado para 0x02 em
vez de 0x01. A versão Float128 passa esse caso e os demais casos acima.

Resultados desta sessão: **3.659.977 verificações/comparações**, somando test,
oracle e precision. Esse total conta assertivas e casos comparados, não entradas
únicas nem uma prova exaustiva da aritmética T16/T32. As execuções dos sanitizadores
não são contadas novamente. GCC -fanalyzer passou no codec e no simulador.

A referência independente e o backend principal compartilham libquadmath;
portanto não são independentes da biblioteca transcendental. Arredondamento
correto para todas as entradas e exatidão do produto escalar não são garantidos.
O produto escalar ainda acumula em long double sequencialmente. A codificação
local, VL=0 na redução e a conversão de largura têm semântica própria documentada;
esses testes não comprovam conformidade integral com RVV.
