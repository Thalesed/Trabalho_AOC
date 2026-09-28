#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

void printBinary(int num) {
  for (int i = sizeof(int) * 8 - 1; i >= 0; i--) {
    printf("%d", (num >> i) & 1);
  }
  printf("\n");
}

typedef uint16_t takum16_t;

static const int16_t C_BIAS[] = {
  [0] = -255, /* 0 000 → D=0, r=7, -2^(r+1)+1 */
  [1] = -127, /* 0 001 → D=0, r=6, -2^(r+1)+1 */
  [2] = -63,  /* 0 010 → D=0, r=5, -2^(r+1)+1 */
  [3] = -31,  /* 0 011 → D=0, r=4, -2^(r+1)+1 */
  [4] = -15,  /* 0 100 → D=0, r=3, -2^(r+1)+1 */
  [5] = -7,   /* 0 101 → D=0, r=2, -2^(r+1)+1 */
  [6] = -3,   /* 0 110 → D=0, r=1, -2^(r+1)+1 */
  [7] = -1,   /* 0 111 → D=0, r=0, -2^(r+1)+1 */
  [8] = 0,    /* 1 000 → D=1, r=0, 2^r-1 */
  [9] = 1,    /* 1 001 → D=1, r=1, 2^r-1 */
  [10] = 3,   /* 1 010 → D=1, r=2, 2^r-1 */
  [11] = 7,   /* 1 011 → D=1, r=3, 2^r-1 */
  [12] = 15,  /* 1 100 → D=1, r=4, 2^r-1 */
  [13] = 31,  /* 1 101 → D=1, r=5, 2^r-1 */
  [14] = 63,  /* 1 110 → D=1, r=6, 2^r-1 */
  [15] = 127, /* 1 111 → D=1, r=7, 2^r-1 */
};


takum16_t takum_add(takum16_t n1, takum16_t n2){
    float f1 = takum_to_float(n1);
    float f2 = takum_to_float(n2);
    takum16_t resultado = float_to_takum(f1 + f2);

    return resultado;
}


takum16_t takum_sub(takum16_t n1, takum16_t n2){
    float f1 = takum_to_float(n1);
    float f2 = takum_to_float(n2);
    takum16_t resultado = float_to_takum(f1 - f2);

    return resultado;
}


float takum_to_float(takum16_t num){
  uint16_t tmp  = num >> 15; // desloca 15 bits: isola bit 15
  uint16_t S = tmp; // S: 1 bit
  tmp = num >> 14; // desloca 14 bits: isola bits 15..14
  tmp -= 2*S; // subtrai 2*S: remove bit 15, sobra D no bit 0
  uint16_t D = tmp; // D: 1 bit
  tmp = num >> 11; // desloca 11 bits: isola bits 15..11
  tmp -= (S*16 + D*8); // subtrai S*16 e D*8: remove S e D, sobra R
  uint16_t R = tmp; // R: 3 bits
  tmp = num >> (11 - R); // desloca (11-R) bits: alinha S,D,R,C no topo
  tmp -= S*(1<<(R+4)); // subtrai S deslocado: remove S
  tmp -= D*(1<<(R+3)); // subtrai D deslocado: remove D
  uint16_t tmp2 = R >> 2; // desloca 2 bits: isola bit 2 de R
  tmp -= tmp2*(1<<(R+2)); // subtrai bit 2 de R deslocado: remove
  uint16_t tmp3 = (R >> 1) - tmp2*2; // desloca 1 bit e subtrai: isola bit 1 de R
  tmp -= tmp3*(1<<(R+1)); // subtrai bit 1 de R deslocado: remove
  uint16_t tmp4 = R - tmp2*4 - tmp3*2; // isola bit 0 de R
  tmp -= tmp4*(1<<R); // subtrai bit 0 de R deslocado: remove
  uint16_t C = tmp; // C: R bits
  tmp = num - S*(1<<15); // subtrai S no bit 15: zera S
  tmp -= D*(1<<14); // subtrai D no bit 14: zera D
  tmp -= (tmp2*(1<<13) + tmp3*(1<<12) + tmp4*(1<<11)); // subtrai R: zera R
  tmp2 = C; // copia C para tmp2
  for(uint16_t i=0;i<R;i++){ // itera R vezes
    tmp3 = tmp2 >> (R-1-i); // desloca (R-1-i) bits: isola bit i de C
    tmp -= tmp3*(1<<(10-i)); // subtrai bit i de C na posicao 10-i
    tmp2 -= tmp3*(1<<(R-1-i)); // remove bit i de C de tmp2
  }
  uint16_t F = tmp; // F: 11-R bits

  float resultado;  // 32 bits

  int16_t c; // 16 bits

  resultado = (F*1.0) / (1<<(11-R));  // F / 2^(11-R)
  
  if(D == 0){
    R  = 7 - R; // inverte R: 0..7 vira 7..0
  }
  c = C + C_BIAS[R+D*(1<<3)]; // soma C com C_BIAS[R + D*8]
  resultado += c;                      // soma c

  resultado     = pow(2, resultado); // 2^resultado
  if(S == 1){
    resultado  *= (-1); // muda o sinal se S=1
  }

  return resultado;
}


takum16_t float_to_takum(float x){
  if(x == 0.0){
    return 0x0000;  // zero: todos bits 0
  }
  
  uint16_t S    = (x < 0) ? 1 : 0; // S: 1 bit
  float abs_x   = fabs(x);  // 32 bits
  
  float l = log2(abs_x); // 32 bits
  
  uint16_t D = (abs_x >= 1.0) ? 1 : 0; // D: 1 bit
  
  int16_t c = (int16_t)floor(l); // 16 bits
  float m = l - c;  // 32 bits
  
  uint16_t R_adj = 0; // 16 bits
  uint16_t R = 0; // 16 bits
  uint16_t C = 0; // 16 bits
  bool found = false; // 1 bit
  
  for(R_adj = 0; R_adj <= 7; R_adj++){ // itera 0..7
    int16_t bias = C_BIAS[R_adj + D*8]; // 16 bits
    
    if(D == 1){
      R  = R_adj;
    } else {
      R = 7 - R_adj; 
    }
    
    int16_t C_candidate = c - bias; // 16 bits
    int16_t max_C = (1 << R) - 1;   // 16 bits
    
    if(C_candidate >= 0 && C_candidate <= max_C){
      C = (uint16_t)C_candidate; // 16 bits
      found = true; // 1 bit
      break;
    }
  }
  
  if(!found){
    return 0xFFFF; // overflow: todos bits 1
  }
  
  uint16_t F = (uint16_t)round(m * (1 << (11 - R))); // 16 bits
  
  if(F >= (1 << (11 - R))){
    F = 0; // zera F
    C++; // incrementa C
  }
  
  uint16_t num = 0;                       // 16 bits
  num |= (S & 1) << 15; // bit 15 = S
  num |= (D & 1) << 14; // bit 14 = D
  
  uint16_t R_field; // 16 bits
  if(D == 1){
    R_field = R_adj; // R_field = R_adj
  } else {
    R_field = 7 - R_adj; // R_field = 7 - R_adj
  }
  num |= (R_field & 7) << 11; // bits 13..11 = R_field
  num |= (C & ((1 << R) - 1)) << (11 - R); // bits 10..(11-R) = C
  num |= (F & ((1 << (11 - R)) - 1)); // bits (10-R)..0 = F
  
  return num;
}


takum16_t takum_multiply(takum16_t a, takum16_t b){
  uint16_t sign_a = (a >> 15) & 1; // isola bit 15 de a
  uint16_t sign_b = (b >> 15) & 1; // isola bit 15 de b
  uint16_t result_sign = sign_a ^ sign_b; // XOR: sinal do resultado
  
  float fa = takum_to_float(a); // 32 bits
  float fb = takum_to_float(b); // 32 bits
  
  float l_a = 2.0 * log(fabs(fa)); // 32 bits
  float l_b = 2.0 * log(fabs(fb)); // 32 bits
  
  float l_result = l_a + l_b; // 32 bits
  
  float f_result = exp(l_result / 2.0);   // 32 bits
  if(result_sign){
    f_result *= (-1); // muda o sinal se result_sign=1
  }
  return float_to_takum(f_result);
}

int main(){
  takum16_t a, b, r;
  float fa, fb, fr;

  printf(" takum_to_float / float_to_takum \n");
  float valores[] = {0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 3.14159, 100.0, -100.0, 1e-10, 1e10};
  int n = sizeof(valores)/sizeof(valores[0]);
  for(int i=0;i<n;i++){
    a  = float_to_takum(valores[i]);
    fa = takum_to_float(a);
    printf("%f -> 0x%04X -> %f\n", valores[i], a, fa);
  }

  printf("\n takum_add \n");
  a  = float_to_takum(2.5);
  b  = float_to_takum(3.0);
  r  = takum_add(a, b);
  fr = takum_to_float(r);
  printf("2.5 + 3.0 = %f (esperado 5.5)\n", fr);

  a  = float_to_takum(-4.0);
  b  = float_to_takum(1.5);
  r  = takum_add(a, b);
  fr = takum_to_float(r);
  printf("-4.0 + 1.5 = %f (esperado -2.5)\n", fr);

  a  = float_to_takum(0.25);
  b  = float_to_takum(0.75);
  r  = takum_add(a, b);
  fr = takum_to_float(r);
  printf("0.25 + 0.75 = %f (esperado 1.0)\n", fr);

  printf("\n takum_sub \n");
  a  = float_to_takum(5.0);
  b  = float_to_takum(3.0);
  r  = takum_sub(a, b);
  fr = takum_to_float(r);
  printf("5.0 - 3.0 = %f (esperado 2.0)\n", fr);

  a  = float_to_takum(2.0);
  b  = float_to_takum(7.0);
  r  = takum_sub(a, b);
  fr = takum_to_float(r);
  printf("2.0 - 7.0 = %f (esperado -5.0)\n", fr);

  printf("\n takum_multiply \n");
  a  = float_to_takum(2.5);
  b  = float_to_takum(3.0);
  r  = takum_multiply(a, b);
  fr = takum_to_float(r);
  printf("2.5 * 3.0 = %f (esperado 7.5)\n", fr);

  a  = float_to_takum(-2.0);
  b  = float_to_takum(4.0);
  r  = takum_multiply(a, b);
  fr = takum_to_float(r);
  printf("-2.0 * 4.0 = %f (esperado -8.0)\n", fr);

  a  = float_to_takum(-3.0);
  b  = float_to_takum(-3.0);
  r  = takum_multiply(a, b);
  fr = takum_to_float(r);
  printf("-3.0 * -3.0 = %f (esperado 9.0)\n", fr);

  a  = float_to_takum(0.5);
  b  = float_to_takum(0.5);
  r  = takum_multiply(a, b);
  fr = takum_to_float(r);
  printf("0.5 * 0.5 = %f (esperado 0.25)\n", fr);

  a  = float_to_takum(1e5);
  b  = float_to_takum(1e5);
  r  = takum_multiply(a, b);
  fr = takum_to_float(r);
  printf("1e5 * 1e5 = %f (esperado 1e10)\n", fr);

  return 0;
}
