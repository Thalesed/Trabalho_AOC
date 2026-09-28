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

double takum_to_double(takum16_t num){
  uint16_t tmp = num >> 15;
  uint16_t S = tmp; // S -> 1 bit
  tmp = num >> 14 ;
  tmp -= 2*S;
  uint16_t D = tmp; // D -> 1 bit
  tmp = num >> 11;
  tmp -= (S*16 + D*8);
  uint16_t R = tmp; // R -> 3 bits

  tmp = num >> (11 - R);
  tmp -= S*(1<<(R+4));
  tmp -= D*(1<<(R+3));
  uint16_t tmp2 = R >> 2; // bit 3 de R
  tmp -= tmp2*(1<<(R+2));
  uint16_t tmp3 = (R >> 1) - tmp2*2; // bit 2 de R
  tmp -= tmp3*(1<<R+1);
  uint16_t tmp4 = R - tmp2*4 - tmp3*2; // bit 1 de R
  tmp -= tmp4*(1<<R);
  uint16_t C = tmp; // C -> R bits
  tmp = num - S*(1<<15);
  tmp -= D*(1<<14);
  tmp -= (tmp2*(1<<13) + tmp3*(1<<12) + tmp4*(1<<11) );
  tmp2 = C;
  for(uint16_t i=0;i<R;i++){
    tmp3 = tmp2 >> (R-1-i);
    tmp -= tmp3*(1<<(10-i));
    tmp2 -= tmp3*(1<<(R-1-i));
  }
  uint16_t F = tmp;

  // printBinary(S);
  // printBinary(D);
  // printBinary(R);
  // printBinary(C);
  // printBinary(F);

  double resultado;

  int16_t c;

  resultado = (F*1.0) / (1<<(11-R));
  
  if(D == 0){
    R = 7 - R;
  }
  c = C + C_BIAS[R+D*(1<<3)]; // + c-bias
  // printf("%d\n", c);
  resultado += c;

  resultado = pow(2, resultado);
  if(S == 1){
    resultado *= (-1);
  }

  return resultado;
}

takum16_t double_to_takum(double x){
  if(x == 0.0){
    return 0x0000; 
  }
  
  uint16_t S = (x < 0) ? 1 : 0;
  double abs_x = fabs(x);
  
  double l = log2(abs_x);
  
  uint16_t D = (abs_x >= 1.0) ? 1 : 0;
  
  int16_t c = (int16_t)floor(l);
  double m = l - c;   // m ∈ [0, 1)
  
  uint16_t R_adj = 0;
  uint16_t R = 0;
  uint16_t C = 0;
  bool found = false;
  
  for(R_adj = 0; R_adj <= 7; R_adj++){
    int16_t bias = C_BIAS[R_adj + D*8];
    
    if(D == 1){
      R = R_adj;
    } else {
      R = 7 - R_adj;
    }
    
    int16_t C_candidate = c - bias;
    int16_t max_C = (1 << R) - 1;
    
    if(C_candidate >= 0 && C_candidate <= max_C){
      C = (uint16_t)C_candidate;
      found = true;
      break;
    }
  }
  
  if(!found){
    return 0xFFFF;
  }
  
  uint16_t F = (uint16_t)round(m * (1 << (11 - R)));
  
  if(F >= (1 << (11 - R))){
    F = 0;
    C++;
  }
  
  uint16_t num = 0;
  num |= (S & 1) << 15;
  num |= (D & 1) << 14;
  
  uint16_t R_field;
  if(D == 1){
    R_field = R_adj;
  } else {
    R_field = 7 - R_adj;
  }
  num |= (R_field & 7) << 11;
  
  num |= (C & ((1 << R) - 1)) << (11 - R);
  
  num |= (F & ((1 << (11 - R)) - 1));
  
  
  return num;
}

takum16_t takum_multiply(takum16_t a, takum16_t b) {
  uint16_t sign_a = (a >> 15) & 1;
  uint16_t sign_b = (b >> 15) & 1;
  uint16_t result_sign = sign_a ^ sign_b; 
  
  double fa = takum_to_double(a);
  double fb = takum_to_double(b);
  
  double l_a = 2.0 * log(fabs(fa));
  double l_b = 2.0 * log(fabs(fb));
  
  double l_result = l_a + l_b;
  
  double f_result = exp(l_result / 2.0);
  if (result_sign) {
    f_result = -f_result;
  }
  return double_to_takum(f_result);
}

int main(){
  takum16_t ta, tb;
  
  ta = double_to_takum(2.5);
  tb = double_to_takum(3.0);
  
  double a = takum_to_double(takum_multiply(ta, tb));
  
  printf("%f", a);
  
  return 0;
}