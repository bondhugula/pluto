/* 4-cascaded IIR biquad filter processing 64 points */
/* Modified to use arrays - SMP */

//#include "traps.h"

#define NPOINTS 64
#define NSECTIONS 4

float input[NPOINTS];
float output[NPOINTS];
float coefficient[NSECTIONS][NSECTIONS][NSECTIONS];
float internal_state[8][NSECTIONS][2];

void iir(float input[NPOINTS], float output[NPOINTS],
         float coefficient[NSECTIONS][NSECTIONS][NSECTIONS],
         float internal_state[8][NSECTIONS][2]);

main()
{
  int i;

  float *pcoef = coefficient[0][0];
  float *pint = internal_state[0][0];

  for(i=0;i<NPOINTS;i++)
  {
    pcoef[i]=1;
    pint[i]=1;
  }
  //input_dsp(input,NPOINTS,0);

  iir(input, output, coefficient, internal_state);

  //output_dsp(input,NPOINTS,0);
  //output_dsp(coefficient,NPOINTS,0);
  //output_dsp(internal_state,NPOINTS,0);
  //output_dsp(output,NPOINTS,0);
}

void iir(float input[NPOINTS], float output[NPOINTS],
         float coefficient[NSECTIONS][NSECTIONS][NSECTIONS],
         float internal_state[8][NSECTIONS][2])
/* input:           input sample array */
/* output:          output sample array */
/* coefficient:     coefficient array */
/* internal_state:  internal state array */
{
  int i, imod8, imodNSECTIONS;
  int j;

  float state_2, state_1;
  float coef_a21, coef_a11, coef_b21, coef_b11;
  float sum;

  for (i = 0; i < NPOINTS; ++i) {

    imod8 = i % 8;
    imodNSECTIONS = i % NSECTIONS;

#pragma scop
    sum = input[i];

    for (j = 0; j < NSECTIONS; ++j) {

      state_2 = internal_state[imod8][j][0];
      state_1 = internal_state[imod8][j][1];

      sum -= internal_state[imod8][j][0] * coefficient[imodNSECTIONS][j][0] +
		internal_state[imod8][j][1] * coefficient[imodNSECTIONS][j][1];

      internal_state[imod8][j][0] = internal_state[imod8][j][1];
      internal_state[imod8][j][1] = sum;

      sum += state_2 * coefficient[imodNSECTIONS][j][2] +
				state_1 * coefficient[imodNSECTIONS][j][3];

    }

    output[i] = sum;
#pragma endscop

  }
}
