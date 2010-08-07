
/*@ begin PerfTuning (
  import spec align_unroll;
) @*/

  int i;

  /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
  /*@ begin Loop ( 
      transform Unroll(ufactor=UF) 
        for (i = 0; i <= n-1; i++)
          y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
  ) @*/
  for (i = 0; i < n; i++)
     y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
  /*@ end @*/
  /*@ end @*/
 
/*@ end @*/
