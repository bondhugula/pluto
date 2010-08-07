
/*@ begin PerfTuning (
  import spec align_unroll;
) @*/

  int i;

  /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
  #pragma disjoint (*x1,*x2,*x3,*x4,*y) 
  if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
    __alignx(16,x1); 
    __alignx(16,x2); 
    __alignx(16,x3); 
    __alignx(16,x4); 
    __alignx(16,y); 

    /*@ begin Loop ( 
        transform Unroll(ufactor=UF) 
          for (i = 0; i <= n-1; i++)
            y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
    ) @*/  {
      for (i=0; i<=n-8; i=i+8) {
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
      }
      for (; i<=n-1; i=i+1) 
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
    }
  /*@ end @*/
    
  } else {

    /*@ begin Loop ( 
        transform Unroll(ufactor=UF) 
          for (i = 0; i <= n-1; i++)
            y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
    ) @*/  {
      for (i=0; i<=n-8; i=i+8) {
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
      }
      for (; i<=n-1; i=i+1) 
        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
    }
  /*@ end @*/
    
  } 
  /*@ end @*/
 
/*@ end @*/
