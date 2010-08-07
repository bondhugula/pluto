
	#define S1(zT0,zT1,zT2,zT3,zT4,zT5,k,j)	a[k][j]=a[k][j]/a[k][k];
	#define S2(zT0,zT1,zT2,zT3,zT4,zT5,k,i,j)	a[i][j]=a[i][j]-a[i][k]*a[k][j];

		int t1, t2, t3, t4, t5, t6, t7, t8, t9;

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.22s. */
if (N >= 2) {
  for (t1=0;t1<=floord(57*N-82,6400);t1++) {
	lb1=max(ceild(32*t1-24,57),ceild(256*t1-N+2,256));
	ub1=min(floord(N-1,200),t1);
#pragma omp parallel for shared(t1,lb1,ub1) private(t2,t3,t4,t5,t6,t7,t8,t9)
	for (t2=lb1; t2<=ub1; t2++) {
      for (t3=t1-t2;t3<=floord(N-1,256);t3++) {
        for (t4=16*t1-16*t2;t4<=min(min(floord(25*t2+24,2),floord(N-2,16)),16*t1-16*t2+15);t4++) {
          for (t5=max(ceild(4*t4-24,25),2*t2);t5<=min(floord(N-1,100),2*t2+1);t5++) {
            for (t6=max(16*t3,t4);t6<=min(floord(N-1,16),16*t3+15);t6++) {
              if ((t1 == t2+t3) && (t4 == t6)) {
                for (t7=16*t4;t7<=min(min(N-2,16*t4+14),100*t5+98);t7++) {
                  for (t8=max(100*t5,t7+1);t8<=min(N-1,100*t5+99);t8++) {
                    if (t1 >= 0) {
                      S1(t1-t2,t2,t1-t2,t4,t5,t4,t7,t8);
                    }
                    for (t9=t7+1;t9<=min(N-1,16*t4+15);t9++) {
                      S2(t1-t2,t2,t1-t2,t4,t5,t4,t7,t9,t8);
                    }
                  }
                }
              }
              if ((t1 == t2+t3) && (t4 == t6) && (t4 <= min(floord(25*t5+20,4),floord(N-17,16)))) {
                for (t8=max(100*t5,16*t4+16);t8<=min(N-1,100*t5+99);t8++) {
                  S1(t1-t2,t2,t1-t2,t4,t5,t4,16*t4+15,t8);
                }
              }
              if (t4 <= t6-1) {
                for (t7=16*t4;t7<=min(16*t4+15,100*t5+98);t7++) {
                  for (t8=max(100*t5,t7+1);t8<=min(N-1,100*t5+99);t8++) {
                    for (t9=16*t6;t9<=min(N-1,16*t6+15);t9++) {
                      S2(t1-t2,t2,t3,t4,t5,t6,t7,t9,t8);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
