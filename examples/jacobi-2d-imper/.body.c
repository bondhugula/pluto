	#define S1(zT0,zT1,zT2,t,i,j)	b[i][j]=0.2*(a[i][j]+a[i][j-1]+a[i][1+j]+a[1+i][j]+a[i-1][j]);
	#define S2(zT0,zT1,zT2,t,i,j)	a[i][j]=b[i][j];

		int t1, t2, t3, t4, t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.47s. */
for (t1=0;t1<=floord(T-1,32);t1++) {
  for (t2=2*t1;t2<=min(floord(2*T+N-3,32),floord(64*t1+N+61,32));t2++) {
    for (t3=max(ceild(32*t2-N-27,32),2*t1);t3<=min(min(floord(2*T+N-3,32),floord(64*t1+N+61,32)),floord(32*t2+N+27,32));t3++) {
      if ((t1 <= floord(32*t3-N+1,64)) && (t2 <= t3-1)) {
        if ((N+1)%2 == 0) {
          for (t5=max(32*t2,32*t3-N+4);t5<=32*t2+31;t5++) {
            S2(t1,t2,t3,(32*t3-N+1)/2,-32*t3+t5+N-2,N-2);
          }
        }
      }
      if ((t1 <= floord(32*t2-N+1,64)) && (t2 >= t3)) {
        if ((N+1)%2 == 0) {
          for (t6=max(32*t3,32*t2-N+4);t6<=min(32*t2,32*t3+31);t6++) {
            S2(t1,t2,t3,(32*t2-N+1)/2,N-2,-32*t2+t6+N-2);
          }
        }
      }
      for (t4=max(max(ceild(32*t2-N+2,2),ceild(32*t3-N+2,2)),32*t1);t4<=min(min(min(floord(32*t2-N+32,2),floord(32*t3-N+32,2)),T-1),32*t1+31);t4++) {
        for (t5=32*t2;t5<=2*t4+N-2;t5++) {
          for (t6=32*t3;t6<=2*t4+N-2;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
          S2(t1,t2,t3,t4,-2*t4+t5-1,N-2);
        }
        for (t6=32*t3;t6<=2*t4+N-1;t6++) {
          S2(t1,t2,t3,t4,N-2,-2*t4+t6-1);
        }
      }
      for (t4=max(max(ceild(32*t2-N+33,2),ceild(32*t3-N+2,2)),32*t1);t4<=min(min(min(floord(32*t3-N+32,2),T-1),32*t1+31),16*t2-2);t4++) {
        for (t5=32*t2;t5<=32*t2+31;t5++) {
          for (t6=32*t3;t6<=2*t4+N-2;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
          S2(t1,t2,t3,t4,-2*t4+t5-1,N-2);
        }
      }
      for (t4=max(max(ceild(32*t2-N+2,2),ceild(32*t3-N+33,2)),32*t1);t4<=min(min(min(floord(32*t2-N+32,2),T-1),32*t1+31),16*t3-2);t4++) {
        for (t5=32*t2;t5<=2*t4+N-2;t5++) {
          for (t6=32*t3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
        for (t6=32*t3;t6<=32*t3+31;t6++) {
          S2(t1,t2,t3,t4,N-2,-2*t4+t6-1);
        }
      }
      for (t4=max(max(ceild(32*t2-N+33,2),ceild(32*t3-N+33,2)),32*t1);t4<=min(min(min(T-1,32*t1+31),16*t2-2),16*t3-2);t4++) {
        for (t5=32*t2;t5<=32*t2+31;t5++) {
          for (t6=32*t3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
      }
      for (t4=max(max(ceild(32*t3-N+2,2),32*t1),16*t2-1);t4<=min(min(min(floord(32*t3-N+32,2),T-1),32*t1+31),16*t2+14);t4++) {
        for (t6=32*t3;t6<=2*t4+N-2;t6++) {
          S1(t1,t2,t3,t4,2,-2*t4+t6);
        }
        for (t5=2*t4+3;t5<=32*t2+31;t5++) {
          for (t6=32*t3;t6<=2*t4+N-2;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
          S2(t1,t2,t3,t4,-2*t4+t5-1,N-2);
        }
      }
      for (t4=max(max(ceild(32*t3-N+33,2),32*t1),16*t2-1);t4<=min(min(min(T-1,32*t1+31),16*t2+14),16*t3-2);t4++) {
        for (t6=32*t3;t6<=32*t3+31;t6++) {
          S1(t1,t2,t3,t4,2,-2*t4+t6);
        }
        for (t5=2*t4+3;t5<=32*t2+31;t5++) {
          for (t6=32*t3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
      }
      for (t4=max(max(ceild(32*t2-N+2,2),32*t1),16*t3-1);t4<=min(min(min(floord(32*t2-N+32,2),T-1),32*t1+31),16*t3+14);t4++) {
        for (t5=32*t2;t5<=2*t4+N-2;t5++) {
          S1(t1,t2,t3,t4,-2*t4+t5,2);
          for (t6=2*t4+3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
        for (t6=2*t4+3;t6<=32*t3+31;t6++) {
          S2(t1,t2,t3,t4,N-2,-2*t4+t6-1);
        }
      }
      for (t4=max(max(ceild(32*t2-N+33,2),32*t1),16*t3-1);t4<=min(min(min(T-1,32*t1+31),16*t2-2),16*t3+14);t4++) {
        for (t5=32*t2;t5<=32*t2+31;t5++) {
          S1(t1,t2,t3,t4,-2*t4+t5,2);
          for (t6=2*t4+3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
      }
      for (t4=max(max(32*t1,16*t2-1),16*t3-1);t4<=min(min(min(T-1,32*t1+31),16*t2+14),16*t3+14);t4++) {
        for (t6=2*t4+2;t6<=32*t3+31;t6++) {
          S1(t1,t2,t3,t4,2,-2*t4+t6);
        }
        for (t5=2*t4+3;t5<=32*t2+31;t5++) {
          S1(t1,t2,t3,t4,-2*t4+t5,2);
          for (t6=2*t4+3;t6<=32*t3+31;t6++) {
            S1(t1,t2,t3,t4,-2*t4+t5,-2*t4+t6);
            S2(t1,t2,t3,t4,-2*t4+t5-1,-2*t4+t6-1);
          }
        }
      }
    }
  }
}
/* End of CLooG code */
