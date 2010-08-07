	#define S1(zT0,zT1,zT2,t,i,j)	a[i][j]=(a[i-1][j-1]+a[i-1][j]+a[i-1][j+1]+a[i][j-1]+a[i][j]+a[i][j+1]+a[i+1][j-1]+a[i+1][j]+a[i+1][j+1])/9.0;

		int t1, t2, t3, t4, t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.08s. */
if ((N >= 3) && (T >= 1)) {
  for (t1=0;t1<=floord(T-1,10);t1++) {
    for (t2=ceild(2*t1-2,3);t2<=min(floord(T+N-3,15),floord(10*t1+N+7,15));t2++) {
      for (t3=max(max(ceild(4*t1-2,3),ceild(2*t1+3*t2-2,3)),ceild(30*t2-N-11,15));t3<=min(min(min(min(floord(2*T+2*N-6,15),floord(20*t1+2*N+14,15)),floord(30*t2+N+25,15)),floord(10*t1+15*t2+N+21,15)),floord(15*t2+T+N+11,15));t3++) {
        for (t4=max(max(max(ceild(15*t3-2*N+4,2),10*t1),15*t2-N+2),-15*t2+15*t3-N-12);t4<=min(min(min(min(floord(15*t3+12,2),T-1),10*t1+9),15*t2+13),-15*t2+15*t3+13);t4++) {
          for (t5=max(max(15*t2,t4+1),15*t3-t4-N+2);t5<=min(min(15*t2+14,15*t3-t4+13),t4+N-2);t5++) {
            for (t6=max(15*t3,t4+t5+1);t6<=min(15*t3+14,t4+t5+N-2);t6++) {
              S1(t1,t2,t3,t4,-t4+t5,-t4-t5+t6);
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
