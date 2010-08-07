
	#define S1(zT0,zT1,zT2,i,j,k)	c[j][k]+=a[i][j]*a[i][k];

		int t1, t2, t3, t4, t4t, newlb_t4, newub_t4, t5, t5t, newlb_t5, newub_t5, t6;

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.02s. */
if (N >= 1) {
	lb1=0;
	ub1=floord(N-1,32);
#pragma omp parallel for shared(lb1,ub1) private(t1,t2,t3,t4,t5,t6)
	for (t1=lb1; t1<=ub1; t1++) {
    for (t2=t1;t2<=floord(N-1,32);t2++) {
      for (t3=0;t3<=floord(N-1,32);t3++) {
/*@ begin Loop(
	transform RegTile(loops=['t4','t5'], ufactors=[8,8])
        for (t4=32*t1;t4<=min(N-1,32*t1+31);t4++) 
          for (t5=32*t3;t5<=min(N-1,32*t3+31);t5++) 
{
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
              S1(t1,t2,t3,t5,t4,t6);
            }
}
}
) @*/{
  for (t4t=32*t1; t4t<=min(N-1,32*t1+31)-7; t4t=t4t+8) {
    for (t5t=32*t3; t5t<=min(N-1,32*t3+31)-7; t5t=t5t+8) {
      for (t4=t4t; t4<=t4t+7; t4=t4+1) 
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t5t,t4,t6);
          S1(t1,t2,t3,(t5t+1),t4,t6);
          S1(t1,t2,t3,(t5t+2),t4,t6);
          S1(t1,t2,t3,(t5t+3),t4,t6);
          S1(t1,t2,t3,(t5t+4),t4,t6);
          S1(t1,t2,t3,(t5t+5),t4,t6);
          S1(t1,t2,t3,(t5t+6),t4,t6);
          S1(t1,t2,t3,(t5t+7),t4,t6);
        }
}
    }
    for (t5=t5t; t5<=min(N-1,32*t3+31); t5=t5+1) {
      for (t4=t4t; t4<=t4t+7; t4=t4+1) 
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t5,t4,t6);
        }
}
    }
  }
  for (t4=t4t; t4<=min(N-1,32*t1+31); t4=t4+1) {
    for (t5t=32*t3; t5t<=min(N-1,32*t3+31)-7; t5t=t5t+8) {
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
        S1(t1,t2,t3,t5t,t4,t6);
        S1(t1,t2,t3,(t5t+1),t4,t6);
        S1(t1,t2,t3,(t5t+2),t4,t6);
        S1(t1,t2,t3,(t5t+3),t4,t6);
        S1(t1,t2,t3,(t5t+4),t4,t6);
        S1(t1,t2,t3,(t5t+5),t4,t6);
        S1(t1,t2,t3,(t5t+6),t4,t6);
        S1(t1,t2,t3,(t5t+7),t4,t6);
      }
}
    }
    for (t5=t5t; t5<=min(N-1,32*t3+31); t5=t5+1) {
{
	lbv=max(32*t2,t4); 	ubv=min(N-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
        S1(t1,t2,t3,t5,t4,t6);
      }
}
    }
  }
}
/*@ end @*/
      }
    }
  }
}
/* End of CLooG code */
