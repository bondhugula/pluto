	#define S1(zT0,zT1,zT2,i,j,k)	b[j][k]+=a[i][k]*b[j][i];

		int t1, t2, t3, t4, t4t, newlb_t4, newub_t4, t5, t5t, newlb_t5, newub_t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.02s. */
if (N >= 2) {
  for (t1=0;t1<=floord(N-1,32);t1++) {
    for (t2=0;t2<=floord(N-1,32);t2++) {
      for (t3=0;t3<=min(floord(N-2,32),t2);t3++) {
/*@ begin Loop(
	transform RegTile(loops=['t4','t5'], ufactors=[8,8])
        for (t4=32*t1;t4<=min(N-1,32*t1+31);t4++) 
          for (t5=max(32*t2,32*t3+1);t5<=min(N-1,32*t2+31);t5++) 
{
{
	lbv=32*t3; 	ubv=min(32*t3+31,t5-1);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
              S1(t1,t2,t3,t5,t4,t6);
            }
}
}
) @*/{
  for (t4t=32*t1; t4t<=min(N-1,32*t1+31)-7; t4t=t4t+8) {
    for (t5t=max(32*t2,32*t3+1); t5t<=min(N-1,32*t2+31)-7; t5t=t5t+8) {
      for (t5=t5t; t5<=t5t+7; t5=t5+1) 
{
	lbv=32*t3; 	ubv=min(32*t3+31,t5-1);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t5,t4t,t6);
          S1(t1,t2,t3,t5,(t4t+1),t6);
          S1(t1,t2,t3,t5,(t4t+2),t6);
          S1(t1,t2,t3,t5,(t4t+3),t6);
          S1(t1,t2,t3,t5,(t4t+4),t6);
          S1(t1,t2,t3,t5,(t4t+5),t6);
          S1(t1,t2,t3,t5,(t4t+6),t6);
          S1(t1,t2,t3,t5,(t4t+7),t6);
        }
}
    }
    for (t5=t5t; t5<=min(N-1,32*t2+31); t5=t5+1) {
{
	lbv=32*t3; 	ubv=min(32*t3+31,t5-1);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
        S1(t1,t2,t3,t5,t4t,t6);
        S1(t1,t2,t3,t5,(t4t+1),t6);
        S1(t1,t2,t3,t5,(t4t+2),t6);
        S1(t1,t2,t3,t5,(t4t+3),t6);
        S1(t1,t2,t3,t5,(t4t+4),t6);
        S1(t1,t2,t3,t5,(t4t+5),t6);
        S1(t1,t2,t3,t5,(t4t+6),t6);
        S1(t1,t2,t3,t5,(t4t+7),t6);
      }
}
    }
  }
  for (t4=t4t; t4<=min(N-1,32*t1+31); t4=t4+1) {
    for (t5t=max(32*t2,32*t3+1); t5t<=min(N-1,32*t2+31)-7; t5t=t5t+8) {
      for (t5=t5t; t5<=t5t+7; t5=t5+1) 
{
	lbv=32*t3; 	ubv=min(32*t3+31,t5-1);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t5,t4,t6);
        }
}
    }
    for (t5=t5t; t5<=min(N-1,32*t2+31); t5=t5+1) {
{
	lbv=32*t3; 	ubv=min(32*t3+31,t5-1);
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
