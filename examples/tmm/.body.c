	#define S1(zT0,zT1,zT2,i,j,k)	C[i][j]+=A[i][k]*B[k][j];

		int t1, t2, t3, t4, t4t, newlb_t4, newub_t4, t5, t5t, newlb_t5, newub_t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.04s. */
if ((Ni >= 1) && (Nj >= 1) && (Nk >= 1)) {
  for (t1=0;t1<=min(min(floord(Ni-1,32),floord(Nj-1,32)),floord(Nk-1,32));t1++) {
    for (t2=t1;t2<=floord(Nj-1,32);t2++) {
      for (t3=t1;t3<=floord(Nk-1,32);t3++) {
/*@ begin Loop(
	transform RegTile(loops=['t4','t5'], ufactors=[8,8])
        for (t4=32*t1;t4<=min(min(min(Ni-1,Nj-1),Nk-1),32*t1+31);t4++) 
          for (t5=max(32*t3,t4);t5<=min(Nk-1,32*t3+31);t5++) 
{
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
              S1(t1,t2,t3,t4,t6,t5);
            }
}
}
) @*/{
  for (t4t=32*t1; t4t<=min(min(min(Ni-1,Nj-1),Nk-1),32*t1+31)-7; t4t=t4t+8) {
    newlb_t5=-2147483648;
    newub_t5=min(Nk-1,32*t3+31);
    for (t4=t4t; t4<=t4t+7; t4=t4+1) 
      newlb_t5=max(newlb_t5,max(32*t3,t4));
    for (t4=t4t; t4<=t4t+7; t4=t4+1) 
      for (t5=max(32*t3,t4); t5<=newlb_t5-1; t5=t5+1) {
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t4,t6,t5);
        }
}
      }
    for (t5t=newlb_t5; t5t<=newub_t5-7; t5t=t5t+8) {
      for (t4=t4t; t4<=t4t+7; t4=t4+1) 
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t4,t6,t5t);
          S1(t1,t2,t3,t4,t6,(t5t+1));
          S1(t1,t2,t3,t4,t6,(t5t+2));
          S1(t1,t2,t3,t4,t6,(t5t+3));
          S1(t1,t2,t3,t4,t6,(t5t+4));
          S1(t1,t2,t3,t4,t6,(t5t+5));
          S1(t1,t2,t3,t4,t6,(t5t+6));
          S1(t1,t2,t3,t4,t6,(t5t+7));
        }
}
    }
    for (t5=t5t; t5<=newub_t5; t5=t5+1) {
      for (t4=t4t; t4<=t4t+7; t4=t4+1) 
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t4,t6,t5);
        }
}
    }
    for (t4=t4t; t4<=t4t+7; t4=t4+1) 
      for (t5=newub_t5+1; t5<=min(Nk-1,32*t3+31); t5=t5+1) {
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
          S1(t1,t2,t3,t4,t6,t5);
        }
}
      }
  }
  for (t4=t4t; t4<=min(min(min(Ni-1,Nj-1),Nk-1),32*t1+31); t4=t4+1) {
    for (t5t=max(32*t3,t4); t5t<=min(Nk-1,32*t3+31)-7; t5t=t5t+8) {
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
        S1(t1,t2,t3,t4,t6,t5t);
        S1(t1,t2,t3,t4,t6,(t5t+1));
        S1(t1,t2,t3,t4,t6,(t5t+2));
        S1(t1,t2,t3,t4,t6,(t5t+3));
        S1(t1,t2,t3,t4,t6,(t5t+4));
        S1(t1,t2,t3,t4,t6,(t5t+5));
        S1(t1,t2,t3,t4,t6,(t5t+6));
        S1(t1,t2,t3,t4,t6,(t5t+7));
      }
}
    }
    for (t5=t5t; t5<=min(Nk-1,32*t3+31); t5=t5+1) {
{
	lbv=max(32*t2,t4); 	ubv=min(Nj-1,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t6=lbv; t6<=ubv; t6++) {
        S1(t1,t2,t3,t4,t6,t5);
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
