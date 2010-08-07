	#define S1(zT0,zT1,zT2,zT3,zT4,zT5,i,j,k)	C[i][j]=beta*C[i][j]+alpha*A[i][k]*B[k][j];

		int t1, t2, t3, t4, t5, t6, t7, t7t, newlb_t7, newub_t7, t8, t8t, newlb_t8, newub_t8, t9;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.05s. */
if ((K >= 1) && (M >= 1) && (N >= 1)) {
  for (t1=0;t1<=floord(M-1,128);t1++) {
    for (t2=0;t2<=floord(N-1,256);t2++) {
      for (t3=0;t3<=floord(K-1,128);t3++) {
        for (t4=16*t1;t4<=min(floord(M-1,8),16*t1+15);t4++) {
          for (t5=2*t2;t5<=min(floord(N-1,128),2*t2+1);t5++) {
            for (t6=16*t3;t6<=min(floord(K-1,8),16*t3+15);t6++) {
/*@ begin Loop(
	transform RegTile(loops=['t7','t8'], ufactors=[8,8])
              for (t7=8*t4;t7<=min(M-1,8*t4+7);t7++) 
                for (t8=8*t6;t8<=min(K-1,8*t6+7);t8++) 
{
{
	lbv=128*t5; 	ubv=min(N-1,128*t5+127);
#pragma ivdep
#pragma vector always
	for (t9=lbv; t9<=ubv; t9++) {
                    S1(t1,t2,t3,t4,t5,t6,t7,t9,t8);
                  }
}
}
) @*/{
  for (t7t=8*t4; t7t<=min(M-1,8*t4+7)-7; t7t=t7t+8) {
    for (t8t=8*t6; t8t<=min(K-1,8*t6+7)-7; t8t=t8t+8) {
{
	lbv=128*t5; 	ubv=min(N-1,128*t5+127);
#pragma ivdep
#pragma vector always
	for (t9=lbv; t9<=ubv; t9++) {
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,(t8t+7));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,(t8t+7));
      }
}
    }
    for (t8=t8t; t8<=min(K-1,8*t6+7); t8=t8+1) {
{
	lbv=128*t5; 	ubv=min(N-1,128*t5+127);
#pragma ivdep
#pragma vector always
	for (t9=lbv; t9<=ubv; t9++) {
        S1(t1,t2,t3,t4,t5,t6,t7t,t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+1),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+2),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+3),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+4),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+5),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+6),t9,t8);
        S1(t1,t2,t3,t4,t5,t6,(t7t+7),t9,t8);
      }
}
    }
  }
  for (t7=t7t; t7<=min(M-1,8*t4+7); t7=t7+1) {
    for (t8t=8*t6; t8t<=min(K-1,8*t6+7)-7; t8t=t8t+8) {
{
	lbv=128*t5; 	ubv=min(N-1,128*t5+127);
#pragma ivdep
#pragma vector always
	for (t9=lbv; t9<=ubv; t9++) {
        S1(t1,t2,t3,t4,t5,t6,t7,t9,t8t);
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+1));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+2));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+3));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+4));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+5));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+6));
        S1(t1,t2,t3,t4,t5,t6,t7,t9,(t8t+7));
      }
}
    }
    for (t8=t8t; t8<=min(K-1,8*t6+7); t8=t8+1) {
{
	lbv=128*t5; 	ubv=min(N-1,128*t5+127);
#pragma ivdep
#pragma vector always
	for (t9=lbv; t9<=ubv; t9++) {
        S1(t1,t2,t3,t4,t5,t6,t7,t9,t8);
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
    }
  }
}
/* End of CLooG code */
