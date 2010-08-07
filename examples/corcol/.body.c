	#define S1(zT0,zT1,zT2,i,j)	data[i][j]-=mean[j];
	#define S2(zT0,zT1,zT2,i,j)	data[i][j]/=sqrt(float_n)*stddev[j];
	#define S3(zT0,zT1,zT2,j1)	symmat[j1][j1]=1.0;
	#define S4(zT0,zT1,zT2,j1,j2)	symmat[j1][j2]=0.0;
	#define S5(zT0,zT1,zT2,j1,j2,i)	symmat[j1][j2]+=(data[i][j1]*data[i][j2]);
	#define S6(zT0,zT1,zT2,j1,j2)	symmat[j2][j1]=symmat[j1][j2];

		int t1, t2, t3, t4, t5, t5t, newlb_t5, newub_t5, t6, t6t, newlb_t6, newub_t6, t7, t8, t9;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.09s. */
for (t2=0;t2<=floord(m-1,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t5=max(1,32*t2);t5<=min(min(m-1,32*t2+31),32*t3+30);t5++) {
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S4(t2,t3,0,t5,t7);
      }
}
    }
  }
}
for (t2=0;t2<=floord(m-1,32);t2++) {
  for (t5=max(1,32*t2);t5<=min(m-1,32*t2+31);t5++) {
    S3(t2,0,0,t5);
  }
}
for (t2=0;t2<=floord(n,32);t2++) {
  for (t3=0;t3<=floord(m,32);t3++) {
    for (t5=max(1,32*t2);t5<=min(n,32*t2+31);t5++) {
{
	lbv=max(1,32*t3); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S1(t2,t3,0,t5,t7);
        S2(t2,t3,0,t5,t7);
      }
}
    }
  }
}
for (t2=0;t2<=floord(m-1,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t4=0;t4<=floord(n,32);t4++) {
/*@ begin Loop(
	transform RegTile(loops=['t5','t6'], ufactors=[8,8])
      for (t5=max(1,32*t2);t5<=min(min(m-1,32*t2+31),32*t3+30);t5++) 
        for (t6=max(1,32*t4);t6<=min(n,32*t4+31);t6++) 
{
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
            S5(t2,t3,t4,t5,t7,t6);
          }
}
}
) @*/{
  for (t5t=max(1,32*t2); t5t<=min(min(m-1,32*t2+31),32*t3+30)-7; t5t=t5t+8) {
    for (t6t=max(1,32*t4); t6t<=min(n,32*t4+31)-7; t6t=t6t+8) {
      for (t5=t5t; t5<=t5t+7; t5=t5+1) 
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
          S5(t2,t3,t4,t5,t7,t6t);
          S5(t2,t3,t4,t5,t7,(t6t+1));
          S5(t2,t3,t4,t5,t7,(t6t+2));
          S5(t2,t3,t4,t5,t7,(t6t+3));
          S5(t2,t3,t4,t5,t7,(t6t+4));
          S5(t2,t3,t4,t5,t7,(t6t+5));
          S5(t2,t3,t4,t5,t7,(t6t+6));
          S5(t2,t3,t4,t5,t7,(t6t+7));
        }
}
    }
    for (t6=t6t; t6<=min(n,32*t4+31); t6=t6+1) {
      for (t5=t5t; t5<=t5t+7; t5=t5+1) 
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
          S5(t2,t3,t4,t5,t7,t6);
        }
}
    }
  }
  for (t5=t5t; t5<=min(min(m-1,32*t2+31),32*t3+30); t5=t5+1) {
    for (t6t=max(1,32*t4); t6t<=min(n,32*t4+31)-7; t6t=t6t+8) {
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S5(t2,t3,t4,t5,t7,t6t);
        S5(t2,t3,t4,t5,t7,(t6t+1));
        S5(t2,t3,t4,t5,t7,(t6t+2));
        S5(t2,t3,t4,t5,t7,(t6t+3));
        S5(t2,t3,t4,t5,t7,(t6t+4));
        S5(t2,t3,t4,t5,t7,(t6t+5));
        S5(t2,t3,t4,t5,t7,(t6t+6));
        S5(t2,t3,t4,t5,t7,(t6t+7));
      }
}
    }
    for (t6=t6t; t6<=min(n,32*t4+31); t6=t6+1) {
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S5(t2,t3,t4,t5,t7,t6);
      }
}
    }
  }
}
/*@ end @*/
    }
  }
}
for (t2=0;t2<=floord(m-1,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t5=max(1,32*t2);t5<=min(min(m-1,32*t2+31),32*t3+30);t5++) {
{
	lbv=max(32*t3,t5+1); 	ubv=min(m,32*t3+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S6(t2,t3,0,t5,t7);
      }
}
    }
  }
}
/* End of CLooG code */
