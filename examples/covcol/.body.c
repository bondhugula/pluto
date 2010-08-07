	#define S1(zT0,zT1,zT2,j)	mean[j]=0.0;
	#define S2(zT0,zT1,zT2,j,i)	mean[j]+=data[i][j];
	#define S3(zT0,zT1,zT2,j)	mean[j]/=float_n;
	#define S4(zT0,zT1,zT2,i,j)	data[i][j]-=mean[j];
	#define S5(zT0,zT1,zT2,j1,j2)	symmat[j1][j2]=0.0;
	#define S6(zT0,zT1,zT2,j1,j2,i)	symmat[j1][j2]+=data[i][j1]*data[i][j2];
	#define S7(zT0,zT1,zT2,j1,j2)	symmat[j2][j1]=symmat[j1][j2];

		int t1, t2, t3, t4, t5, t6, t7;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 0.12s. */
for (t2=0;t2<=floord(m,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t6=max(1,32*t3);t6<=min(m,32*t3+31);t6++) {
{
	lbv=max(1,32*t2); 	ubv=min(t6,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S5(t2,t3,0,t7,t6);
      }
}
    }
  }
}
for (t2=0;t2<=floord(m,32);t2++) {
{
	lbv=max(1,32*t2); 	ubv=min(m,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
    S1(t2,0,0,t7);
  }
}
}
for (t2=0;t2<=floord(m,32);t2++) {
  for (t3=0;t3<=floord(n,32);t3++) {
    for (t6=max(1,32*t3);t6<=min(n,32*t3+31);t6++) {
{
	lbv=max(1,32*t2); 	ubv=min(m,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S2(t2,t3,0,t7,t6);
      }
}
    }
  }
}
for (t2=0;t2<=floord(m,32);t2++) {
{
	lbv=max(1,32*t2); 	ubv=min(m,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
    S3(t2,0,0,t7);
  }
}
}
for (t2=0;t2<=floord(n,32);t2++) {
  for (t3=0;t3<=floord(m,32);t3++) {
    for (t6=max(1,32*t3);t6<=min(m,32*t3+31);t6++) {
{
	lbv=max(1,32*t2); 	ubv=min(n,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S4(t2,t3,0,t7,t6);
      }
}
    }
  }
}
for (t2=0;t2<=floord(m,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t4=0;t4<=floord(n,32);t4++) {
      for (t5=max(1,32*t4);t5<=min(n,32*t4+31);t5++) {
        for (t6=max(1,32*t3);t6<=min(m,32*t3+31);t6++) {
{
	lbv=max(1,32*t2); 	ubv=min(t6,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
            S6(t2,t3,t4,t7,t6,t5);
          }
}
        }
      }
    }
  }
}
for (t2=0;t2<=floord(m,32);t2++) {
  for (t3=t2;t3<=floord(m,32);t3++) {
    for (t6=max(1,32*t3);t6<=min(m,32*t3+31);t6++) {
{
	lbv=max(1,32*t2); 	ubv=min(t6,32*t2+31);
#pragma ivdep
#pragma vector always
	for (t7=lbv; t7<=ubv; t7++) {
        S7(t2,t3,0,t7,t6);
      }
}
    }
  }
}
/* End of CLooG code */
