

#include <stdio.h>
#include <sys/time.h>
#include <math.h>



#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define T TVAL
#define N NVAL
#define coeff1 0.5
#define coeff2 0.7
double h[N];
double e[N+1];

void init_arrays()
{
  int i1;
  for (i1=0; i1<N; i1++)
   h[i1] = (i1) % 5 + 1;
  for (i1=0; i1<N+1; i1++)
   e[i1] = (i1) % 5 + 1;
}

double rtclock()
{
  struct timezone tzp;
  struct timeval tp;
  int stat;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main()
{
  init_arrays();

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  for (annot_i=0; annot_i<REPS; annot_i++)
  {
    annot_t_start = rtclock();
    
  




int t, i, ii, j, k, l;

	#define S1(zT0,zT1,t,i)	{e[i]=e[i]-coeff1*(h[i]-h[i-1]);}
	#define S2(zT0,zT1,t,i)	{h[i]=h[i]-coeff2*(e[1+i]-e[i]);}

	int c1, c2, c3, c4, c5;

	register int lbv, ubv;

for (c1=0;c1<=floord(T,512);c1++) {
  for (c2=max(ceild(256*c1-255,256),0);c2<=min(floord(512*c1+N+511,512),floord(N+T,512));c2++) {
 if ((c1 <= floord(512*c2-N,512)) && (c2 >= ceild(N+1,512))) {
      S2(c1,-c1+c2,512*c2-N,N-1) ;
    }
    for (c3=max(max(512*c2-N+1,512*c1),1);c3<=min(min(T,512*c1+511),512*c2-N+511);c3++) {
      for (c4=max(c3+1,512*c2);c4<=c3+N-1;c4++) {
        S1(c1,-c1+c2,c3,-c3+c4) ;
        S2(c1,-c1+c2,c3,-c3+c4-1) ;
      }
      S2(c1,-c1+c2,c3,N-1) ;
    }
/*@ begin Loop(
  transform Composite(
    tile = [('c3',T1,'ii')],
    unrolljam = [('c3',U1),('c4',U2)],
    vector = (VEC, ['ivdep','vector always'])
    )                             
    for (c3=max(max(1,512*c1),512*c2-N+512);c3<=min(min(512*c1+511,T),512*c2+510);c3++) 
      for (c4=max(c3+1,512*c2);c4<=512*c2+511;c4++) 
{
        S1(c1,-c1+c2,c3,-c3+c4) ;
        S2(c1,-c1+c2,c3,-c3+c4-1) ;
}
) @*/{
  for (c3=max(max(1,512*c1),512*c2-N+512); c3<=min(min(512*c1+511,T),512*c2+510)-7; c3=c3+8) {
    register int cbv_1, cbv_2;
    cbv_1=max(c3+1,512*c2);
    cbv_2=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_1; c4<=cbv_2; c4=c4+1) {
      S1(c1,-c1+c2,c3,-c3+c4);
      S2(c1,-c1+c2,c3,-c3+c4-1);
    }
    register int cbv_3, cbv_4;
    cbv_3=max(c3+2,512*c2);
    cbv_4=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_3; c4<=cbv_4; c4=c4+1) {
      S1(c1,-c1+c2,(c3+1),-(c3+1)+c4);
      S2(c1,-c1+c2,(c3+1),-(c3+1)+c4-1);
    }
    register int cbv_5, cbv_6;
    cbv_5=max(c3+3,512*c2);
    cbv_6=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_5; c4<=cbv_6; c4=c4+1) {
      S1(c1,-c1+c2,(c3+2),-(c3+2)+c4);
      S2(c1,-c1+c2,(c3+2),-(c3+2)+c4-1);
    }
    register int cbv_7, cbv_8;
    cbv_7=max(c3+4,512*c2);
    cbv_8=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_7; c4<=cbv_8; c4=c4+1) {
      S1(c1,-c1+c2,(c3+3),-(c3+3)+c4);
      S2(c1,-c1+c2,(c3+3),-(c3+3)+c4-1);
    }
    register int cbv_9, cbv_10;
    cbv_9=max(c3+5,512*c2);
    cbv_10=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_9; c4<=cbv_10; c4=c4+1) {
      S1(c1,-c1+c2,(c3+4),-(c3+4)+c4);
      S2(c1,-c1+c2,(c3+4),-(c3+4)+c4-1);
    }
    register int cbv_11, cbv_12;
    cbv_11=max(c3+6,512*c2);
    cbv_12=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_11; c4<=cbv_12; c4=c4+1) {
      S1(c1,-c1+c2,(c3+5),-(c3+5)+c4);
      S2(c1,-c1+c2,(c3+5),-(c3+5)+c4-1);
    }
    register int cbv_13, cbv_14;
    cbv_13=max(c3+7,512*c2);
    cbv_14=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_13; c4<=cbv_14; c4=c4+1) {
      S1(c1,-c1+c2,(c3+6),-(c3+6)+c4);
      S2(c1,-c1+c2,(c3+6),-(c3+6)+c4-1);
    }
    register int cbv_15, cbv_16;
    cbv_15=max(c3+8,512*c2);
    cbv_16=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_15; c4<=cbv_16; c4=c4+1) {
      S1(c1,-c1+c2,(c3+7),-(c3+7)+c4);
      S2(c1,-c1+c2,(c3+7),-(c3+7)+c4-1);
    }
  }
  for (; c3<=min(min(512*c1+511,T),512*c2+510); c3=c3+1) {
    register int cbv_17, cbv_18;
    cbv_17=max(c3+1,512*c2);
    cbv_18=512*c2+511;
#pragma ivdep
#pragma vector always
    for (c4=cbv_17; c4<=cbv_18; c4=c4+1) {
      S1(c1,-c1+c2,c3,-c3+c4);
      S2(c1,-c1+c2,c3,-c3+c4-1);
    }
  }
}
/*@ end @*/
  }
}


    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return 1;
}
                                                                   
