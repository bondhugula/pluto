

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define tmax T
#define nx N
#define ny N
double ex[nx][ny +1];
double ey[nx +1][ny];
double hz[nx][ny];

void init_arrays()
{
    int i, j;
    for (i=0; i<nx+1; i++)  {
        for (j=0; j<ny; j++)  {
            ey[i][j] = 0;
        }
    }
    for (i=0; i<nx; i++)  {
        for (j=0; j<ny+1; j++)  {
            ex[i][j] = 0;
        }
    }
    for (j=0; j<ny; j++)  {
        ey[0][j] = ((double)j)/ny;
    }
    for (i=0; i<nx; i++)    {
        for (j=0; j<ny; j++)  {
            hz[i][j] = 0;
        }
    }
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
    

                                    

register int i,j,k,t;  
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;  
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,  
  newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;  
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,  
  newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;  


#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

				
	int c1, c2, c3, c4, c5, c6, c7, c8, c9;

	register int lb, ub, lb1, ub1, lb2, ub2;
/* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 1.73s. */
for (c1=-1;c1<=floord(3*tmax+ny-3,64);c1++) {
	lb1=max(max(ceild(32*c1-63,96),ceild(32*c1-tmax+1,32)),0);
	ub1=min(min(floord(32*c1+31,32),floord(32*c1+ny+31,96)),floord(tmax+ny-1,64));
#pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6,c7,c8,c9)
	for (c2=lb1; c2<=ub1; c2++) {
    for (c3=max(max(max(max(max(max(ceild(32*c1-32*c2-7,8),0),ceild(-16*c1-80*c2-119,12)),ceild(4256*c1-8416*c2-3*tmax-4482,520)),ceild(64*c2-ny-6,8)),ceild(8*c1-456*c2-483,56)),ceild(424*c1-872*c2-483,56));c3<=min(min(floord(32*c1-32*c2+nx+31,8),floord(tmax+nx-1,8)),floord(64*c2+nx+62,8));c3++) {
      if ((c1 <= floord(32*c2+8*c3-nx,32)) && (c2 <= floord(8*c3-nx+ny,64)) && (c3 >= ceild(nx,8))) {
        for (c8=max(64*c2,8*c3-nx+1);c8<=min(8*c3-nx+ny,64*c2+63);c8++) {
          hz[nx-1][-8*c3+c8+nx-1]=hz[nx-1][-8*c3+c8+nx-1]-((double)(7))/10*(ey[1+nx-1][-8*c3+c8+nx-1]+ex[nx-1][1+-8*c3+c8+nx-1]-ex[nx-1][-8*c3+c8+nx-1]-ey[nx-1][-8*c3+c8+nx-1]) ;
        }
      }
      if ((c1 <= floord(96*c2-ny,32)) && (c2 >= max(ceild(8*c3-nx+ny+1,64),ceild(ny,64)))) {
        for (c9=max(64*c2-ny+1,8*c3);c9<=min(8*c3+7,64*c2+nx-ny);c9++) {
          hz[-64*c2+c9+ny-1][ny-1]=hz[-64*c2+c9+ny-1][ny-1]-((double)(7))/10*(ey[1+-64*c2+c9+ny-1][ny-1]+ex[-64*c2+c9+ny-1][1+ny-1]-ex[-64*c2+c9+ny-1][ny-1]-ey[-64*c2+c9+ny-1][ny-1]) ;
        }
      }
      if (c1 == c2+c3) {
        for (c7=max(max(32*c3,64*c2-ny+1),0);c7<=min(min(64*c2-1,8*c3+6),64*c2-ny+63);c7++) {
          for (c8=64*c2;c8<=c7+ny-1;c8++) {
            ey[0][-c7+c8]=c7 ;
            ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
            for (c9=c7+1;c9<=8*c3+7;c9++) {
              ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
              ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
              hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
            }
          }
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
          }
        }
      }
      if (c1 == c2+c3) {
        for (c7=max(max(64*c2,32*c3),0);c7<=min(8*c3+6,64*c2-ny+63);c7++) {
          ey[0][0]=c7 ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
          }
          for (c8=c7+1;c8<=c7+ny-1;c8++) {
            ey[0][-c7+c8]=c7 ;
            ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
            for (c9=c7+1;c9<=8*c3+7;c9++) {
              ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
              ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
              hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
            }
          }
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
          }
        }
      }
      if (c1 == c2+c3) {
        for (c7=max(max(0,32*c3),64*c2-ny+64);c7<=min(8*c3+6,64*c2-1);c7++) {
          for (c8=64*c2;c8<=64*c2+63;c8++) {
            ey[0][-c7+c8]=c7 ;
            ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
            for (c9=c7+1;c9<=8*c3+7;c9++) {
              ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
              ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
              hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
            }
          }
        }
      }
      if (c1 == c2+c3) {
        for (c7=max(max(max(64*c2,0),32*c3),64*c2-ny+64);c7<=min(8*c3+6,64*c2+62);c7++) {
          ey[0][0]=c7 ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
          }
          for (c8=c7+1;c8<=64*c2+63;c8++) {
            ey[0][-c7+c8]=c7 ;
            ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
            for (c9=c7+1;c9<=8*c3+7;c9++) {
              ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
              ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
              hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
            }
          }
        }
      }
      for (c7=max(max(max(64*c2-ny+1,8*c3-nx+1),0),32*c1-32*c2);c7<=min(min(min(min(64*c2-1,8*c3-nx+7),tmax-1),32*c1-32*c2+31),64*c2-ny+63);c7++) {
        for (c8=64*c2;c8<=c7+ny-1;c8++) {
          for (c9=8*c3;c9<=c7+nx-1;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
          hz[nx-1][-c7+c8-1]=hz[nx-1][-c7+c8-1]-((double)(7))/10*(ey[1+nx-1][-c7+c8-1]+ex[nx-1][1+-c7+c8-1]-ex[nx-1][-c7+c8-1]-ey[nx-1][-c7+c8-1]) ;
        }
        for (c9=8*c3;c9<=c7+nx;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
      for (c7=max(max(max(64*c2,8*c3-nx+1),0),32*c1-32*c2);c7<=min(min(min(8*c3-nx+7,tmax-1),32*c1-32*c2+31),64*c2-ny+63);c7++) {
        for (c9=8*c3;c9<=c7+nx-1;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=c7+ny-1;c8++) {
          for (c9=8*c3;c9<=c7+nx-1;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
          hz[nx-1][-c7+c8-1]=hz[nx-1][-c7+c8-1]-((double)(7))/10*(ey[1+nx-1][-c7+c8-1]+ex[nx-1][1+-c7+c8-1]-ex[nx-1][-c7+c8-1]-ey[nx-1][-c7+c8-1]) ;
        }
        for (c9=8*c3;c9<=c7+nx;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
      for (c7=max(max(max(32*c1-32*c2,0),8*c3-nx+1),64*c2-ny+64);c7<=min(min(min(32*c1-32*c2+31,tmax-1),64*c2-1),8*c3-nx+7);c7++) {
        for (c8=64*c2;c8<=64*c2+63;c8++) {
          for (c9=8*c3;c9<=c7+nx-1;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
          hz[nx-1][-c7+c8-1]=hz[nx-1][-c7+c8-1]-((double)(7))/10*(ey[1+nx-1][-c7+c8-1]+ex[nx-1][1+-c7+c8-1]-ex[nx-1][-c7+c8-1]-ey[nx-1][-c7+c8-1]) ;
        }
      }
      for (c7=max(max(max(8*c3-nx+8,64*c2-ny+1),0),32*c1-32*c2);c7<=min(min(min(min(64*c2-1,8*c3-1),tmax-1),32*c1-32*c2+31),64*c2-ny+63);c7++) {
        for (c8=64*c2;c8<=c7+ny-1;c8++) {
          for (c9=8*c3;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
        for (c9=8*c3;c9<=8*c3+7;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
      for (c7=max(max(max(max(64*c2,32*c1-32*c2),0),8*c3-nx+1),64*c2-ny+64);c7<=min(min(min(32*c1-32*c2+31,64*c2+62),tmax-1),8*c3-nx+7);c7++) {
        for (c9=8*c3;c9<=c7+nx-1;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=64*c2+63;c8++) {
          for (c9=8*c3;c9<=c7+nx-1;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
          hz[nx-1][-c7+c8-1]=hz[nx-1][-c7+c8-1]-((double)(7))/10*(ey[1+nx-1][-c7+c8-1]+ex[nx-1][1+-c7+c8-1]-ex[nx-1][-c7+c8-1]-ey[nx-1][-c7+c8-1]) ;
        }
      }
      for (c7=max(max(max(64*c2,8*c3-nx+8),0),32*c1-32*c2);c7<=min(min(min(8*c3-1,tmax-1),32*c1-32*c2+31),64*c2-ny+63);c7++) {
        for (c9=8*c3;c9<=8*c3+7;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=c7+ny-1;c8++) {
          for (c9=8*c3;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
        for (c9=8*c3;c9<=8*c3+7;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
      for (c7=max(max(8*c3,64*c2-ny+1),32*c1-32*c2);c7<=min(min(min(min(min(64*c2-1,tmax-1),32*c1-32*c2+31),8*c3+6),32*c3-1),64*c2-ny+63);c7++) {
        for (c8=64*c2;c8<=c7+ny-1;c8++) {
          ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
        for (c9=c7+1;c9<=8*c3+7;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
/*@ begin Loop(
transform Composite(
permut = [['c7', 'c9', 'c8']],
  regtile = (['c7', 'c8', 'c9'],[3, 1, 1]),
  scalarreplace = (True, 'double'),
  vector = (True, ['ivdep','vector always']))

      for (c7=max(max(max(32*c1-32*c2,0),8*c3-nx+8),64*c2-ny+64);c7<=min(min(min(tmax-1,32*c1-32*c2+31),64*c2-1),8*c3-1);c7++) {
        for (c8=64*c2;c8<=64*c2+63;c8++) {
          for (c9=8*c3;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
      }

) @*/{
  for (c7t=max(max(max(32*c1-32*c2,0),8*c3-nx+8),64*c2-ny+64); c7t<=min(min(min(tmax-1,32*c1-32*c2+31),64*c2-1),8*c3-1)-2; c7t=c7t+3) {
    for (c9=8*c3; c9<=8*c3+7; c9++ ) {
      register int cbv_1, cbv_2;
      cbv_1=64*c2;
      cbv_2=64*c2+63;
#pragma ivdep
#pragma vector always
      for (c8=cbv_1; c8<=cbv_2; c8++ ) {
        double scv_1, scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8;
        double scv_9, scv_10, scv_11, scv_12;
        scv_1=hz[-(c7t+2)+c9-1][-(c7t+2)+c8-1];
        scv_2=hz[-c7t+c9][-c7t+c8];
        scv_3=ey[-(c7t+2)+c9][-(c7t+2)+c8];
        scv_4=hz[-(c7t+2)+c9][-(c7t+2)+c8];
        scv_5=hz[-c7t+c9-1][-c7t+c8-1];
        scv_6=ex[-c7t+c9][-c7t+c8];
        scv_7=ex[-(c7t+1)+c9][-(c7t+1)+c8];
        scv_8=ey[-c7t+c9][-c7t+c8];
        scv_9=ex[-(c7t+2)+c9][-(c7t+2)+c8];
        scv_10=ey[-(c7t+1)+c9][-(c7t+1)+c8];
        scv_11=hz[-(c7t+1)+c9][-(c7t+1)+c8];
        scv_12=hz[-(c7t+1)+c9-1][-(c7t+1)+c8-1];
        scv_8=scv_8-((double)(1))/2*(scv_2-hz[-c7t+c9-1][-c7t+c8]);
        scv_10=scv_10-((double)(1))/2*(scv_11-hz[-(c7t+1)+c9-1][-(c7t+1)+c8]);
        scv_3=scv_3-((double)(1))/2*(scv_4-hz[-(c7t+2)+c9-1][-(c7t+2)+c8]);
        scv_6=scv_6-((double)(1))/2*(scv_2-hz[-c7t+c9][-c7t+c8-1]);
        scv_7=scv_7-((double)(1))/2*(scv_11-hz[-(c7t+1)+c9][-(c7t+1)+c8-1]);
        scv_9=scv_9-((double)(1))/2*(scv_4-hz[-(c7t+2)+c9][-(c7t+2)+c8-1]);
        scv_5=scv_5-((double)(7))/10*(ey[1+-c7t+c9-1][-c7t+c8-1]+ex[-c7t+c9-1][1+-c7t+c8-1]-ex[-c7t+c9-1][-c7t+c8-1]-ey[-c7t+c9-1][-c7t+c8-1]);
        scv_12=scv_12-((double)(7))/10*(ey[1+-(c7t+1)+c9-1][-(c7t+1)+c8-1]+ex[-(c7t+1)+c9-1][1+-(c7t+1)+c8-1]-ex[-(c7t+1)+c9-1][-(c7t+1)+c8-1]-ey[-(c7t+1)+c9-1][-(c7t+1)+c8-1]);
        scv_1=scv_1-((double)(7))/10*(ey[1+-(c7t+2)+c9-1][-(c7t+2)+c8-1]+ex[-(c7t+2)+c9-1][1+-(c7t+2)+c8-1]-ex[-(c7t+2)+c9-1][-(c7t+2)+c8-1]-ey[-(c7t+2)+c9-1][-(c7t+2)+c8-1]);
        hz[-(c7t+2)+c9-1][-(c7t+2)+c8-1]=scv_1;
        ey[-(c7t+2)+c9][-(c7t+2)+c8]=scv_3;
        hz[-c7t+c9-1][-c7t+c8-1]=scv_5;
        ex[-c7t+c9][-c7t+c8]=scv_6;
        ex[-(c7t+1)+c9][-(c7t+1)+c8]=scv_7;
        ey[-c7t+c9][-c7t+c8]=scv_8;
        ex[-(c7t+2)+c9][-(c7t+2)+c8]=scv_9;
        ey[-(c7t+1)+c9][-(c7t+1)+c8]=scv_10;
        hz[-(c7t+1)+c9-1][-(c7t+1)+c8-1]=scv_12;
      }
    }
  }
  for (c7=c7t; c7<=min(min(min(tmax-1,32*c1-32*c2+31),64*c2-1),8*c3-1); c7=c7+1) {
    for (c9=8*c3; c9<=8*c3+7; c9++ ) {
      register int cbv_3, cbv_4;
      cbv_3=64*c2;
      cbv_4=64*c2+63;
#pragma ivdep
#pragma vector always
      for (c8=cbv_3; c8<=cbv_4; c8++ ) {
        double scv_13, scv_14, scv_15, scv_16;
        scv_13=ey[-c7+c9][-c7+c8];
        scv_14=hz[-c7+c9][-c7+c8];
        scv_15=ex[-c7+c9][-c7+c8];
        scv_16=hz[-c7+c9-1][-c7+c8-1];
        scv_13=scv_13-((double)(1))/2*(scv_14-hz[-c7+c9-1][-c7+c8]);
        scv_15=scv_15-((double)(1))/2*(scv_14-hz[-c7+c9][-c7+c8-1]);
        scv_16=scv_16-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]);
        ey[-c7+c9][-c7+c8]=scv_13;
        ex[-c7+c9][-c7+c8]=scv_15;
        hz[-c7+c9-1][-c7+c8-1]=scv_16;
      }
    }
  }
}
/*@ end @*/

      for (c7=max(max(64*c2,8*c3),32*c1-32*c2);c7<=min(min(min(min(tmax-1,32*c1-32*c2+31),8*c3+6),32*c3-1),64*c2-ny+63);c7++) {
        for (c9=c7+1;c9<=8*c3+7;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=c7+ny-1;c8++) {
          ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
        for (c9=c7+1;c9<=8*c3+7;c9++) {
          hz[-c7+c9-1][ny-1]=hz[-c7+c9-1][ny-1]-((double)(7))/10*(ey[1+-c7+c9-1][ny-1]+ex[-c7+c9-1][1+ny-1]-ex[-c7+c9-1][ny-1]-ey[-c7+c9-1][ny-1]) ;
        }
      }
      for (c7=max(max(max(max(64*c2,8*c3-nx+8),32*c1-32*c2),0),64*c2-ny+64);c7<=min(min(min(tmax-1,32*c1-32*c2+31),64*c2+62),8*c3-1);c7++) {
        for (c9=8*c3;c9<=8*c3+7;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=64*c2+63;c8++) {
          for (c9=8*c3;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
      }
      for (c7=max(max(32*c1-32*c2,8*c3),64*c2-ny+64);c7<=min(min(min(min(32*c3-1,8*c3+6),32*c1-32*c2+31),tmax-1),64*c2-1);c7++) {
        for (c8=64*c2;c8<=64*c2+63;c8++) {
          ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
      }
      for (c7=max(max(max(64*c2,32*c1-32*c2),8*c3),64*c2-ny+64);c7<=min(min(min(min(32*c3-1,8*c3+6),64*c2+62),32*c1-32*c2+31),tmax-1);c7++) {
        for (c9=c7+1;c9<=8*c3+7;c9++) {
          ey[-c7+c9][0]=ey[-c7+c9][0]-((double)(1))/2*(hz[-c7+c9][0]-hz[-c7+c9-1][0]) ;
        }
        for (c8=c7+1;c8<=64*c2+63;c8++) {
          ex[0][-c7+c8]=ex[0][-c7+c8]-((double)(1))/2*(hz[0][-c7+c8]-hz[0][-c7+c8-1]) ;
          for (c9=c7+1;c9<=8*c3+7;c9++) {
            ey[-c7+c9][-c7+c8]=ey[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9-1][-c7+c8]) ;
            ex[-c7+c9][-c7+c8]=ex[-c7+c9][-c7+c8]-((double)(1))/2*(hz[-c7+c9][-c7+c8]-hz[-c7+c9][-c7+c8-1]) ;
            hz[-c7+c9-1][-c7+c8-1]=hz[-c7+c9-1][-c7+c8-1]-((double)(7))/10*(ey[1+-c7+c9-1][-c7+c8-1]+ex[-c7+c9-1][1+-c7+c8-1]-ex[-c7+c9-1][-c7+c8-1]-ey[-c7+c9-1][-c7+c8-1]) ;
          }
        }
      }
      if ((-c1 == -c2-c3) && (c1 <= floord(72*c3-57,64))) {
        ey[0][0]=64*c1-64*c3+63 ;
        for (c9=64*c1-64*c3+64;c9<=8*c3+7;c9++) {
          ey[-64*c1+64*c3+c9-63][0]=ey[-64*c1+64*c3+c9-63][0]-((double)(1))/2*(hz[-64*c1+64*c3+c9-63][0]-hz[-64*c1+64*c3+c9-63 -1][0]) ;
        }
      }
      if ((-c1 == -c2-c3) && (c1 >= ceild(72*c2-7,8)) && (c1 <= floord(72*c2+55,8))) {
        ey[0][0]=8*c1-8*c2+7 ;
        for (c8=8*c1-8*c2+8;c8<=min(64*c2+63,8*c1-8*c2+ny+6);c8++) {
          ey[0][-8*c1+8*c2+c8-7]=8*c1-8*c2+7 ;
          ex[0][-8*c1+8*c2+c8-7]=ex[0][-8*c1+8*c2+c8-7]-((double)(1))/2*(hz[0][-8*c1+8*c2+c8-7]-hz[0][-8*c1+8*c2+c8-7 -1]) ;
        }
      }
      if ((-c1 == -c2-c3) && (c1 <= 9*c2-1)) {
        for (c8=64*c2;c8<=min(64*c2+63,8*c1-8*c2+ny+6);c8++) {
          ey[0][-8*c1+8*c2+c8-7]=8*c1-8*c2+7 ;
          ex[0][-8*c1+8*c2+c8-7]=ex[0][-8*c1+8*c2+c8-7]-((double)(1))/2*(hz[0][-8*c1+8*c2+c8-7]-hz[0][-8*c1+8*c2+c8-7 -1]) ;
        }
      }
      if ((-c1 == -9*c2-7) && (-8*c1 == -9*c3+7)) {
        if ((64*c1+119)%9 == 0) {
          if ((64*c1+119)%9 == 0) {
            if ((64*c1+119)%9 == 0) {
              if ((8*c1+7)%9 == 0) {
                if ((-7*c1-14)%9 == 0) {
                  if ((8*c1+7)%9 == 0) {
                    if ((-7*c1-14)%9 == 0) {
                      if ((64*c1+119)%9 == 0) {
                        ey[0][0]=(64*c1+119)/9 ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      if ((c1 >= 3*c2+1) && (c2 <= min(min(floord(8*c3-57,64),floord(c3-2,2)),floord(tmax-64,64)))) {
        for (c9=max(64*c2+64,8*c3);c9<=min(8*c3+7,64*c2+nx+62);c9++) {
          ey[-64*c2+c9-63][0]=ey[-64*c2+c9-63][0]-((double)(1))/2*(hz[-64*c2+c9-63][0]-hz[-64*c2+c9-63 -1][0]) ;
        }
      }
      if ((c1 >= ceild(4*c2+c3-3,4)) && (c2 >= ceild(8*c3-55,64)) && (c3 >= 1) && (c3 <= floord(tmax-8,8))) {
        for (c8=max(64*c2,8*c3+8);c8<=min(64*c2+63,8*c3+ny+6);c8++) {
          ex[0][-8*c3+c8-7]=ex[0][-8*c3+c8-7]-((double)(1))/2*(hz[0][-8*c3+c8-7]-hz[0][-8*c3+c8-7 -1]) ;
        }
      }
    }
  }
}



    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);

  return ((int) hz[0][0]); 

}
