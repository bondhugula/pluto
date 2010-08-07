	#define S1(zT0,zT1,zT2,t,j)	ey[0][j]=t;
	#define S2(zT0,zT1,zT2,t,i,j)	ey[i][j]=ey[i][j]-0.5*(hz[i][j]-hz[i-1][j]);
	#define S3(zT0,zT1,zT2,t,i,j)	ex[i][j]=ex[i][j]-0.5*(hz[i][j]-hz[i][j-1]);
	#define S4(zT0,zT1,zT2,t,i,j)	hz[i][j]=hz[i][j]-0.7*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);

		int t1, t2, t3, t4, t5, t6;

	register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 2.74s. */
if ((nx >= 2) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+31,32),floord(tmax+ny-1,32));t2++) {
      for (t3=max(ceild(32*t2-ny-30,32),t1);t3<=min(min(floord(32*t1+nx+31,32),floord(32*t2+nx+30,32)),floord(tmax+nx-1,32));t3++) {
        if ((t1 == t3) && (t1 <= floord(32*t2-ny,32))) {
          for (t6=32*t2-ny+1;t6<=min(32*t1+31,32*t2-ny+nx);t6++) {
            S4(t1,t2,t1,32*t2-ny,-32*t2+t6+ny-1,ny-1);
          }
        }
        if ((t1 <= min(floord(32*t2-ny,32),t3-1)) && (t2 >= ceild(32*t3+ny-nx+1,32))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t2-ny+nx);t6++) {
            S4(t1,t2,t3,32*t2-ny,-32*t2+t6+ny-1,ny-1);
          }
        }
        if ((t1 <= min(floord(32*t2-ny,32),floord(32*t2-ny+nx-32,32))) && (32*t2 == 32*t3+ny-nx)) {
          if ((31*ny+nx)%32 == 0) {
            S4(t1,t2,(32*t2-ny+nx)/32,32*t2-ny,nx-1,ny-1);
          }
        }
        if ((t1 <= floord(32*t3-nx,32)) && (t2 <= floord(32*t3+ny-nx-1,32))) {
          for (t5=max(32*t2,32*t3-nx+1);t5<=min(32*t2+31,32*t3+ny-nx);t5++) {
            S4(t1,t2,t3,32*t3-nx,nx-1,-32*t3+t5+nx-1);
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(min(tmax-1,32*t1-nx+31),32*t1-ny+31);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t6=t4+1;t6<=t4+nx-1;t6++) {
              S2(t1,t1,t1,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                S2(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t1,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t1,t1,t4,nx-1,-t4+t5-1);
            }
            for (t6=t4+1;t6<=t4+nx;t6++) {
              S4(t1,t1,t1,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-ny+32);t4<=min(tmax-1,32*t1-nx+31);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t6=t4+1;t6<=t4+nx-1;t6++) {
              S2(t1,t1,t1,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                S2(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t1,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t1,t1,t4,nx-1,-t4+t5-1);
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-nx+32);t4<=min(tmax-1,32*t1-ny+31);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              S2(t1,t1,t1,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                S2(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t1,t4,-t4+t6-1,-t4+t5-1);
              }
            }
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              S4(t1,t1,t1,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t1-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              S2(t1,t1,t1,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                S2(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t1,t4,-t4+t6-1,-t4+t5-1);
              }
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+1);t4<=min(min(tmax-1,32*t1-nx+31),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                S2(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t1,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t2,t1,t4,nx-1,-t4+t5-1);
            }
            for (t6=t4+1;t6<=t4+nx;t6++) {
              S4(t1,t2,t1,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+32);t4<=min(tmax-1,32*t1-nx+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                S2(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t1,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t2,t1,t4,nx-1,-t4+t5-1);
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t2-ny+1);t4<=min(min(32*t1+30,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                S2(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t1,t4,-t4+t6-1,-t4+t5-1);
              }
            }
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              S4(t1,t2,t1,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t2-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                S2(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t1,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t1,t4,-t4+t6-1,-t4+t5-1);
              }
            }
          }
        }
        if ((t1 == t3) && (t1 <= min(floord(tmax-32,32),t2-1))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+30);t5++) {
            S1(t1,t2,t1,32*t1+31,-32*t1+t5-31);
            S3(t1,t2,t1,32*t1+31,0,-32*t1+t5-31);
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          S1(t1,t1,t1,32*t1+31,0);
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+1);t4<=min(min(tmax-1,32*t1-ny+31),32*t3-nx+31);t4++) {
            for (t6=32*t3;t6<=t4+nx-1;t6++) {
              S2(t1,t1,t3,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                S2(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t3,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t1,t3,t4,nx-1,-t4+t5-1);
            }
            for (t6=32*t3;t6<=t4+nx;t6++) {
              S4(t1,t1,t3,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(max(32*t1,32*t1-ny+32),32*t3-nx+1);t4<=min(min(32*t1+30,tmax-1),32*t3-nx+31);t4++) {
            for (t6=32*t3;t6<=t4+nx-1;t6++) {
              S2(t1,t1,t3,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                S2(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t3,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t1,t3,t4,nx-1,-t4+t5-1);
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+32);t4<=min(tmax-1,32*t1-ny+31);t4++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              S2(t1,t1,t3,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                S2(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t3,t4,-t4+t6-1,-t4+t5-1);
              }
            }
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              S4(t1,t1,t3,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(max(32*t1,32*t1-ny+32),32*t3-nx+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              S2(t1,t1,t3,t4,-t4+t6,0);
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                S2(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t1,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t1,t3,t4,-t4+t6-1,-t4+t5-1);
              }
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+1),32*t3-nx+1);t4<=min(min(min(32*t1+31,tmax-1),32*t2-ny+31),32*t3-nx+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                S2(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t3,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t2,t3,t4,nx-1,-t4+t5-1);
            }
            for (t6=32*t3;t6<=t4+nx;t6++) {
              S4(t1,t2,t3,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+32),32*t3-nx+1);t4<=min(min(32*t1+31,tmax-1),32*t3-nx+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                S2(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t3,t4,-t4+t6-1,-t4+t5-1);
              }
              S4(t1,t2,t3,t4,nx-1,-t4+t5-1);
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+1),32*t3-nx+32);t4<=min(min(32*t1+31,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                S2(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t3,t4,-t4+t6-1,-t4+t5-1);
              }
            }
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              S4(t1,t2,t3,t4,-t4+t6-1,ny-1);
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+32),32*t3-nx+32);t4<=min(32*t1+31,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                S2(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S3(t1,t2,t3,t4,-t4+t6,-t4+t5);
                S4(t1,t2,t3,t4,-t4+t6-1,-t4+t5-1);
              }
            }
          }
        }
        if ((t1 == t2) && (t1 <= min(floord(tmax-32,32),t3-1))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+30);t6++) {
            S2(t1,t1,t3,32*t1+31,-32*t1+t6-31,0);
          }
        }
      }
    }
  }
}
if ((nx >= 2) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(tmax,32),t1+1);t2++) {
      for (t3=t2;t3<=min(min(floord(32*t1+nx+31,32),floord(32*t2+nx+30,32)),floord(tmax+nx-1,32));t3++) {
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(32*t1+30,tmax-1);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t6=t4+1;t6<=min(32*t1+31,t4+nx-1);t6++) {
              S2(t1,t1,t1,t4,-t4+t6,0);
            }
            for (t6=t4+1;t6<=min(32*t1+31,t4+nx);t6++) {
              S4(t1,t1,t1,t4,-t4+t6-1,0);
            }
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          S1(t1,t1,t1,32*t1+31,0);
        }
        if ((t1 == t2) && (t1 <= floord(32*t3-nx,32))) {
          S4(t1,t1,t3,32*t3-nx,nx-1,0);
        }
        if ((t1 == t2-1) && (32*t1 == 32*t3-nx-31)) {
          if ((nx+31)%32 == 0) {
            S4(t1,t1+1,(32*t1+nx+31)/32,32*t1+31,nx-1,0);
          }
        }
        if ((t1 == t2-1) && (t1 >= ceild(32*t3-nx-30,32))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+31);t6++) {
            S4(t1,t1+1,t3,32*t1+31,-32*t1+t6-32,0);
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+1);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t6=32*t3;t6<=min(32*t3+31,t4+nx-1);t6++) {
              S2(t1,t1,t3,t4,-t4+t6,0);
            }
            for (t6=32*t3;t6<=min(32*t3+31,t4+nx);t6++) {
              S4(t1,t1,t3,t4,-t4+t6-1,0);
            }
          }
        }
        if ((t1 == t2) && (t1 <= min(floord(tmax-32,32),t3-1))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+30);t6++) {
            S2(t1,t1,t3,32*t1+31,-32*t1+t6-31,0);
          }
        }
      }
    }
  }
}
if ((nx == 1) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+31,32),floord(tmax+ny-1,32));t2++) {
      for (t3=max(ceild(32*t2-ny-30,32),t1);t3<=min(min(floord(tmax,32),t2),t1+1);t3++) {
        if ((t1 == t3) && (t1 <= floord(32*t2-ny,32))) {
          S4(t1,t2,t1,32*t2-ny,0,ny-1);
        }
        if ((t1 == t3-1) && (32*t1 == 32*t2-ny-31)) {
          if ((ny+31)%32 == 0) {
            S4(t1,(32*t1+ny+31)/32,t1+1,32*t1+31,0,ny-1);
          }
        }
        if ((t1 == t3-1) && (t1 >= ceild(32*t2-ny-30,32))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+31);t5++) {
            S4(t1,t2,t1+1,32*t1+31,0,-32*t1+t5-32);
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(tmax-1,32*t1-ny+31);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              S4(t1,t1,t1,t4,0,-t4+t5-1);
            }
            S4(t1,t1,t1,t4,0,ny-1);
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            S1(t1,t1,t1,t4,0);
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              S1(t1,t1,t1,t4,-t4+t5);
              S3(t1,t1,t1,t4,0,-t4+t5);
              S4(t1,t1,t1,t4,0,-t4+t5-1);
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+1);t4<=min(min(32*t1+30,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              S4(t1,t2,t1,t4,0,-t4+t5-1);
            }
            S4(t1,t2,t1,t4,0,ny-1);
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              S1(t1,t2,t1,t4,-t4+t5);
              S3(t1,t2,t1,t4,0,-t4+t5);
              S4(t1,t2,t1,t4,0,-t4+t5-1);
            }
          }
        }
        if ((t1 == t3) && (t1 <= min(floord(tmax-32,32),t2-1))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+30);t5++) {
            S1(t1,t2,t1,32*t1+31,-32*t1+t5-31);
            S3(t1,t2,t1,32*t1+31,0,-32*t1+t5-31);
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          S1(t1,t1,t1,32*t1+31,0);
        }
      }
    }
  }
}
if ((nx <= 0) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+30,32),floord(tmax+ny-2,32));t2++) {
      for (t4=max(32*t1,32*t2-ny+1);t4<=min(32*t1+31,tmax-1);t4++) {
        for (t5=max(32*t2,t4);t5<=min(32*t2+31,t4+ny-1);t5++) {
          S1(t1,t2,t1,t4,-t4+t5);
        }
      }
    }
  }
}
if ((nx == 1) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(tmax,32),t1+1);t2++) {
      if (t1 == t2) {
        for (t4=32*t1;t4<=min(32*t1+30,tmax-1);t4++) {
          S1(t1,t1,t1,t4,0);
          S4(t1,t1,t1,t4,0,0);
        }
      }
      if ((t1 == t2) && (t1 <= floord(tmax-32,32))) {
        S1(t1,t1,t1,32*t1+31,0);
      }
      if (t1 == t2-1) {
        S4(t1,t1+1,t1+1,32*t1+31,0,0);
      }
    }
  }
}
if ((nx <= 0) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t4=32*t1;t4<=min(32*t1+31,tmax-1);t4++) {
      S1(t1,t1,t1,t4,0);
    }
  }
}
/* End of CLooG code */
