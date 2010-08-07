/* Generated from ../../../git/cloog/test/./reservoir/mg-interp.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 1.01s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j,k) { hash(1); hash(i); hash(j); hash(k); }
#define S2(i,j,k) { hash(2); hash(i); hash(j); hash(k); }
#define S3(i,j,k) { hash(3); hash(i); hash(j); hash(k); }
#define S4(i,j,k) { hash(4); hash(i); hash(j); hash(k); }
#define S5(i,j,k) { hash(5); hash(i); hash(j); hash(k); }
#define S6(i,j,k) { hash(6); hash(i); hash(j); hash(k); }
#define S7(i,j,k) { hash(7); hash(i); hash(j); hash(k); }
#define S8(i,j,k) { hash(8); hash(i); hash(j); hash(k); }
#define S9(i,j,k) { hash(9); hash(i); hash(j); hash(k); }
#define S10(i,j,k) { hash(10); hash(i); hash(j); hash(k); }
#define S11(i,j,k) { hash(11); hash(i); hash(j); hash(k); }
#define S12(i,j,k) { hash(12); hash(i); hash(j); hash(k); }
#define S13(i,j,k) { hash(13); hash(i); hash(j); hash(k); }
#define S14(i,j,k) { hash(14); hash(i); hash(j); hash(k); }
#define S15(i,j,k) { hash(15); hash(i); hash(j); hash(k); }

void test(int M, int N, int O, int P, int Q, int R, int S, int T, int U)
{
  /* Scattering iterators. */
  int c2, c4, c6;
  /* Original iterators. */
  int i, j, k;
  if ((M >= 2) && (N >= 4)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c6=1;c6<=M;c6++) {
        S1(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S6(c2,1,c6) ;
        S7(c2,1,c6) ;
      }
      for (c6=1;c6<=M;c6++) {
        S3(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S1(c2,2,c6) ;
      }
      S1(c2,2,M) ;
      for (c6=1;c6<=M-1;c6++) {
        S6(c2,2,c6) ;
        S7(c2,2,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S11(c2,1,c6) ;
      }
      for (c4=3;c4<=2*N-5;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          if ((c4+1)%2 == 0) {
            j = (c4-1)/2 ;
            S10(c2,(c4-1)/2,c6) ;
          }
        }
        for (c6=1;c6<=M;c6++) {
          if ((c4+1)%2 == 0) {
            j = (c4+1)/2 ;
            S3(c2,(c4+1)/2,c6) ;
          }
        }
        for (c6=1;c6<=M-1;c6++) {
          if (c4%2 == 0) {
            j = (c4+2)/2 ;
            S6(c2,(c4+2)/2,c6) ;
            S7(c2,(c4+2)/2,c6) ;
          }
          if ((c4+1)%2 == 0) {
            j = (c4+3)/2 ;
            S1(c2,(c4+3)/2,c6) ;
          }
        }
        if ((c4+1)%2 == 0) {
          j = (c4+3)/2 ;
          S1(c2,(c4+3)/2,M) ;
        }
        for (c6=1;c6<=M-1;c6++) {
          if (c4%2 == 0) {
            S11(c2,c4/2,c6) ;
          }
        }
      }
      c4 = 2*N-4 ;
      for (c6=1;c6<=M-1;c6++) {
        j = N-1 ;
        S6(c2,N-1,c6) ;
        S7(c2,N-1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        j = N-2 ;
        S11(c2,N-2,c6) ;
      }
      c4 = 2*N-3 ;
      for (c6=1;c6<=M-1;c6++) {
        j = N-2 ;
        S10(c2,N-2,c6) ;
      }
      for (c6=1;c6<=M;c6++) {
        j = N-1 ;
        S3(c2,N-1,c6) ;
      }
      c4 = 2*N-2 ;
      for (c6=1;c6<=M-1;c6++) {
        j = N-1 ;
        S11(c2,N-1,c6) ;
      }
      c4 = 2*N-1 ;
      for (c6=1;c6<=M-1;c6++) {
        j = N-1 ;
        S10(c2,N-1,c6) ;
      }
    }
  }
  if ((M >= 2) && (N == 3)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c6=1;c6<=M;c6++) {
        S1(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S6(c2,1,c6) ;
        S7(c2,1,c6) ;
      }
      for (c6=1;c6<=M;c6++) {
        S3(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S1(c2,2,c6) ;
      }
      S1(c2,2,M) ;
      for (c6=1;c6<=M-1;c6++) {
        S6(c2,2,c6) ;
        S7(c2,2,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S11(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S10(c2,1,c6) ;
      }
      for (c6=1;c6<=M;c6++) {
        S3(c2,2,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S11(c2,2,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S10(c2,2,c6) ;
      }
    }
  }
  if ((M >= 2) && (N == 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c6=1;c6<=M;c6++) {
        S1(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S6(c2,1,c6) ;
        S7(c2,1,c6) ;
      }
      for (c6=1;c6<=M;c6++) {
        S3(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S11(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S10(c2,1,c6) ;
      }
    }
  }
  if ((M == 1) && (N >= 3)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=-1;c4<=0;c4++) {
        if ((c4+1)%2 == 0) {
          j = (c4+3)/2 ;
          S1(c2,(c4+3)/2,1) ;
        }
      }
      for (c4=1;c4<=2*N-5;c4++) {
        if ((c4+1)%2 == 0) {
          j = (c4+1)/2 ;
          S3(c2,(c4+1)/2,1) ;
        }
        if ((c4+1)%2 == 0) {
          j = (c4+3)/2 ;
          S1(c2,(c4+3)/2,1) ;
        }
      }
      for (c4=2*N-4;c4<=2*N-3;c4++) {
        if ((c4+1)%2 == 0) {
          j = (c4+1)/2 ;
          S3(c2,(c4+1)/2,1) ;
        }
      }
    }
  }
  if ((M == 1) && (N == 2)) {
    for (c2=1;c2<=O-1;c2++) {
      S1(c2,1,1) ;
      S3(c2,1,1) ;
    }
  }
  if ((M >= 2) && (N >= 3)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c6=1;c6<=M;c6++) {
        S2(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S8(c2,1,c6) ;
      }
      for (c4=3;c4<=2*N-2;c4++) {
        for (c6=1;c6<=M;c6++) {
          if (c4%2 == 0) {
            S2(c2,c4/2,c6) ;
          }
        }
        for (c6=1;c6<=M-1;c6++) {
          if (c4%2 == 0) {
            S8(c2,c4/2,c6) ;
          }
        }
        for (c6=1;c6<=M-1;c6++) {
          if ((c4+1)%2 == 0) {
            j = (c4-1)/2 ;
            S9(c2,(c4-1)/2,c6) ;
          }
        }
      }
      c4 = 2*N-1 ;
      for (c6=1;c6<=M-1;c6++) {
        j = N-1 ;
        S9(c2,N-1,c6) ;
      }
    }
  }
  if ((M >= 2) && (N == 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c6=1;c6<=M;c6++) {
        S2(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S8(c2,1,c6) ;
      }
      for (c6=1;c6<=M-1;c6++) {
        S9(c2,1,c6) ;
      }
    }
  }
  if ((M == 1) && (N >= 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=2;c4<=2*N-2;c4++) {
        if (c4%2 == 0) {
          S2(c2,c4/2,1) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S4(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= 2)) {
    for (c2=1;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S5(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= P+1) && (N >= Q+1)) {
    for (c2=R;c2<=O-1;c2++) {
      for (c4=Q;c4<=N-1;c4++) {
        for (c6=P;c6<=M-1;c6++) {
          S12(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= Q+1)) {
    for (c2=R;c2<=O-1;c2++) {
      for (c4=Q;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S13(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= P+1) && (N >= 2)) {
    for (c2=R;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=P;c6<=M-1;c6++) {
          S14(c2,c4,c6) ;
        }
      }
    }
  }
  if ((M >= 2) && (N >= 2)) {
    for (c2=R;c2<=O-1;c2++) {
      for (c4=1;c4<=N-1;c4++) {
        for (c6=1;c6<=M-1;c6++) {
          S15(c2,c4,c6) ;
        }
      }
    }
  }
}
