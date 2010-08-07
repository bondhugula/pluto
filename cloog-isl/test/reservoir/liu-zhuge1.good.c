/* Generated from ../../../git/cloog/test/./reservoir/liu-zhuge1.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.05s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j) { hash(2); hash(i); hash(j); }
#define S3(i,j) { hash(3); hash(i); hash(j); }

void test(int M, int N)
{
  /* Scattering iterators. */
  int c2, c4;
  /* Original iterators. */
  int i, j;
  if ((M >= 0) && (N >= 0)) {
    for (c2=-4;c2<=min(-1,3*M+N-4);c2++) {
      for (c4=max(0,c2-3*M+4);c4<=min(c2+4,N);c4++) {
        if ((c2+2*c4+1)%3 == 0) {
          i = (c2-c4+4)/3 ;
          S1((c2-c4+4)/3,c4) ;
        }
      }
    }
  }
  if ((M <= 1) && (M >= 0)) {
    for (c2=0;c2<=3*M+N-4;c2++) {
      for (c4=max(c2-3*M,0);c4<=c2;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S2((c2-c4)/3,c4) ;
        }
      }
      for (c4=c2-3*M+4;c4<=min(c2+4,N);c4++) {
        if ((c2+2*c4+1)%3 == 0) {
          i = (c2-c4+4)/3 ;
          S1((c2-c4+4)/3,c4) ;
        }
      }
      for (c4=max(0,c2-3*M);c4<=c2;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S3((c2-c4)/3,c4) ;
        }
      }
    }
  }
  for (c2=0;c2<=min(3*M-4,N-1);c2++) {
    for (c4=0;c4<=c2;c4++) {
      if ((c2+2*c4)%3 == 0) {
        i = (c2-c4)/3 ;
        S2((c2-c4)/3,c4) ;
      }
      if ((c2+2*c4+1)%3 == 0) {
        i = (c2-c4+4)/3 ;
        S1((c2-c4+4)/3,c4) ;
      }
    }
    for (c4=c2+1;c4<=min(c2+4,N);c4++) {
      if ((c2+2*c4+1)%3 == 0) {
        i = (c2-c4+4)/3 ;
        S1((c2-c4+4)/3,c4) ;
      }
    }
    for (c4=0;c4<=c2;c4++) {
      if ((c2+2*c4)%3 == 0) {
        i = (c2-c4)/3 ;
        S3((c2-c4)/3,c4) ;
      }
    }
  }
  if (M >= 2) {
    for (c2=3*M-3;c2<=N-1;c2++) {
      for (c4=max(c2-3*M,0);c4<=c2-3*M+3;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S2((c2-c4)/3,c4) ;
        }
      }
      for (c4=c2-3*M+4;c4<=c2;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S2((c2-c4)/3,c4) ;
        }
        if ((c2+2*c4+1)%3 == 0) {
          i = (c2-c4+4)/3 ;
          S1((c2-c4+4)/3,c4) ;
        }
      }
      for (c4=c2+1;c4<=min(c2+4,N);c4++) {
        if ((c2+2*c4+1)%3 == 0) {
          i = (c2-c4+4)/3 ;
          S1((c2-c4+4)/3,c4) ;
        }
      }
      for (c4=max(0,c2-3*M);c4<=c2;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S3((c2-c4)/3,c4) ;
        }
      }
    }
  }
  if (N >= 0) {
    for (c2=N;c2<=3*M-4;c2++) {
      for (c4=0;c4<=N;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S2((c2-c4)/3,c4) ;
        }
        if ((c2+2*c4+1)%3 == 0) {
          i = (c2-c4+4)/3 ;
          S1((c2-c4+4)/3,c4) ;
        }
      }
      for (c4=0;c4<=N;c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S3((c2-c4)/3,c4) ;
        }
      }
    }
  }
  for (c2=max(3*M-3,N);c2<=3*M+N-4;c2++) {
    for (c4=max(c2-3*M,0);c4<=c2-3*M+3;c4++) {
      if ((c2+2*c4)%3 == 0) {
        i = (c2-c4)/3 ;
        S2((c2-c4)/3,c4) ;
      }
    }
    for (c4=c2-3*M+4;c4<=N;c4++) {
      if ((c2+2*c4)%3 == 0) {
        i = (c2-c4)/3 ;
        S2((c2-c4)/3,c4) ;
      }
      if ((c2+2*c4+1)%3 == 0) {
        i = (c2-c4+4)/3 ;
        S1((c2-c4+4)/3,c4) ;
      }
    }
    for (c4=max(0,c2-3*M);c4<=N;c4++) {
      if ((c2+2*c4)%3 == 0) {
        i = (c2-c4)/3 ;
        S3((c2-c4)/3,c4) ;
      }
    }
  }
  if ((M >= 0) && (N >= 0)) {
    for (c2=max(3*M+N-3,0);c2<=3*M+N;c2++) {
      for (c4=max(0,c2-3*M);c4<=min(c2,N);c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S2((c2-c4)/3,c4) ;
        }
      }
      for (c4=max(0,c2-3*M);c4<=min(c2,N);c4++) {
        if ((c2+2*c4)%3 == 0) {
          i = (c2-c4)/3 ;
          S3((c2-c4)/3,c4) ;
        }
      }
    }
  }
}
