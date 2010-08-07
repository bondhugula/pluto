/* Generated from ../../../git/cloog/test/thomasset.cloog by CLooG 0.14.0-72-gefe2fc2 gmp bits in 0.10s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(i,j) { hash(1); hash(i); hash(j); }
#define S2(i,j,k,p,q) { hash(2); hash(i); hash(j); hash(k); hash(p); hash(q); }

void test(int n)
{
  /* Scattering iterators. */
  int c1, c2;
  /* Original iterators. */
  int i, j, k, p, q;
  for (c1=0;c1<=floord(n-5,3);c1++) {
    for (i=max(3*c1+1,1);i<=3*c1+3;i++) {
      S1(i,c1) ;
    }
  }
  if (n == 1) {
    S1(1,0) ;
    for (k=0;k<=min(0,0);k++) {
      for (p=max(0,ceild(-3*k-1,3));p<=min(floord(-3*k+1,3),0);p++) {
        q = -k-p ;
        S2(1,1,k,p,-k-p) ;
      }
    }
  }
  if (n >= 2) {
    for (c1=max(0,ceild(n-4,3));c1<=0;c1++) {
      S1(1,c1) ;
      for (j=1;j<=min(n,3*c1-n+5);j++) {
        for (k=0;k<=floord(3*c1-j-n+4,3);k++) {
          for (p=ceild(n-2,3);p<=floord(3*c1-j-3*k+2,3);p++) {
            q = c1-k-p ;
            S2(1,j,k,p,c1-k-p) ;
          }
        }
      }
      for (i=2;i<=min(n,3*c1+3);i++) {
        S1(i,c1) ;
      }
      for (c2=1;c2<=n-1;c2++) {
        i = c2+1 ;
        for (j=1;j<=min(3*c1-n+5,n);j++) {
          for (k=0;k<=floord(3*c1-j-n+4,3);k++) {
            for (p=ceild(n-2,3);p<=floord(3*c1-j-3*k+2,3);p++) {
              q = c1-k-p ;
              S2(c2+1,j,k,p,c1-k-p) ;
            }
          }
        }
      }
    }
  }
  for (c1=max(1,ceild(n-4,3));c1<=floord(n-1,3);c1++) {
    for (j=1;j<=3*c1-n+5;j++) {
      for (k=0;k<=min(floord(3*c1-j-n+4,3),0);k++) {
        for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(3*c1-j-3*k+2,3),floord(n,3));p++) {
          q = c1-k-p ;
          S2(1,j,k,p,c1-k-p) ;
        }
      }
    }
    for (i=3*c1+1;i<=min(n,3*c1+3);i++) {
      S1(i,c1) ;
    }
    for (c2=1;c2<=n-1;c2++) {
      i = c2+1 ;
      for (j=1;j<=3*c1-n+5;j++) {
        for (k=0;k<=min(floord(3*c1-j-n+4,3),0);k++) {
          for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(3*c1-j-3*k+2,3),floord(n,3));p++) {
            q = c1-k-p ;
            S2(c2+1,j,k,p,c1-k-p) ;
          }
        }
      }
    }
  }
  if (n >= 1) {
    for (c1=ceild(n,3);c1<=floord(2*n+1,3);c1++) {
      for (c2=0;c2<=n-1;c2++) {
        i = c2+1 ;
        for (j=max(1,3*c1-n-1);j<=min(n,3*c1-n+5);j++) {
          for (k=max(ceild(3*c1-j-n,3),0);k<=min(floord(3*c1-j-n+4,3),0);k++) {
            for (p=max(ceild(n-2,3),ceild(3*c1-j-3*k,3));p<=min(floord(3*c1-j-3*k+2,3),floord(n,3));p++) {
              q = c1-k-p ;
              S2(c2+1,j,k,p,c1-k-p) ;
            }
          }
        }
      }
    }
  }
}
