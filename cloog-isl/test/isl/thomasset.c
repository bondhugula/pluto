/* Generated from ../../../git/cloog/test/thomasset.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.18s. */
if (n >= 1) {
  for (c1=0;c1<=floord(n-4,3);c1++) {
    for (i=3*c1+1;i<=3*c1+3;i++) {
      S1(i,c1) ;
    }
  }
  if ((n >= 2) && (n <= 3)) {
    S1(1,0) ;
    for (j=1;j<=-n+4;j++) {
      for (p=ceild(n-2,3);p<=floord(-j+2,3);p++) {
        S2(1,j,0,p,-p) ;
      }
    }
    for (i=2;i<=n;i++) {
      S1(i,0) ;
    }
    for (c2=1;c2<=n-1;c2++) {
      for (j=1;j<=-n+4;j++) {
        for (p=ceild(n-2,3);p<=floord(-j+2,3);p++) {
          S2(c2+1,j,0,p,-p) ;
        }
      }
    }
  }
  if (n == 1) {
    S1(1,0) ;
    S2(1,1,0,0,0) ;
  }
  for (c1=max(1,ceild(n-3,3));c1<=floord(n-1,3);c1++) {
    for (j=1;j<=3*c1-n+4;j++) {
      for (p=ceild(n-2,3);p<=floord(3*c1-j+2,3);p++) {
        S2(1,j,0,p,c1-p) ;
      }
    }
    for (i=3*c1+1;i<=n;i++) {
      S1(i,c1) ;
    }
    for (c2=1;c2<=n-1;c2++) {
      for (j=1;j<=3*c1-n+4;j++) {
        for (p=ceild(n-2,3);p<=floord(3*c1-j+2,3);p++) {
          S2(c2+1,j,0,p,c1-p) ;
        }
      }
    }
  }
  for (c1=ceild(n,3);c1<=floord(2*n,3);c1++) {
    for (c2=0;c2<=n-1;c2++) {
      for (j=max(1,3*c1-n);j<=min(n,3*c1-n+4);j++) {
        for (p=max(ceild(3*c1-j,3),ceild(n-2,3));p<=min(floord(n,3),floord(3*c1-j+2,3));p++) {
          S2(c2+1,j,0,p,c1-p) ;
        }
      }
    }
  }
}
