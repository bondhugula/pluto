/* Generated from ../../../git/cloog/test/classen.cloog by CLooG 0.14.0-76-gfd78716 gmp bits in 1.73s. */
extern void hash(int);

/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#define S1(coordT1,coordP1,other1,other2) { hash(1); hash(coordT1); hash(coordP1); hash(other1); hash(other2); }
#define S2(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(2); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S3(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(3); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S4(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(4); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S5(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(5); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S6(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(6); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S7(coordT1,coordP1,other1,other2,other3,other4,other5,other6) { hash(7); hash(coordT1); hash(coordP1); hash(other1); hash(other2); hash(other3); hash(other4); hash(other5); hash(other6); }
#define S8(coordT1,coordP1) { hash(8); hash(coordT1); hash(coordP1); }

void test(int m)
{
  /* Scattering iterators. */
  int glT1, rp1, local1, local2;
  /* Original iterators. */
  int coordT1, coordP1, other1, other2, other3, other4, other5, other6;
  if (m >= 2) {
    S1(0,1,1,1) ;
    S2(0,1,1,1,1,1,2,1) ;
    S3(0,1,1,2,1,1,1,2) ;
    S4(0,1,2,2,1,1,2,2) ;
    S8(0,1) ;
  }
  if (m == 1) {
    S1(0,1,1,1) ;
    S8(0,1) ;
  }
  if (m >= 3) {
    S5(0,1,1,1,1,1,2,1) ;
    S1(1,1,2,1) ;
    S2(1,1,2,1,2,1,3,1) ;
    S3(1,1,2,2,2,1,2,2) ;
    S4(1,1,3,2,2,1,3,2) ;
    S6(0,1,1,2,1,1,1,2) ;
    S7(0,1,2,2,1,1,2,2) ;
    S1(1,2,1,2) ;
    S2(1,2,2,2,1,2,2,2) ;
    S3(1,2,2,3,1,2,1,3) ;
    S4(1,2,3,3,1,2,2,3) ;
    for (coordP1=1;coordP1<=2;coordP1++) {
      S8(1,coordP1) ;
    }
  }
  for (glT1=2;glT1<=m-2;glT1++) {
    coordT1 = glT1-1 ;
    other5 = glT1+1 ;
    S5(glT1-1,1,glT1,1,glT1,1,glT1+1,1) ;
    other1 = glT1+1 ;
    S1(glT1,1,glT1+1,1) ;
    local1 = glT1+1 ;
    other1 = glT1+1 ;
    other3 = glT1+1 ;
    other5 = glT1+2 ;
    S2(glT1,1,glT1+1,1,glT1+1,1,glT1+2,1) ;
    other3 = glT1+1 ;
    other5 = glT1+1 ;
    S3(glT1,1,glT1+1,2,glT1+1,1,glT1+1,2) ;
    other1 = glT1+2 ;
    other3 = glT1+1 ;
    other5 = glT1+2 ;
    S4(glT1,1,glT1+2,2,glT1+1,1,glT1+2,2) ;
    for (rp1=2;rp1<=glT1;rp1++) {
      local1 = glT1-rp1+1 ;
      coordT1 = glT1-1 ;
      other3 = glT1-rp1+1 ;
      other5 = glT1-rp1+2 ;
      S5(glT1-1,rp1,glT1,rp1,glT1-rp1+1,rp1,glT1-rp1+2,rp1) ;
      local1 = glT1-rp1+2 ;
      local2 = rp1-1 ;
      coordT1 = glT1-1 ;
      coordP1 = rp1-1 ;
      other3 = glT1-rp1+2 ;
      other4 = rp1-1 ;
      other5 = glT1-rp1+2 ;
      S6(glT1-1,rp1-1,glT1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+2,rp1) ;
      other1 = glT1+1 ;
      other3 = glT1-rp1+2 ;
      other4 = rp1-1 ;
      other5 = glT1-rp1+3 ;
      S7(glT1-1,rp1-1,glT1+1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+3,rp1) ;
      other1 = glT1-rp1+2 ;
      S1(glT1,rp1,glT1-rp1+2,rp1) ;
      local1 = glT1-rp1+2 ;
      other1 = glT1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+3 ;
      S2(glT1,rp1,glT1+1,rp1,glT1-rp1+2,rp1,glT1-rp1+3,rp1) ;
      other2 = rp1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+2 ;
      other6 = rp1+1 ;
      S3(glT1,rp1,glT1+1,rp1+1,glT1-rp1+2,rp1,glT1-rp1+2,rp1+1) ;
      other1 = glT1+2 ;
      other2 = rp1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+3 ;
      other6 = rp1+1 ;
      S4(glT1,rp1,glT1+2,rp1+1,glT1-rp1+2,rp1,glT1-rp1+3,rp1+1) ;
    }
    rp1 = glT1+1 ;
    coordT1 = glT1-1 ;
    other2 = glT1+1 ;
    other6 = glT1+1 ;
    S6(glT1-1,glT1,glT1,glT1+1,1,glT1,1,glT1+1) ;
    other1 = glT1+1 ;
    other2 = glT1+1 ;
    other6 = glT1+1 ;
    S7(glT1-1,glT1,glT1+1,glT1+1,1,glT1,2,glT1+1) ;
    coordP1 = glT1+1 ;
    other2 = glT1+1 ;
    S1(glT1,glT1+1,1,glT1+1) ;
    local2 = glT1+1 ;
    coordP1 = glT1+1 ;
    other1 = glT1+1 ;
    other2 = glT1+1 ;
    other4 = glT1+1 ;
    other6 = glT1+1 ;
    S2(glT1,glT1+1,glT1+1,glT1+1,1,glT1+1,2,glT1+1) ;
    other2 = glT1+2 ;
    other4 = glT1+1 ;
    other6 = glT1+2 ;
    S3(glT1,glT1+1,glT1+1,glT1+2,1,glT1+1,1,glT1+2) ;
    other1 = glT1+2 ;
    other2 = glT1+2 ;
    other4 = glT1+1 ;
    other6 = glT1+2 ;
    S4(glT1,glT1+1,glT1+2,glT1+2,1,glT1+1,2,glT1+2) ;
    for (coordP1=1;coordP1<=glT1+1;coordP1++) {
      S8(glT1,coordP1) ;
    }
  }
  if (m >= 3) {
    glT1 = m-1 ;
    local1 = m-1 ;
    coordT1 = m-2 ;
    other1 = m-1 ;
    other3 = m-1 ;
    S5(m-2,1,m-1,1,m-1,1,m,1) ;
    coordT1 = m-1 ;
    S1(m-1,1,m,1) ;
    coordT1 = m-1 ;
    S3(m-1,1,m,2,m,1,m,2) ;
    for (rp1=2;rp1<=m-1;rp1++) {
      local1 = -rp1+m ;
      coordT1 = m-2 ;
      other1 = m-1 ;
      other3 = -rp1+m ;
      other5 = -rp1+m+1 ;
      S5(m-2,rp1,m-1,rp1,-rp1+m,rp1,-rp1+m+1,rp1) ;
      local1 = -rp1+m+1 ;
      local2 = rp1-1 ;
      coordT1 = m-2 ;
      coordP1 = rp1-1 ;
      other1 = m-1 ;
      other3 = -rp1+m+1 ;
      other4 = rp1-1 ;
      other5 = -rp1+m+1 ;
      S6(m-2,rp1-1,m-1,rp1,-rp1+m+1,rp1-1,-rp1+m+1,rp1) ;
      other3 = -rp1+m+1 ;
      other4 = rp1-1 ;
      other5 = -rp1+m+2 ;
      S7(m-2,rp1-1,m,rp1,-rp1+m+1,rp1-1,-rp1+m+2,rp1) ;
      coordT1 = m-1 ;
      other1 = -rp1+m+1 ;
      S1(m-1,rp1,-rp1+m+1,rp1) ;
      local1 = -rp1+m+1 ;
      coordT1 = m-1 ;
      other3 = -rp1+m+1 ;
      other5 = -rp1+m+2 ;
      S2(m-1,rp1,m,rp1,-rp1+m+1,rp1,-rp1+m+2,rp1) ;
      other2 = rp1+1 ;
      other3 = -rp1+m+1 ;
      other5 = -rp1+m+1 ;
      other6 = rp1+1 ;
      S3(m-1,rp1,m,rp1+1,-rp1+m+1,rp1,-rp1+m+1,rp1+1) ;
      other1 = m+1 ;
      other2 = rp1+1 ;
      other3 = -rp1+m+1 ;
      other5 = -rp1+m+2 ;
      other6 = rp1+1 ;
      S4(m-1,rp1,m+1,rp1+1,-rp1+m+1,rp1,-rp1+m+2,rp1+1) ;
    }
    local2 = m-1 ;
    coordT1 = m-2 ;
    coordP1 = m-1 ;
    other1 = m-1 ;
    other4 = m-1 ;
    S6(m-2,m-1,m-1,m,1,m-1,1,m) ;
    other4 = m-1 ;
    S7(m-2,m-1,m,m,1,m-1,2,m) ;
    coordT1 = m-1 ;
    S1(m-1,m,1,m) ;
    coordT1 = m-1 ;
    S2(m-1,m,m,m,1,m,2,m) ;
    coordT1 = m-1 ;
    for (coordP1=1;coordP1<=m;coordP1++) {
      S8(m-1,coordP1) ;
    }
  }
  for (glT1=m;glT1<=2*m-4;glT1++) {
    rp1 = glT1-m+2 ;
    local1 = m-1 ;
    local2 = glT1-m+2 ;
    coordT1 = glT1-1 ;
    coordP1 = glT1-m+2 ;
    other2 = glT1-m+2 ;
    other3 = m-1 ;
    other4 = glT1-m+2 ;
    other6 = glT1-m+2 ;
    S5(glT1-1,glT1-m+2,glT1,glT1-m+2,m-1,glT1-m+2,m,glT1-m+2) ;
    local2 = glT1-m+1 ;
    coordT1 = glT1-1 ;
    coordP1 = glT1-m+1 ;
    other2 = glT1-m+2 ;
    other4 = glT1-m+1 ;
    other6 = glT1-m+2 ;
    S6(glT1-1,glT1-m+1,glT1,glT1-m+2,m,glT1-m+1,m,glT1-m+2) ;
    coordP1 = glT1-m+2 ;
    other2 = glT1-m+2 ;
    S1(glT1,glT1-m+2,m,glT1-m+2) ;
    local2 = glT1-m+2 ;
    coordP1 = glT1-m+2 ;
    other1 = glT1+1 ;
    other2 = glT1-m+3 ;
    other4 = glT1-m+2 ;
    other6 = glT1-m+3 ;
    S3(glT1,glT1-m+2,glT1+1,glT1-m+3,m,glT1-m+2,m,glT1-m+3) ;
    for (rp1=glT1-m+3;rp1<=m-1;rp1++) {
      local1 = glT1-rp1+1 ;
      coordT1 = glT1-1 ;
      other3 = glT1-rp1+1 ;
      other5 = glT1-rp1+2 ;
      S5(glT1-1,rp1,glT1,rp1,glT1-rp1+1,rp1,glT1-rp1+2,rp1) ;
      local1 = glT1-rp1+2 ;
      local2 = rp1-1 ;
      coordT1 = glT1-1 ;
      coordP1 = rp1-1 ;
      other3 = glT1-rp1+2 ;
      other4 = rp1-1 ;
      other5 = glT1-rp1+2 ;
      S6(glT1-1,rp1-1,glT1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+2,rp1) ;
      other1 = glT1+1 ;
      other3 = glT1-rp1+2 ;
      other4 = rp1-1 ;
      other5 = glT1-rp1+3 ;
      S7(glT1-1,rp1-1,glT1+1,rp1,glT1-rp1+2,rp1-1,glT1-rp1+3,rp1) ;
      other1 = glT1-rp1+2 ;
      S1(glT1,rp1,glT1-rp1+2,rp1) ;
      local1 = glT1-rp1+2 ;
      other1 = glT1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+3 ;
      S2(glT1,rp1,glT1+1,rp1,glT1-rp1+2,rp1,glT1-rp1+3,rp1) ;
      other2 = rp1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+2 ;
      other6 = rp1+1 ;
      S3(glT1,rp1,glT1+1,rp1+1,glT1-rp1+2,rp1,glT1-rp1+2,rp1+1) ;
      other1 = glT1+2 ;
      other2 = rp1+1 ;
      other3 = glT1-rp1+2 ;
      other5 = glT1-rp1+3 ;
      other6 = rp1+1 ;
      S4(glT1,rp1,glT1+2,rp1+1,glT1-rp1+2,rp1,glT1-rp1+3,rp1+1) ;
    }
    local1 = glT1-m+1 ;
    coordT1 = glT1-1 ;
    other3 = glT1-m+1 ;
    other5 = glT1-m+2 ;
    S5(glT1-1,m,glT1,m,glT1-m+1,m,glT1-m+2,m) ;
    local1 = glT1-m+2 ;
    local2 = m-1 ;
    coordT1 = glT1-1 ;
    coordP1 = m-1 ;
    other3 = glT1-m+2 ;
    other4 = m-1 ;
    other5 = glT1-m+2 ;
    S6(glT1-1,m-1,glT1,m,glT1-m+2,m-1,glT1-m+2,m) ;
    other1 = glT1+1 ;
    other3 = glT1-m+2 ;
    other4 = m-1 ;
    other5 = glT1-m+3 ;
    S7(glT1-1,m-1,glT1+1,m,glT1-m+2,m-1,glT1-m+3,m) ;
    other1 = glT1-m+2 ;
    S1(glT1,m,glT1-m+2,m) ;
    local1 = glT1-m+2 ;
    other1 = glT1+1 ;
    other3 = glT1-m+2 ;
    other5 = glT1-m+3 ;
    S2(glT1,m,glT1+1,m,glT1-m+2,m,glT1-m+3,m) ;
    for (coordP1=glT1-m+2;coordP1<=m;coordP1++) {
      S8(glT1,coordP1) ;
    }
  }
  if (m >= 3) {
    glT1 = 2*m-3 ;
    rp1 = m-1 ;
    local1 = m-1 ;
    local2 = m-1 ;
    coordT1 = 2*m-4 ;
    coordP1 = m-1 ;
    other1 = 2*m-3 ;
    other2 = m-1 ;
    other3 = m-1 ;
    other4 = m-1 ;
    other6 = m-1 ;
    S5(2*m-4,m-1,2*m-3,m-1,m-1,m-1,m,m-1) ;
    local2 = m-2 ;
    coordT1 = 2*m-4 ;
    coordP1 = m-2 ;
    other1 = 2*m-3 ;
    other2 = m-1 ;
    other4 = m-2 ;
    other6 = m-1 ;
    S6(2*m-4,m-2,2*m-3,m-1,m,m-2,m,m-1) ;
    coordT1 = 2*m-3 ;
    coordP1 = m-1 ;
    other2 = m-1 ;
    S1(2*m-3,m-1,m,m-1) ;
    local2 = m-1 ;
    coordT1 = 2*m-3 ;
    coordP1 = m-1 ;
    other1 = 2*m-2 ;
    other4 = m-1 ;
    S3(2*m-3,m-1,2*m-2,m,m,m-1,m,m) ;
    local1 = m-2 ;
    coordT1 = 2*m-4 ;
    other1 = 2*m-3 ;
    other3 = m-2 ;
    other5 = m-1 ;
    S5(2*m-4,m,2*m-3,m,m-2,m,m-1,m) ;
    local1 = m-1 ;
    local2 = m-1 ;
    coordT1 = 2*m-4 ;
    coordP1 = m-1 ;
    other1 = 2*m-3 ;
    other3 = m-1 ;
    other4 = m-1 ;
    other5 = m-1 ;
    S6(2*m-4,m-1,2*m-3,m,m-1,m-1,m-1,m) ;
    other1 = 2*m-2 ;
    other3 = m-1 ;
    other4 = m-1 ;
    S7(2*m-4,m-1,2*m-2,m,m-1,m-1,m,m) ;
    coordT1 = 2*m-3 ;
    other1 = m-1 ;
    S1(2*m-3,m,m-1,m) ;
    local1 = m-1 ;
    coordT1 = 2*m-3 ;
    other1 = 2*m-2 ;
    other3 = m-1 ;
    S2(2*m-3,m,2*m-2,m,m-1,m,m,m) ;
    coordT1 = 2*m-3 ;
    for (coordP1=m-1;coordP1<=m;coordP1++) {
      S8(2*m-3,coordP1) ;
    }
  }
  if (m == 2) {
    S5(0,1,1,1,1,1,2,1) ;
    S1(1,1,2,1) ;
    S3(1,1,2,2,2,1,2,2) ;
    S6(0,1,1,2,1,1,1,2) ;
    S7(0,1,2,2,1,1,2,2) ;
    S1(1,2,1,2) ;
    S2(1,2,2,2,1,2,2,2) ;
    for (coordP1=1;coordP1<=2;coordP1++) {
      S8(1,coordP1) ;
    }
  }
  if (m >= 2) {
    glT1 = 2*m-2 ;
    local1 = m-1 ;
    coordT1 = 2*m-3 ;
    other1 = 2*m-2 ;
    other3 = m-1 ;
    S5(2*m-3,m,2*m-2,m,m-1,m,m,m) ;
    local2 = m-1 ;
    coordT1 = 2*m-3 ;
    coordP1 = m-1 ;
    other1 = 2*m-2 ;
    other4 = m-1 ;
    S6(2*m-3,m-1,2*m-2,m,m,m-1,m,m) ;
    coordT1 = 2*m-2 ;
    S1(2*m-2,m,m,m) ;
    coordT1 = 2*m-2 ;
    S8(2*m-2,m) ;
  }
}
