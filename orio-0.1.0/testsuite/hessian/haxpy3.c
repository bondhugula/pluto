for (i=0; i<=N-1; i=i+1)
  for (j=0; j<=N-1; j=j+1) 
    Y[i][j]=
      a0*X0[i][j]+a1*X1[i][j]+a2*X2[i][j]
      +2.0*b00*u0[i]*u0[j]
      +2.0*b11*u1[i]*u1[j]
      +2.0*b22*u2[i]*u2[j]
      +b01*(u0[i]*u1[j]+u1[i]*u0[j])
      +b02*(u0[i]*u2[j]+u2[i]*u0[j])
      +b12*(u1[i]*u2[j]+u2[i]*u1[j]);
