for (t=1; t<=T; t++)
  {
    for (i=1; i<=N-1; i++)
      e[i] = e[i] - coeff1*(h[i]-h[i-1]);
    for (i=0; i<=N-1; i++)
      h[i] = h[i] - coeff2*(e[i+1]-e[i]);
  }
