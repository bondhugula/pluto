int n;

void main(void)
{
  int i,j,k,l;
  double **a,s;
  
#pragma scop
  for(k = 1; k <= n; k++)
    {
      s = 1/a[k][k];

      for(l = k+1; l <= n; l++)
	a[l][k] = a[l][k]*s;

      for(j = k+1; j <= n+1; j++)
	for(i = k+1; i <= n; i++)
	  a[i][j] = a[i][j] - a[k][j] * a[i][k];
    }
#pragma endscop
}
