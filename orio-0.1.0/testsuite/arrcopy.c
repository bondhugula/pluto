
/*@ begin Loop(
  transform ArrCopy(aref='C[i][j]', dimsizes=[32,32], suffix='_copy')
  transform ArrCopy(aref='A[i][k]', dimsizes=[32,32])
  transform ArrCopy(aref='B[k][j]', dimsizes=[32,32], dtype='double')
  for (iii=0; iii<=cbv_1; iii=iii+128)
   for (jjj=0; jjj<=N-1; jjj=jjj+512)
    for (kkk=0; kkk<=K-1; kkk=kkk+64)
     for (ii=iii; ii<=min(M-1,iii+96); ii=ii+32)
      for (jj=jjj; jj<=min(N-1,jjj+480); jj=jj+32)
       for (kk=kkk; kk<=min(K-1,kkk+32); kk=kk+32)
        for (i=ii; i<=min(M-1,ii+31); i=i+1)
	 for (j=jj; j<=min(N-1,jj+31); j=j+1)
          for (k=kk; k<=min(K-1,kk+31); k=k+1)
           C[i][j]=C[i][j]+A[i][k]*B[k][j];
) @*/

/*@ end @*/


