

/*@ begin Loop (
transform Composite(
 tile = [('i',32,'ii'),('j',64,'jj'),('k',128,'kk')],
 permut = [(['ii'],['jj'],['kk'],'i','j','k')],
 arrcopy = [('C[i][j]',[32,64]),('A[i][k]',[32,128]),('B[k][j]',[128,64])],
 unrolljam = [('i',2),('j',2),('k',2)],
 scalarreplace = (True, 'double', 'itv_'),
 boundreplace = (True, 'lbv_', 'ubv_'),
 pragma = [(['dummy'], ['ivdep', 'vector always'])],
 openmp = (True, 'omp parallel for private(ii,jj,kk,i,j,k,A_buf,B_buf,C_buf)'),
 vector = (True, ['ivdep','vector always'])
)
for (i = 0; i <= M-1; i++)
  for (j = 0; j <= N-1; j++)
    for (k = 0; k <= O-1; k++)
      {
        C[i][j] = C[i][j] + A[i][k]*B[k][j];
      }
) @*/

for (i = 0; i <= M-1; i++)
    for (j = 0; j <= N-1; j++)
        for (k = 0; k <= O-1; k++) {
            C[i][j] = C[i][j] + A[i][k]*B[k][j];
        }

/*@ end @*/



