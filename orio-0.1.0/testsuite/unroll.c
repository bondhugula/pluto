
/*@ begin Loop (
for (i = 0; i <= M-1; i++)
  transform Unroll(ufactor=2)
  for (j = 0; j <= N-1; j++)
    S(i,j,k);
) @*/

for (i = 0; i <= M-1; i++)
  for (j = 0; j <= N-1; j++)
    S(i,j,k);

/*@ end @*/
