
/*@ begin Loop (
transform UnrollJam(ufactor=2)
for (i = 0; i <= M-1; i++)
  transform UnrollJam(ufactor=2)
  for (j = 0; j <= N-1; j++)
    transform UnrollJam(ufactor=2)
    for (k = 0; k <= O-1; k++)
      {
        S(i,j,k);
      }
) @*/

for (i = 0; i <= M-1; i++)
    for (j = 0; j <= N-1; j++)
        for (k = 0; k <= O-1; k++) {
            S(i,j,k);
        }

/*@ end @*/
