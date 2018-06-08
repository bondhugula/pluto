
/*@ begin Loop (
for (it = 0; it <= m - 1; it = it + 32)
 for (jt = 0; jt <= n - 1; jt = jt + 32)
  for (i = it; i <= min(m - 1, it + 31); i = i + 1)
   transform Unroll(ufactor=4)
   for (j = jt; j <= min(n - 1, jt + 31); j = j + 1)
    S(i, j);
) @*/

for (it = 0; it <= m - 1; it = it + 32)
    for (jt = 0; jt <= n - 1; jt = jt + 32)
        for (i = it; i <= min(m - 1, it + 31); i = i + 1)
            for (j = jt; j <= min(n - 1, jt + 31); j = j + 1)
                S(i, j);

/*@ end @*/
