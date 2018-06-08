
/*@ begin Loop(
 transform Permut(order=(k,(x),j,i))
   for (i = 1; i <= p; i += 1)
    for (j = 1; j <= q; j += 1)
     for (k = 1; k <= r; k += 1)
      S(i,j,k);
) @*/
for (i = 1; i <= p; i += 1)
    for (j = 1; j <= q; j += 1)
        for (k = 1; k <= r; k += 1)
            S(i,j,k);
/*@ end @*/

