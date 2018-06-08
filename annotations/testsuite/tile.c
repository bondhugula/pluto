
/*@ begin Loop(
 transform Tile(tsize=32, tindex=it)
 for (i = 1; i <= m; i++)
   transform Tile(tsize=32, tindex=jt)
   for (j = 1; j <= n; j++)
     S(i,j);
) @*/
for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
        S(i,j);
/*@ end @*/


/*@ begin Loop(
 transform Permut(order=(it,jt,i,j))
 transform Tile(tsize=16, tindex=it)
 for (i = 1; i <= m; i++)
   transform Tile(tsize=16, tindex=jt)
   for (j = 1; j <= n; j++)
     S(i,j);
) @*/
for (i = 1; i <= m; i++)
    for (j = 1; j <= n; j++)
        S(i,j);
/*@ end @*/
