
/*@ begin Loop(
 transform Tile(tsize=32, tindex='ii')
 for (i = 0; i <= m-1; i++)
   transform Tile(tsize=32, tindex='jj')
   for (j = 0; j <= n-1; j++)
     S(i,j);
) @*/

for (i = 0; i <= m-1; i++)
    for (j = 0; j <= n-1; j++)
        S(i,j);

/*@ end @*/


