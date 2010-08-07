 /*@ begin Loop( 
 transform UnrollJam(ufactor=4)
 for (i = 1; i <= m; i++)
   transform Unroll(ufactor=2)
   for (j = 1; j <= n; j++)
     {
       S1(i,j);
       S2(i,j);
     }
 ) @*/
 for (i = 1; i <= m; i++)
   for (j = 1; j <= n; j++)
     {
       S1(i,j);
       S2(i,j);
     }
 /*@ end @*/



