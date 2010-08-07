
/*@ begin Loop( 
 transform RegTile(loops=['i','j','k'], ufactors=[2,2,2])
 for (i=0; i<=M-1; i++)
 {
   for (j=0; j<=N-1; j++)
   {
     for (k=i; k<=O-1; k++)
     {
       S(i,j,k);
     }
   }
 }
) @*/ 

/*@ end @*/


