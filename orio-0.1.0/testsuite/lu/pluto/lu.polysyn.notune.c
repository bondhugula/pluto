

register int i,j,k;
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6, 
  newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;    
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6, 
  newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;    


/*@ begin PolySyn( 
 l1_tiles = [32,16,16]; 
 l2_tiles = [8,1,4]; 
 hotspot_permut = [0,1,2]; 
 unroll_factors = [4,4,1]; 
 parallelize = True; 
 scalar_replace = False; 
 icc_vectorize = True; 
 ) @*/  

/* pluto start (N) */ 
for (k=0; k<=N-1; k++) 
  {
    for (j=k+1; j<=N-1; j++)
      A[k][j] = A[k][j]/A[k][k];
  /*@ begin Loop(
    transform Composite(permut = [('i','j')])
    for(i=k+1; i<=N-1; i++)
      for (j=k+1; j<=N-1; j++)   
	A[i][j] = A[i][j] - A[i][k]*A[k][j];
  ) @*/
  /*@ end @*/
  }
/* pluto end */

/*@ end @*/


