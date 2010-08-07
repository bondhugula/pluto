
/*@ begin Loop (
 transform UnrollJam(ufactor=4)
 for (c3=max(32*c1,0);c3<=min(32*c1+31,N-1);c3++) 
   {
     transform Unroll(ufactor=2)
     for (c4=max(0,32*c2);c4<=min(32*c2+31,N-1);c4++) 
       {
        S1(c1,c2,c3,c4);
        S2(c2,c1,c4,c3);
       }
   }
) @*/
for (c3=max(32*c1,0);c3<=min(32*c1+31,N-1);c3++) 
  {
    for (c4=max(0,32*c2);c4<=min(32*c2+31,N-1);c4++) 
      {
	S1(c1,c2,c3,c4);
	S2(c2,c1,c4,c3);
      }
  }
/*@ end @*/

