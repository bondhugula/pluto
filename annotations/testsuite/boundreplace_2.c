
/*@ begin Loop(
  transform BoundReplace(dtype="int")
  for (c7=0; c7<=X1; c7++)
    for (c9=0; c9<=X2; c9++)
      for (c8=max(c7+1,300*c5); c8<=min(N-1,300*c5+299); c8++)
        S2(c1,c3,c2,c4,c6,c5,c7,c9,c8) ;
) @*/

/*@ end @*/

