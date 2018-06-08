
/*@ begin Loop(
  transform BoundReplace()
  for (c1=0; c1<=(N-1)/512; c1++)
    for (c2=0; c2<=(N-1)/512; c2++)
      for (c3=0; c3<=(N-1)/1024; c3++)
        for (c4=16*c1; c4<=16*c1+15; c4++)
          for (c5=16*c2; c5<=16*c2+15; c5++)
            for (c6=8*c3; c6<=8*c3+7; c6++)
              {
                for (i=32*c4; i<=32*c4+31-3; i=i+4)
                  {
                    for (j=32*c5; j<=32*c5+31-3; j=j+4)
                      for (k=128*c6; k<=128*c6+127; k++)
                        S(i,j,k);
                    for (; j<=32*c5+31; j=j+1)
                      for (k=128*c6; k<=128*c6+127; k++)
                        S(i,j,k);
                  }
                for (; i<=32*c4+31; i=i+1)
                  {
                    for (j=32*c5; j<=32*c5+31-3; j=j+4)
                      for (k=128*c6; k<=128*c6+127; k++)
                        S(i,j,k);
                    for (; j<=32*c5+31; j=j+1)
                      for (k=128*c6; k<=128*c6+127; k++)
                        S(i,j,k);
                  }
              }
) @*/

/*@ end @*/

