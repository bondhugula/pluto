
/*@ begin Loop(
  transform ScalarReplace(prefix='it')
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
                        { 
                          C[i][j]=A[i][k]*B[k][j]+C[i][j]; 
                          C[i][j+1]=A[i][k]*B[k][j+1]+C[i][j+1]; 
                          C[i][j+2]=A[i][k]*B[k][j+2]+C[i][j+2]; 
                          C[i][j+3]=A[i][k]*B[k][j+3]+C[i][j+3]; 
                          C[i+1][j]=A[i+1][k]*B[k][j]+C[i+1][j]; 
                          C[i+1][j+1]=A[i+1][k]*B[k][j+1]+C[i+1][j+1]; 
                          C[i+1][j+2]=A[i+1][k]*B[k][j+2]+C[i+1][j+2]; 
                          C[i+1][j+3]=A[i+1][k]*B[k][j+3]+C[i+1][j+3]; 
                          C[i+2][j]=A[i+2][k]*B[k][j]+C[i+2][j]; 
                          C[i+2][j+1]=A[i+2][k]*B[k][j+1]+C[i+2][j+1]; 
                          C[i+2][j+2]=A[i+2][k]*B[k][j+2]+C[i+2][j+2]; 
                          C[i+2][j+3]=A[i+2][k]*B[k][j+3]+C[i+2][j+3]; 
                          C[i+3][j]=A[i+3][k]*B[k][j]+C[i+3][j]; 
                          C[i+3][j+1]=A[i+3][k]*B[k][j+1]+C[i+3][j+1]; 
                          C[i+3][j+2]=A[i+3][k]*B[k][j+2]+C[i+3][j+2]; 
                          C[i+3][j+3]=A[i+3][k]*B[k][j+3]+C[i+3][j+3]; 
                        } 
                    for (; j<=32*c5+31; j=j+1) 
                      for (k=128*c6; k<=128*c6+127; k++) 
                        { 
                          C[i][j]=A[i][k]*B[k][j]+C[i][j]; 
                          C[i+1][j]=A[i+1][k]*B[k][j]+C[i+1][j]; 
                          C[i+2][j]=A[i+2][k]*B[k][j]+C[i+2][j]; 
                          C[i+3][j]=A[i+3][k]*B[k][j]+C[i+3][j]; 
                        } 
                  } 
                for (; i<=32*c4+31; i=i+1) 
                  { 
                    for (j=32*c5; j<=32*c5+31-3; j=j+4) 
                      for (k=128*c6; k<=128*c6+127; k++) 
                        { 
                          C[i][j]=A[i][k]*B[k][j]+C[i][j]; 
                          C[i][j+1]=A[i][k]*B[k][j+1]+C[i][j+1]; 
                          C[i][j+2]=A[i][k]*B[k][j+2]+C[i][j+2]; 
                          C[i][j+3]=A[i][k]*B[k][j+3]+C[i][j+3]; 
                        }
                    for (; j<=32*c5+31; j=j+1) 
                      for (k=128*c6; k<=128*c6+127; k++) 
                        C[i][j]=A[i][k]*B[k][j]+C[i][j]; 
                  } 
              } 
) @*/

/*@ end @*/

