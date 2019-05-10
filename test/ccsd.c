#pragma scop
for (a = 0; a < V; a++)
  for (b = 0; b < V; b++)
    for (c = 0; c < V; c++)
      for (e = 0; e < V; e++)
        for (i = 0; i < O; i++)
          for (j = 0; j < O; j++)
            for (k = 0; k < O; k++)
              for (m = 0; m < O; m++)
                X[a][b][c][i][j][k] = X[a][b][c][i][j][k] +
                                      T2[a][b][k][m] * O1[c][m][i][j] +
                                      T2[c][e][i][j] * O2[a][b][e][k];

for (a = 0; a < V; a++)
  for (b = 0; b < V; b++)
    for (c = 0; c < V; c++)
      for (i = 0; i < O; i++)
        for (j = 0; j < O; j++)
          for (k = 0; k < O; k++)
            Y[a][b][c][i][j][k] = T1[c][k] * O3[a][b][i][j];

for (a = 0; a < V; a++)
  for (b = 0; b < V; b++)
    for (c = 0; c < V; c++)
      for (i = 0; i < O; i++)
        for (j = 0; j < O; j++)
          for (k = 0; k < O; k++)
            e1 = e1 + X[a][b][c][i][j][k] * X[a][b][c][i][j][k];

for (a = 0; a < V; a++)
  for (b = 0; b < V; b++)
    for (c = 0; c < V; c++)
      for (i = 0; i < O; i++)
        for (j = 0; j < O; j++)
          for (k = 0; k < O; k++)
            e2 = e2 + X[a][b][c][i][j][k] * Y[a][b][c][i][j][k];
#pragma endscop
