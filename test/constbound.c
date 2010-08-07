
#pragma scop
for (i = 0 ; i < 4 ; i++)
for (j = 0 ; j < 4 ; j++)
a[i] = b[i][j] + a[i];
#pragma endscop
