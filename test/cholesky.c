#pragma scop
for (i = 1; i < nSlip; i++) {
  fdot[i] = 0.0;
  for (k = 0; k < i; k++)
    fdot[i] += a[i][k] * a[k][i];
  a[i][i] = a[i][i] - fdot[i];
  for (j = i + 1; j < nSlip; j++) {
    fdot[i] = 0.0;
    for (k = 0; k < i; k++)
      fdot[i] += a[i][k] * a[k][j];
    a[i][j] = a[i][j] - fdot[i];
    fdot[i] = 0.0;
    for (k = 0; k < i; k++)
      fdot[i] += a[j][k] * a[k][i];
    a[j][i] = (a[j][i] - fdot[i]) / a[i][i];
  }
}
#pragma endscop
