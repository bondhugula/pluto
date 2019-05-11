
#pragma scop
QR[k][k] = QR[k][k] + 1;
for (j = k + 1; j < N; j++) {
  s[k][j] = 0.0;
  for (i = k; i < N; i++) {
    QR[i][j] = QR[i][j] + s[k][j] * QR[i][k];
  }
}
#pragma endscop
