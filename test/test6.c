constant M, N;

for (i = 0; i < M; i++) {
  a[i] = i;
  for (j = 0; j < N; j++) {
    b[j] = (b[j] - a[i]) / 2;
  }
}
