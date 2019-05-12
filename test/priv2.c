

float a[1000];
float b[1000][1000];

main(int argc, char *argv[]) {
  int i, j, k;
  long N;

  float dist, dist1;

  // for (i=0; i<100; i++) {
  // a[i] = i;
  // }

  N = atoi(argv[0]);
  N = 100;

#pragma scop
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      dist = 0.0;
      dist1 = 0.0;
      for (k = 0; k < N; k++) {
        dist += a[10 * i + j] * a[10 * i + j];
        dist1 += a[10 * i + j] + a[10 * i + j];
      }
      b[i][j] = dist;
      b[i][j] += dist1;
    }
  }
#pragma endscop

  return (int)round(b[0][50]);
}
