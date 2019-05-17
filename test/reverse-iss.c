int main() {
  int i;
  int a[1000];

  long N = 1000;

#pragma scop
  for (i = 0; i < N; i++) {
    a[i] = 2 * a[N - 1 - i];
  }
#pragma endscop

  return 0;
}
