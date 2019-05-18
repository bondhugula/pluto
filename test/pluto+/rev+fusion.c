// CHECK: T(S1): (i, 0)
// CHECK: T(S2): (N-i-1, 1)
int main() {
  long N = 1000;
  int a[N];
  int b[N];
  int c[N];
  int i, j;

#pragma scop
  for (i = 0; i < N; i++) {
    a[i] = b[i - 1];
  }
  for (i = 0; i < N; i++) {
    c[i] = a[N - 1 - i];
  }
#pragma endscop
}

