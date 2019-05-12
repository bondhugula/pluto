#define length 100

// dynamic programming
int main() {
  int i, j, k;
  int W[length][length]; //:input;
  int sum_c[length][length][length];
  int c[length][length];
  int out[1]; //:output;

  for (i = 0; i <= length - 2; i++) {
    for (j = length - 1 - i; j <= length - 1; j++) {
      sum_c[length - 2 - i][j][length - 2 - i] = 0;
      for (k = length - 2 - i + 1; k <= j - 1; k++) {
        sum_c[length - 2 - i][j][k] =
            sum_c[length - 2 - i][j][k - 1] + c[length - 2 - i][k] + c[k][j];
      }
      c[length - 2 - i][j] =
          sum_c[length - 2 - i][j][j - 1] + W[length - 2 - i][j];
    }
  }
  out[0] = c[0][length - 1];

  return 0;
}
