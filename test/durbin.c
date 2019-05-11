#define N 100

// durbin
int main(int argc, char *argv[]) {
  int k, i;
  int y[N][N];
  int sum[N][N];
  int beta[N];
  int alpha[N];
  int r[N];   // input
  int out[N]; // output

  y[0][0] = r[0];
  beta[0] = 1;
  alpha[0] = r[0];

  for (k = 1; k <= N - 1; k++) {
    beta[k] = beta[k - 1] - alpha[k - 1] * alpha[k - 1] * beta[k - 1];
    sum[0][k] = r[k];
    for (i = 0; i <= k - 1; i++) {
      sum[i + 1][k] = sum[i][k] + r[k - i - 1] * y[i][k - 1];
    }
    alpha[k] = -sum[k][k] * beta[k];
    for (i = 0; i <= k - 1; i++) {
      y[i][k] = y[i][k - 1] + alpha[k] * y[k - i - 1][k - 1];
    }
    y[k][k] = alpha[k];
  }
  for (i = 0; i <= N - 1; i++) {
    out[i] = y[i][N - 1];
  }
}
