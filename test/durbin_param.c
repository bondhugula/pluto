int N;

// durbin
void main() {
  int k, i;
  //  int N=100;//this parameter is to be changed
  int **y;
  int **sum;
  int *beta;
  int *alpha;
  int *r;
  int *out;

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
