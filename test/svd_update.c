#define n 100

// SVD Update algorithm
void main() {
  int i, j, l;
  // int n=100;//n is the parameter to be modified
  int R[n][n][2 * n - 1];
  int V[n][n][n];
  int A[n][2 * n];
  int a[n]; //:input;
  int alpha[n];
  int theta[n];
  int phi[n];
  int Rout[n][n]; //:output;
  int Vout[n][n]; //:output;
  for (j = 0; j <= n - 1; j++) {
    A[j][0] = 0;
    for (i = 0; i <= n - 1; i++) {
      A[j][i + 1] = A[j][i] + a[i] * V[i][j][n - 1];
    }
  }

  for (i = 0; i <= n - 1; i++) {
    alpha[i] = A[i][i + n] - R[i][i][2 * n - 2];
    for (j = 0; j <= n - 1; j++) {
      if (j < i) {
        R[i][j][0] = 0;
        A[j][i + n + 1] = A[j][i + n];
      }
      if (j >= i) {
        R[i][j][0] = R[i][j][2 * n - 2] * alpha[i] + A[j][i + n] * alpha[i];
        A[j][i + n + 1] =
            R[i][j][2 * n - 2] * alpha[i] + A[j][i + n] * alpha[i];
      }
      V[i][j][0] = V[i][j][n - 1];
    }
  }

  for (l = 0; l <= n - 2; l++) {

    theta[l] = R[l + 1][l + 1][2 * l] * R[l][l + 1][2 * l] +
               R[l][l][2 * l] * R[l][l][2 * l] -
               R[l + 1][l + 1][2 * l] * R[l + 1][l + 1][2 * l] +
               R[l][l + 1][2 * l] * R[l][l + 1][2 * l];

    phi[l] =
        R[l + 1][l + 1][2 * l] * theta[l] + R[l][l + 1][2 * l] - R[l][l][2 * l];
    for (i = 0; i <= n - 1; i++) {
      for (j = 0; j <= n - 1; j++) {
        int tmp1 = R[i + 1][j][2 * l];
        int tmp2 = R[i][j][2 * l];
        int tmp3 = R[i - 1][j][2 * l];
        int res;

        if (j >= l) {
          if (i == l) {
            res = tmp1 * theta[l] - tmp2 * theta[l];
          }
          if (i == l + 1) {
            res = tmp2 * theta[l] + tmp3 * theta[l];
          }

          if (i < l) {
            res = tmp2;
          }
          if (i > l + 1) {
            res = tmp2;
          }
        }
        if (j < l) {
          res = tmp2;
        }
        R[i][j][2 * l + 1] = res;
      }
    }

    for (i = 0; i <= n - 1; i++) {
      for (j = 0; j <= n - 1; j++) {
        int tmp1 = R[i][j + 1][2 * l + 1];
        int tmp2 = R[i][j][2 * l + 1];
        int tmp3 = R[i][j - 1][2 * l + 1];
        int res;
        if (i <= l + 1) {
          if (j == l) {
            res = tmp1 * phi[l] - tmp2 * phi[l];
          }
          if (j == l + 1) {
            res = tmp2 * phi[l] + tmp3 * phi[l];
          }
          if (j < l) {
            res = tmp2;
          }
          if (j > l + 1) {
            res = tmp2;
          }
        }
        if (i > l + 1) {
          res = tmp2;
        }

        R[i][j][2 * l + 2] = res;

        tmp1 = V[i][j + 1][l];
        tmp2 = V[i][j][l];
        tmp3 = V[i][j - 1][l];

        if (j == l) {
          res = tmp1 * phi[l] - tmp2 * phi[l];
        }
        if (j == l + 1) {
          res = tmp2 * phi[l] + tmp3 * phi[l];
        }
        if (j < l) {
          res = tmp2;
        }
        if (j > l + 1) {
          res = tmp2;
        }
        V[i][j][l + 1] = res;
      }
    }
  }

  for (i = 0; i <= n - 1; i++) {
    for (j = 0; j <= n - 1; j++) {
      Rout[i][j] = R[i][j][2 * n - 2];
      Vout[i][j] = V[i][j][n - 1];
    }
  }
}
