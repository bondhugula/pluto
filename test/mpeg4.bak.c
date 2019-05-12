constant n, m;

for (y_s = 0; y_s <= n - 1; y_s++) {
  for (x_s = 0; x_s <= n - 1; x_s++) {
    sad[y_s][x_s][0][0] = curr[0][0] + prev[y_s][x_s];
    for (x_p = 1; x_p <= m - 1; x_p++) {
      sad[y_s][x_s][0][x_p] =
          sad[y_s][x_s][0][x_p - 1] + curr[0][x_p] + prev[y_s][x_s + x_p];
    }
    for (y_p = 1; y_p <= m - 1; y_p++) {
      sad[y_s][x_s][y_p][0] =
          sad[y_s][x_s][y_p - 1][m - 1] + curr[y_p][0] + prev[y_s + y_p][x_s];

      for (x_p = 1; x_p <= m - 1; x_p++) {
        sad[y_s][x_s][y_p][x_p] = sad[y_s][x_s][y_p][x_p - 1] + curr[y_p][x_p] +
                                  prev[y_s + y_p][x_s + x_p];
      }
    }
  }
}

for (y_s = 0; y_s <= n - 1; y_s++) {
  for (x_s = 0; x_s <= n - 1; x_s++) {
    result[y_s][x_s] = sad[y_s][x_s][m - 1][m - 1];
  }
}
