constant n, m;

for (y_p = 1; y_p <= m - 1; y_p++) {
  sad[y_p][0] = sad[y_p - 1][m - 1] + curr[y_p][0] + prev[y_p][0];

  for (x_p = 1; x_p <= m - 1; x_p++) {
    sad[y_p][x_p] = sad[y_p][x_p - 1] + curr[y_p][x_p] + prev[y_p][x_p];
  }
}

result = sad[m - 1][m - 1];
