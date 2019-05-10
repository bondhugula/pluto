constant n, m;

for (y_s = 0; y_s < n; y_s++) {
  for (x_s = 0; x_s < n; x_s++) {
    sad[y_s][x_s] = 0 for (y_p = 0; y_p < m; y_p++) {
      for (x_p = 0; x_p < m; x_p++) {
        sad[y_s][x_s] =
            sad[y_s][x_s] + curr[y_p][x_p] - prev[y_s + y_p][x_s + x_p];
      }
    }
  }
}

for (y_s = 0; y_s < n; y_s++) {
  for (x_s = 0; x_s < n; x_s++) {
    result[y_s][x_s] = sad[y_s][x_s];
  }
}
