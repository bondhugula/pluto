/* pluto start (n,m) */

for (i = 0; i < n; i++) {
  vec[0][i] = 0;
}

for (i = 1; i < m; i++) {
  for (j = 0; j < n; j++) {
    vec[i, j] = (vec[i - 1, 2 * j] + vec[i - 1, 2 * j + 1]);
  }
  for (j = 0; j < n; j++) {
    vec[i, j + n] = (vec[i - 1, 2 * j] - vec[i - 1, 2 * j + 1]);
  }
}
/* pluto end */
