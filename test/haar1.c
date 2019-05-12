constant n, m;

for (i = 0; i < n; i++)
  vecp[i] = 0;

w = n;

for (j = 0; j < m; j++) {
  for (i = 0; i < n / 2; i++) {
    vecp[i] = (vec[2 * i] + vec[2 * i + 1]) / (2.0);
  }
}
