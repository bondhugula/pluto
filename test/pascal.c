CONSTANT n, m;

for (i = 0; i <= n; i++) {
  PD[i, 0] = 1 PD[i, i] = 1 for (j = 1; j < i; j++) {
    PD[i, j] = PD[i - 1, j - 1] + PD[i - 1, j]
  }
}

result = PD[n, m]
