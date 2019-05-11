#pragma scop

for (i = 0; i < 10000; i++) {
  for (j = 0; j <= i; j++) {
    a[i] += b[i] * c[j];
  }
}

#pragma endscop
