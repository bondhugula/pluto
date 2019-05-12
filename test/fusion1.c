// CHECK: Output written
#pragma scop
for (i = 1; i < N - 2; i++) {
  A[i] = 0.33 * (In[i - 1] + In[i] + In[i + 1]);
}

for (i = 2; i < N - 3; i++) {
  Out[i] = 0.33 * (A[i - 1] + A[i] + A[i + 1]);
}
#pragma endscop
