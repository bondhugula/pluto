/* pluto start (n) */
#pragma scop
for (i = 1; i < n - 1; i++) {
  a1[i] = a0[i - 1] + a0[i] + a0[i + 1];
}
for (i = 2; i < n - 2; i++) {
  a2[i] = a1[i - 1] + a1[i] + a1[i + 1];
}
for (i = 3; i < n - 3; i++) {
  a3[i] = a2[i - 1] + a2[i] + a2[i + 1];
}
for (i = 4; i < n - 4; i++) {
  a4[i] = a3[i - 1] + a3[i] + a3[i + 1];
}
for (i = 5; i < n - 5; i++) {
  a5[i] = a4[i - 1] + a4[i] + a4[i + 1];
}
#pragma endscop
/* pluto end */
