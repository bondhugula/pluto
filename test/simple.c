
#pragma scop

for (i=0; i<N; i++) {
    a[i] =   b[i] + c[i];
}

#pragma endscop
