#pragma scop
a ^= (b | (c + 2)) % 2 - !(a[n]/2);
#pragma endscop
