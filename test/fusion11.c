// CHECK: Output written
#pragma scop
for (int _i0 = 0; (_i0 < (N - 7371)); _i0++) {
  for (int _i1 = 0; (_i1 <= 7238); _i1++) {
    yf1[_i0] = (yf1[_i0] + (yds[((7238 + _i0) - _i1)] * taps1[_i1]));
  }
}
for (int _i0 = 7238; (_i0 < (N - 14609 + 7238)); _i0++) {
  for (int _i1 = 0; (_i1 <= 7238); _i1++) {
    yf2[_i0 - 7238] = (yf2[_i0 - 7238] + (yf1[((_i0)-_i1)] * taps2[_i1]));
  }
}
#pragma endscop
