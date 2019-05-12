#pragma scop
for (i0 = 0; i0 < 100; i0++) {
  for (i1 = 0; i1 < 100; i1++) {
    for (i2 = 0; i2 < 100; i2++) {
      for (i3 = 0; i3 < 100; i3++) {
        for (i4 = 0; i4 < 100; i4++) {
          for (i5 = 0; i5 < 100; i5++) {
            for (i6 = 0; i6 < 100; i6++) {
              for (i7 = 0; i7 < 100; i7++) {
                for (i8 = 0; i8 < 100; i8++) {
                  for (i9 = 0; i9 < 100; i9++) {
                    s += 1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
#pragma endscop
