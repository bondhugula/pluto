
#pragma scop

for (x = 0; x < 10000; x = x + 1) {
  for (y = 0; y < 10000; y = y + 1) {
    green[x][y] = x + y;
  }
}

for (x = 0; x < 10000; x = x + 1) {
  for (y = 0; y < 10000; y = y + 1) {
    red[x][y] = x + y;
  }
}
for (x = 0; x < 10000; x = x + 1) {
  for (y = 0; y < 10000; y = y + 1) {
    red[x][y] = diff(red[x][y], green[x][y]);
  }
}

for (x = 1; x < 4999; x += 1) {
  for (y = 1; y < 4999; y += 1) {
    red[2 * x + 1][2 * y + 1] =
        RBK_3x3_1(red[2 * x][2 * y], red[2 * x][2 * y + 2],
                  red[2 * x + 2][2 * y], red[2 * x + 2][2 * y + 2]);
  }
}

#pragma endscop
