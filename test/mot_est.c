// motion estimation
void main() {
  int x_s, y_s, x_p, y_p;
  int sad[32][32][16][16];
  int curr[16][16];   //:input;
  int prev[46][46];   //:input;
  int result[32][32]; //:output;

  for (y_s = 0; y_s <= 31; y_s++) {
    for (x_s = 0; x_s <= 31; x_s++) {
      for (y_p = 0; y_p <= 15; y_p++) {
        for (x_p = 0; x_p <= 15; x_p++) {
          if ((x_p == 0) && (y_p == 0)) {
            sad[y_s][x_s][y_p][x_p] =
                curr[y_p][x_p] + prev[y_s + y_p][x_s + x_p];
          }
          if ((x_p == 0) && (y_p > 0)) {
            sad[y_s][x_s][y_p][x_p] = sad[y_s][x_s][y_p - 1][15] +
                                      curr[y_p][x_p] +
                                      prev[y_s + y_p][x_s + x_p];
          }
          if (x_p > 0) {
            sad[y_s][x_s][y_p][x_p] = sad[y_s][x_s][y_p][x_p - 1] +
                                      curr[y_p][x_p] +
                                      prev[y_s + y_p][x_s + x_p];
          }
        }
      }
    }
  }
  for (y_s = 0; y_s <= 31; y_s++) {
    for (x_s = 0; x_s <= 31; x_s++) {
      result[y_s][x_s] = sad[y_s][x_s][15][15];
    }
  }
}
