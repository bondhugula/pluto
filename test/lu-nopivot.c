CONSTANT n;

DO k = 1, n - 1 DO i = k + 1, n a[k, i] = a[k, i] / a[k, k];
END DO DO i = k + 1, n DO j = k + 1,
          n a[i, j] = a[i, j] - a[i, k] * a[k, j] END DO END DO END DO
