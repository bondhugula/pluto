constant N;

DO a = 1, N DO q = 1, N DO r = 1, N DO s = 1, N DO p = 1,
   N T1[a, q, r, s] =
       T1[a, q, r, s] + A[p, q, r, s] *
                            C4[p, a] END DO END DO END DO END DO END DO

                                DO a = 1,
   N DO b = 1, N DO r = 1, N DO s = 1, N DO q = 1,
   N T2[a, b, r, s] =
       T2[a, b, r, s] +
       T1[a, q, r, s] * C3[q, b] END DO END DO END DO END DO END DO
