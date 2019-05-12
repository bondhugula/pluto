
/* pluto start (JS,JE,LS,LE) */
DO J = JS, JE DO M = 1, 5 DO N = 1, 5 DO L = LS,
   LE B(M, N, J, L) = B(M, N, J, L) - A(M, 1, J, L) * B(1, N, J - 1, L) -
                      A(M, 2, J, L) * B(2, N, J - 1, L) -
                      A(M, 3, J, L) * B(3, N, J - 1, L) -
                      A(M, 4, J, L) * B(4, N, J - 1, L) -
                      A(M, 5, J, L) *
                          B(5, N, J - 1, L) END DO END DO END DO END DO
                          /* pluto end */
