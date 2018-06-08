
/* pluto start (n) */
do i = 1, n
           do j = 2, n
                      b[i,j] = a[i,j-1] + 1
                               end do
                                   end do

                                       do i = 1, n-1
                                                  do j = 2, n
                                                             c[i,j] = b[i+1,j] + 1
                                                                     end do
                                                                         end do
                                                                             /* pluto end */
