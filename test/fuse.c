/* pluto start (n) */

do i = 1, n

           do j = 1, n
                      a[2*i] = a[2*i] + b[j]
                               end do

                                   do j = 1, n
                                              a[2*i+1] = a[2*i+1] + b[j]
                                                      end do

                                                          end do

                                                              /* pluto end */
