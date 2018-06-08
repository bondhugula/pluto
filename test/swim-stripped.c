constant n;

do t = 2, n
           do i = 1, n
                      do j = 1, n
                                 a[t,i,j] = a[t-1,i-1,j] + a[t-1,i,j+1] + a[t-1,i+1,j] + a[t-1,i,j-1]
                                            + a[t-1,i,j]
                                            end do
                                                end do

                                                    do j = 1,n
                                                               a[t,0,j] = a[t,n,j]
                                                                       end do

                                                                           do i = 1, n
                                                                                       a[t,i,n+1] = a[t,i,1]
                                                                                               end do

                                                                                                   end do
