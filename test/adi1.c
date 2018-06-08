constant n;

do t = 1, n

           do i = 1, n
                      do j = 1, n
                                 a[i,j] = a[i,j-1] + 1
                                          end do
                                              end do

                                                  do i = 1, n
                                                             do j = 1, n
                                                                         a[i,n-j] = a[i,n-j+1] + 1
                                                                                 end do
                                                                                     end do

                                                                                         do i = 1, n
                                                                                                     do j = 1, n
                                                                                                                 a[i,j] = a[i-1,j] + 1
                                                                                                                         end do
                                                                                                                             end do

                                                                                                                                 do i = 1, n
                                                                                                                                             do j = 1, n
                                                                                                                                                         a[i,n-j] = a[i-1,n-j+1] + 1
                                                                                                                                                                 end do
                                                                                                                                                                     end do

                                                                                                                                                                         end do
