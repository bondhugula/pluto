
/* pluto start (n) */

do i = 4 -n, n + 2
do j = 4 -n, n + 2
 a[i,j] = a[i-1,j] + 2

end do
end do

do i = 4 -n, n + 2
do j = 4 -n, n + 2
  a[i,j] = a[i,j] + 1
end do
end do

/* pluto end */
