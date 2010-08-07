
/* pluto start (n) */
do k = 1, n
do i = 1, n
do j = 1, n
a[i,j] = a[i,j] - a[i,k]*a[k,j]
end do
end do
end do

/* pluto end */
