
/* pluto start (N) */
do i = 2,N
    a[i] = b[i] + c[i]
    b[i+2] = a[i-1] + c[i-1]
    a[i+1] = b[i+3] + 1
end do
/* pluto end */
