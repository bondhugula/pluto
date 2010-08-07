constant n;

do i = 1, n 
do j = 1, n
    a[i,j] = a[i-1,j+1] + 1
end do
end do
