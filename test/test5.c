constant n;

do
  i = 1, n do j = 1, n a[i, j] = 0;
end do end do

    do i = n + 1,
       2 *n do j = n + 1, 2 * n b[i, j] = a[i - n, j - n] + 2;
end do end do
