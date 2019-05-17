constant n;

do
i = 1, n do j = 1,
            n b[i] = b[i] + a[i, j] end do end do

                     do i = 1,
            n do j = 1, n c[i] = c[i] + a[j, i] end do end do
