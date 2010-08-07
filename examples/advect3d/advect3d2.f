      constant nx, ny, nz, nbdy;

      do j = -2,ny+nbdy
         do i = -2,nx+nbdy
            do k = -2,nz+nbdy
               af(k,i,j) = 1
            end do
         end do
      end do


      do j = -2,ny+nbdy
         do i = -2,nx+nbdy
            do k = -2,nz+nbdy
               athird(k,i,j) = af(k+1,i,j)
            end do
         end do
      end do
