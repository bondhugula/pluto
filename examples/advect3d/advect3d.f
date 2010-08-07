      constant nx, ny, nz, T, nbdy;

        do j = 7,ny+nbdy-5
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-6
               abthird(j,i,k) = (f60 * (athird(j-1,i,k) + athird(j,i,k))
                   + f61 * (athird(j-2,i,k) + athird(j+1,i,k))
                   + f62 * (athird(j-3,i,k) + athird(j+2,i,k)))
                   * halfdtbydx * uybthird(j,i,k)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-5
            do k = 7,nz+nbdy-6
               althird(j,i,k) = (f60 * (athird(j,i-1,k) + athird(j,i,k))
                   + f61 * (athird(j,i-2,k) + athird(j,i+1,k))
                   + f62 * (athird(j,i-3,k) + athird(j,i+2,k)))
                   * halfdtbydx * uxlthird(j,i,k)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-5
               afthird(j,i,k) = (f60 * (athird(j,i,k-1) + athird(j,i,k))
                   + f61 * (athird(j,i,k-2) + athird(j,i,k+1))
                   + f62 * (athird(j,i,k-3) + athird(j,i,k+2)))
                   * halfdtbydx * uzfthird(j,i,k)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-6
               ahalf(j,i,k) = a(j,i,k) + (althird(j,i+1,k) 
                   - althird(j,i,k))
                   + (abthird(j+1,i,k) - abthird(j,i,k))
                   + (afthird(j,i,k+1) - afthird(j,i,k))
            end do
         end do
      end do
