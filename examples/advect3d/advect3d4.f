      constant nx, ny, nz, nbdy;

      do j = 4,ny+nbdy-2
         do i = 4,nx+nbdy-3
            do k = 4,nz+nbdy-3
            ab(k,i,j) = (f60 * (a(k,i,j-1) + a(k,i,j)) + f61 * (a(k,i,j-2) + a(k,i,j+1)) + f62 * (a(k,i,j-3) + a(k,i,j+2))) * thirddtbydy * uyb(k,i,j)
            end do
         end do
      end do

      do j = 4,ny+nbdy-3
         do i = 4,nx+nbdy-2
            do k = 4,nz+nbdy-3
               al(k,i,j) = (f60 * (a(k,i-1,j) + a(k,i,j)) + f61 
     * (a(k,i-2,j) + a(k,i+1,j)) + f62 * (a(k,i-3,j) + a(k,i+2,j)))
     * thirddtbydx * uxl(k,i,j)
            end do
         end do
      end do

      do j = 4,ny+nbdy-3
         do i = 4,nx+nbdy-3
            do k = 4,nz+nbdy-2
               af(k,i,j) = (f60 * (a(k-1,i,j) + a(k,i,j)) + f61 
     * (a(k-2,i,j) + a(k+1,i,j)) + f62 * (a(k-3,i,j) + a(k+2,i,j)))
      * thirddtbydz * uzf(k,i,j)
            end do
         end do
      end do


      do j = 4,ny+nbdy-3
         do i = 4,nx+nbdy-3
            do k = 4,nz+nbdy-3
               athird(k,i,j) = a(k,i,j) + (al(k,i+1,j) - al(k,i,j))
      + (ab(k,i,j+1) - ab(k,i,j)) + (af(k+1,i,j) - af(k,i,j))
            end do
         end do
      end do

        do j = 7,ny+nbdy-5
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-6
               abthird(k,i,j) = (f60 * (athird(k,i,j-1) + athird(k,i,j))
                   + f61 * (athird(k,i,j-2) + athird(k,i,j+1))
                   + f62 * (athird(k,i,j-3) + athird(k,i,j+2)))
                   * halfdtbydx * uybthird(k,i,j)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-5
            do k = 7,nz+nbdy-6
               althird(k,i,j) = (f60 * (athird(k,i-1,j) + athird(k,i,j))
                   + f61 * (athird(k,i-2,j) + athird(k,i+1,j))
                   + f62 * (athird(k,i-3,j) + athird(k,i+2,j)))
                   * halfdtbydx * uxlthird(k,i,j)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-5
               afthird(k,i,j) = (f60 * (athird(k-1,i,j) + athird(k,i,j))
                   + f61 * (athird(k-2,i,j) + athird(k+1,i,j))
                   + f62 * (athird(k-3,i,j) + athird(k+2,i,j)))
                   * halfdtbydx * uzfthird(k,i,j)
            end do
         end do
      end do


      do j = 7,ny+nbdy-6
         do i = 7,nx+nbdy-6
            do k = 7,nz+nbdy-6
               ahalf(k,i,j) = a(k,i,j) + (althird(k,i+1,j) 
                   - althird(k,i,j))
                   + (abthird(k,i,j+1) - abthird(k,i,j))
                   + (afthird(k+1,i,j) - afthird(k,i,j))
            end do
         end do
      end do
