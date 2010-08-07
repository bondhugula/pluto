
      /* pluto start (nx,ny,nz,nbdy) */

      do j = 4-nbdy,ny+nbdy-3
         do i = 4-nbdy,nx+nbdy-3
            do k = 4-nbdy,nz+nbdy-2
               af(k,i,j) = (f60 * (a(k-1,i,j) + a(k,i,j)) + f61 
     * (a(k-2,i,j) + a(k+1,i,j)) + f62 * (a(k-3,i,j) + a(k+2,i,j)))
      * thirddtbydz * uzf(k,i,j)
            end do
         end do
      end do


      do j = 4-nbdy,ny+nbdy-3
         do i = 4-nbdy,nx+nbdy-3
            do k = 4-nbdy,nz+nbdy-3
               athird(k,i,j) = a(k,i,j) + (al(k,i+1,j) - al(k,i,j))
      + (ab(k,i,j+1) - ab(k,i,j)) + (af(k+1,i,j) - af(k,i,j))
            end do
         end do
      end do
      /* pluto end */
