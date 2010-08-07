
      subroutine twostages(a, ahalf, uxl, uyb, uzf, uxlthird, uybthird, 
     & uzfthird,
     &     dt, dx, dy, dz, nx, ny, nz, nbdy)

!      implicit none

      dimension a(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension ahalf(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uxl(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uyb(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uzf(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uxlthird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uybthird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension uzfthird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)

c     Local storage follows:
      dimension ab(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension al(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension af(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension athird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension abthird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension althird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)
      dimension afthird(1-nbdy:nz+nbdy, 1-nbdy:nx+nbdy, 1-nbdy:ny+nbdy)

      real a, ahalf, uxl, uyb, uzf, uxlthird, uybthird, uzfthird, ab
      real al, af, athird, abthird, althird, afthird, dt, dx, dy
      real dz, f60, f61, thirddtbydy, halfdtbydy, thirddtbydx
      real halfdtbydx, thirddtbydz, halfdtbydz, f62

      integer nx, ny, nz, nbdy, i, j, k

      f60 = 1.2
      f61 = 1.5
      f62 = 2.2

      thirddtbydy = dt / (3. * dy)
      halfdtbydy = 1.5 * thirddtbydy
      thirddtbydx = dt / (3. * dx)
      halfdtbydx = 1.5 * thirddtbydx
      thirddtbydz = dt / (3. * dz)
      halfdtbydz = 1.5 * thirddtbydz

      do j = 4-nbdy,ny+nbdy-2
         do i = 4-nbdy,nx+nbdy-3
            do k = 4-nbdy,nz+nbdy-3
               ab(k,i,j) = (f60 * (a(k,i,j-1) + a(k,i,j)) + f61 
     & * (a(k,i,j-2) + a(k,i,j+1)) + f62 * (a(k,i,j-3) + a(k,i,j+2))) 
     & * thirddtbydy * uyb(k,i,j)
            enddo
         enddo
      enddo

      do j = 4-nbdy,ny+nbdy-3
         do i = 4-nbdy,nx+nbdy-2
            do k = 4-nbdy,nz+nbdy-3
               al(k,i,j) = (f60 * (a(k,i-1,j) + a(k,i,j)) + f61 
     & * (a(k,i-2,j) + a(k,i+1,j)) + f62 * (a(k,i-3,j) + a(k,i+2,j)))
     & * thirddtbydx * uxl(k,i,j)
            enddo
         enddo
      enddo

      do j = 4-nbdy,ny+nbdy-3
         do i = 4-nbdy,nx+nbdy-3
            do k = 4-nbdy,nz+nbdy-2
               af(k,i,j) = (f60 * (a(k-1,i,j) + a(k,i,j)) + f61 
     & * (a(k-2,i,j) + a(k+1,i,j)) + f62 * (a(k-3,i,j) + a(k+2,i,j)))
     & * thirddtbydz * uzf(k,i,j)
            enddo
         enddo
      enddo


      do j = 4-nbdy,ny+nbdy-3
         do i = 4-nbdy,nx+nbdy-3
            do k = 4-nbdy,nz+nbdy-3
               athird(k,i,j) = a(k,i,j) + (al(k,i+1,j) - al(k,i,j))
     & + (ab(k,i,j+1) - ab(k,i,j)) + (af(k+1,i,j) - af(k,i,j))
            enddo
         enddo
      enddo


      do j = 7-nbdy,ny+nbdy-5
         do i = 7-nbdy,nx+nbdy-6
            do k = 7-nbdy,nz+nbdy-6
               abthird(k,i,j) = (f60 * (athird(k,i,j-1) + athird(k,i,j))
     &              + f61 * (athird(k,i,j-2) + athird(k,i,j+1))
     &              + f62 * (athird(k,i,j-3) + athird(k,i,j+2)))
     &              * halfdtbydx * uybthird(k,i,j)
            enddo
         enddo
      enddo


      do j = 7-nbdy,ny+nbdy-6
         do i = 7-nbdy,nx+nbdy-5
            do k = 7-nbdy,nz+nbdy-6
               althird(k,i,j) = (f60 * (athird(k,i-1,j) + athird(k,i,j))
     &              + f61 * (athird(k,i-2,j) + athird(k,i+1,j))
     &              + f62 * (athird(k,i-3,j) + athird(k,i+2,j)))
     &              * halfdtbydx * uxlthird(k,i,j)
            enddo
         enddo
      enddo


      do j = 7-nbdy,ny+nbdy-6
         do i = 7-nbdy,nx+nbdy-6
            do k = 7-nbdy,nz+nbdy-5
               afthird(k,i,j) = (f60 * (athird(k-1,i,j) + athird(k,i,j))
     &              + f61 * (athird(k-2,i,j) + athird(k+1,i,j))
     &              + f62 * (athird(k-3,i,j) + athird(k+2,i,j)))
     &              * halfdtbydx * uzfthird(k,i,j)
            enddo
         enddo
      enddo


      do j = 7-nbdy,ny+nbdy-6
         do i = 7-nbdy,nx+nbdy-6
            do k = 7-nbdy,nz+nbdy-6
               ahalf(k,i,j) = a(k,i,j) + (althird(k,i+1,j) 
     &              - althird(k,i,j))
     &              + (abthird(k,i,j+1) - abthird(k,i,j))
     &              + (afthird(k+1,i,j) - afthird(k,i,j))
            enddo
         enddo
      enddo
      end






