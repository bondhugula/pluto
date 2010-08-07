      PROGRAM ADVECT3D
      IMPLICIT NONE
      INTEGER nx, ng, tt

      parameter (nx = 256)
      parameter (ng = 9)
      parameter (tt = 1)
        REAL x(1-ng:nx+ng), y(1-ng:nx+ng), z(1-ng:nx+ng)
	REAL dx, dis
	REAL pi, Ad1, Ad2
	REAL  A0(1-ng:nx+1+ng, 1-ng:nx+1+ng)
	REAL  An(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+  ng)
	REAL  Ao(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+  ng)
        REAL  Ai(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+  ng) 
	REAL   u(1-ng:nx  +ng, 1-ng:nx+1+ng,1-ng:nx+  ng)
	REAL   v(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+1+ng)
	REAL   w(1-ng:nx+1+ng, 1-ng:nx+  ng,1-ng:nx+  ng)
	REAL  u0(1-ng:nx  +ng, 1-ng:nx+1+ng,1-ng:nx+  ng)
	REAL  v0(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+1+ng)
	REAL  w0(1-ng:nx+1+ng, 1-ng:nx+  ng,1-ng:nx+  ng)
	REAL  ut(1-ng:nx  +ng, 1-ng:nx+1+ng,1-ng:nx+  ng)
	REAL  vt(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+1+ng)
	REAL  wt(1-ng:nx+1+ng, 1-ng:nx+  ng,1-ng:nx+  ng)
	REAL  uh(1-ng:nx  +ng, 1-ng:nx+1+ng,1-ng:nx+  ng)
	REAL  vh(1-ng:nx  +ng, 1-ng:nx+  ng,1-ng:nx+1+ng)
	REAL  wh(1-ng:nx+1+ng, 1-ng:nx+  ng,1-ng:nx+  ng)
        REAL X0, Y0, Z0, rad
        REAL dt
        REAL maxi
	INTEGER i,j,k,t,m, imgT
        INTEGER iloc, jloc, kloc
	real xl, yl, zl, third, half
        PARAMETER (dt = 0.005)
        PARAMETER (pi = 3.1415926)
        PARAMETER (X0 = 0.00)
	PARAMETER (Y0 = 0.00)
	PARAMETER (Z0 = 0.00)
        PARAMETER (rad = 0.15)
        imgT = Tt/4
        dx = 1./REAL(nx)
        maxi = 0.0
	third = 1./3.
	half = 0.5

	do i=1-ng,nx+ng
	   x(i) = -.5 + dx*REAL(i)
	   y(i) = -.5 + dx*REAL(i)
	   z(i) = -.5 + dx*REAL(i)
	end do
	xl = REAL(nx) * dx * .5
	yl = REAL(nx) * dx * .5
	zl = REAL(nx) * dx * .5

        do j=1,nx
          do i=1,nx
            do k=1,nx
              u0(k,i,j) =  1.0 * (x(k) - 0.5*dx)
              v0(k,i,j) = 0.
              w0(k,i,j) = -1.0 * (z(i) - 0.5*dx)
            end do
          end do
        end do

        Ad1 = 0.
	 do j=1,nx	
           do i=1,nx
	     do k=1,nx	
	     dis = SQRT( (x(i) -.5*dx - X0)**2 +
     +                   (y(j) -.5*dx - Y0)**2 +
     +                   (z(k) -.5*dx - Z0)**2)
                Ao(k,i,j) = (20. * 0.5 * (1. + COS(dis*pi/rad)))*
     +                    (.5+sign(.5,rad-dis))
              Ai(k,i,j) = Ao(k,i,j)
              An(k,i,j) = 0.0
              Ad1 = Ad1 + Ao(k,i,j)
	    end do
	  end do
	end do
	do i=1-ng,nx+ng
           do k=1-ng,nx+ng
              do m = 0, 1-ng
                 Ao(k,i,m) = Ao(k,i,1)
              enddo
              do m = nx + 1, nx + ng
                 Ao(k,i,m) = Ao(k,i,nx)
              enddo
           enddo
        enddo
        do j=1-ng,nx+ng
           do k=1-ng,nx+ng
	   do m = 0, 1-ng
              Ao(k,m,j) = Ao(k,1 ,j)
	   enddo
           do m = nx+1, nx+ng
              Ao(k,m,j) = Ao(k,nx,j)
           enddo
        enddo

           do i=1,nx
	   do k=1-ng,nx+ng
              maxi = MAX(maxi,Ao(k,i,j))
	      if(Ao(k,i,j) .EQ. maxi) then
		 iloc = i
		 jloc = j
		 kloc = k
	      endif
           end do
           end do
        end do
         write(*,*) maxi, 'max val start ',iloc,jloc,kloc
        maxi = 0.0 
	print *, '( nx, ng, tt ) = ', '(', nx, ',', ng, ',', tt, ')'
        write (*,*) Tt
        do t=1, 5
            call twostages(Ao, An, u0, v0, w0, u0, v0, w0,
     &                      dt, dx, dx, dx, nx, nx, nx, ng)

            do j = 1, nx
               do i = 1, nx
                  do k = 1, nx
                     Ao(k,i,j) = An(k,i,j)
                  enddo
               enddo
            enddo
        end do
        write(*,*) An(kloc,jloc,iloc), ' val @ finish for initial max '
        Ad2 = 0.
        do j=1,nx
          do i=1,nx
	   do k=1-ng,nx+ng
              maxi = MAX(maxi,Ao(k,i,j))
	      if(Ao(k,i,j) .EQ. maxi) then
		 iloc = i
		 jloc = j
		 kloc = k
	      endif
	      Ad2 = Ad2 + An(k,i,j)
           end do
          end do
        end do
        write(*,*) maxi, 'max val finish ',iloc,jloc
        maxi = 0.0 
        do j = 1, nx
           do i = 1, nx
              do k = 1, nx
                 Ai(k, i, j) = An(k, i, k) - Ai(k, i, j)
              enddo
           enddo
        enddo
        do j=1,nx
          do i=1,nx
	    do k=1-ng,nx+ng
              maxi = MAX(maxi,Ai(k,i,j))
	      if(Ai(k,i,j) .EQ. maxi) then
		 iloc = i
		 jloc = j
		 kloc = k
	      endif
           end do
          end do
        end do
         write(*,*) 'max error finish: ',maxi,iloc,jloc
         write(*,*) 'Starting amount:  ',Ad1
         write(*,*) 'Finishing amount: ',Ad2
         write(*,*) 'Percent conserved: ',Ad2/Ad1 * 100.

       END 

