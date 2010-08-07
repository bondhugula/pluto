            program main
            DOUBLE precision A(250,250), B(250,250),u1(250)
            DOUBLE precision u2(250), v1(250), v2(250)
            DOUBLE precision w(250), x(250), y(250), z(250)
            DOUBLE precision t_start, t_end

            integer i, j
            m = 250;
            n = 250;

            alpha = 1.0
            beta = 1.0

            do i = 1, 250
            do j = 1, 250
            B(i,j) = 1;
            end do
            u1(i)  = 1;
            u2(i)  = 1;
            v1(i)  = 1;
            v2(i)  = 1;
            w(i)  = 1;
            x(i)  = 1;
            y(i)  = 1;
            z(i)  = 1;
            end do

            call rtclock (t_start)

            call dcopy(m * n, A, 1, B, 1);
            call dger(m, n, 1.0, u1, 1, v1 , 1, B, m);
            call dger(m, n, 1.0, u2, 1, v2 , 1, B, m);
            call dcopy(n,z,1,x,1)
            call dgemv('T', m, n, beta, B, m, y, 1, 1.0, x, 1);
            call dgemv('N', m, n, alpha, B, m, x, 1, 0.0, w, 1);
            call rtclock (t_end)
            write(*, 900) t_end-t_start
900         format (F10.6)
            end
