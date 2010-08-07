            program main
            DOUBLE precision A(2048,2048), B(2048,2048), C(2048,2048),alpha,beta
            integer i, j
            M = 2048
            N = 2048
            K = 2048
            LDA = 2048
            LDB = 2048
            LDC = 2048
            alpha = 1.0
            beta = 1.0

            do i = 1, 2048
            do j = 1, 2048
            a(i,j) = i;
            b(i,j) = j;
            c(i,j) = 0
            end do
            end do

            CALL DGEMM('T','N',M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)

            end
