module diagstuff
    implicit none

    type store
        DOUBLE PRECISION, allocatable :: temp1(:,:), w(:), Mop(:,:), rato
        complex*16, allocatable :: e_vecs(:,:)
        integer :: N, dim_1, dim_2, nzers
    end type store

contains

    subroutine tp(self, a, b)
        Implicit none

        integer :: k, l, m, n, i, j, o, p, t1, t2
        type(store), intent(out) :: self
        double precision, allocatable :: a(:,:), b(:,:), t(:,:)
        k = size(a, dim = 1)
        l = size(a, dim = 2)

        m = size(b, dim = 1)
        n = size(b, dim = 2)
        allocate(t(k*m, l*n))

        do i=0,k-1
            do j=0,l-1
                do o=1,m
                    do p=1,n
                        t( m*i+o , n*j+p ) = a(i+1,j+1)*b(o,p)
                    enddo
                enddo
            enddo
        enddo

        if (allocated(self%temp1)) deallocate(self%temp1)
        allocate(self%temp1(k*m, l*n))
        self%temp1 = t

        self%dim_1 = k*m
        self%dim_2 = l*n

    end subroutine

    subroutine ising(self, lambd, sides)
        Implicit none
      
        integer :: N, LDA, i,j,k,l,m,o,p,q, sides, z, INFO2, LWORK
        logical :: temp1
        INTEGER :: ISEED(4)
        DOUBLE PRECISION, allocatable :: D(:), W(:), RWORK(:), Mopr(:,:)
        double precision, allocatable :: pau_x(:,:), pau_z(:,:), id(:,:), tmp(:,:), tmp2(:,:), spns(:,:)
        double precision :: lambd
        COMPLEX*16, allocatable :: WORK(:), WORK2(:), H(:,:)
        type(store), intent(inout):: self

        allocate(pau_x(2,2))
        allocate(pau_z(2,2))
        allocate(id(2,2))

        pau_x = reshape((/ 0, 1, 1, 0  /), shape(pau_x))
        pau_z = reshape((/ 1, 0, 0, -1  /), shape(pau_z))
        id = reshape((/ 1, 0, 0, 1  /), shape(id))

        allocate(H(2**sides,2**sides))
        allocate(Mopr(2**sides,2**sides))
        Mopr = 0
        H = 0
        !allocate(tmp(sides*sides,sides*sides))
        do i=1,sides
            if (i .eq. 1) then
                allocate(tmp(2,2))
                    tmp = pau_z
            else
                allocate(tmp(2,2))
                tmp = id
            endif

            do j=2,sides
                if (j .eq. i .and. i /= 1) then
                    call tp(self, tmp, pau_z)
                else 
                    call tp(self, tmp, id)
                endif
                deallocate(tmp)
                allocate(tmp(self%dim_1, self%dim_2))
                tmp = self%temp1
            enddo
            Mopr = Mopr + tmp
            H = H + tmp
            deallocate(tmp)
        enddo
        H = lambd*H
        Mopr = Mopr/7

        do i=1,sides-1
            temp1 = .False.
            allocate(tmp2(2,2))
            if (i .eq. 1) then 
                tmp2 = pau_x
            else 
                tmp2 = id
            endif

            do j=2,sides
                if (j .eq. i .and. i /= 1) then
                    call tp(self, tmp2, pau_x)
                    temp1 = .True.
                    
                else if (temp1 .eqv. .True.) then
                    call tp(self, tmp2, pau_x)
                    temp1 = .False.

                else if (i == 1 .and. j == 2) then 
                    call tp(self, tmp2, pau_x)
                    temp1 = .False.

                else 
                    call tp(self, tmp2, id)
                endif

                deallocate(tmp2)
                allocate(tmp2(self%dim_1, self%dim_2))
                tmp2 = self%temp1                
            enddo
            H = H + tmp2
            deallocate(tmp2)
        enddo
        self%nzers = count(H/=0)
        self%rato = ((2**sides*2**sides)-self%nzers)/(2**sides*2**sides)

        LDA = size(H, dim=1)
        N = size(H, dim=2)
        LWORK = -1
        allocate(WORK2(3000))
        allocate(D(N))
        allocate(WORK(2*N))
        allocate(RWORK(3*N-2))
        allocate(W(N))

        call zheev('V','L', N, H, LDA, W, WORK2, LWORK, RWORK, INFO2)
        LWORK = int(WORK2(1))
        deallocate(WORK2)
        allocate(WORK2(LWORK))
        call zheev('V','L', N, H, LDA, W, WORK2, LWORK, RWORK, INFO2)
        self%W = W
        self%e_vecs = H
        self%Mop = Mopr
        !print *, W

    end subroutine ising

    subroutine txt_stuff(self)
        Implicit none

        integer :: i
        type(store), intent(in) :: self
        !character (len = *) :: filename
        open (1, file = 'eigenvalues.txt', status = 'old')
        write(1,*) self%W
        close(1)
        open (2, file = 'eigenvectors.txt', status = 'old')
        write(2,*) Real(self%e_vecs)
        close(2)

    end subroutine txt_stuff

    subroutine txt_changelambda(self)
        Implicit none
        
        integer :: i, sides, j
        double precision, allocatable :: tmp(:,:), M(:,:)
        type(store), intent(inout) :: self
        !double precision :: lamda(5) = (/0.0 ,0.1 ,1.0 ,3.0 , 10.0/)
        double precision :: lamda !lamda(300) = (/((i*1.0)/30, i=0,299, 1)/)
        !character (len = *) :: filename
        open (1, file = 'magnetization.txt', status = 'old')
        open (2, file = 'gs.txt', status = 'old')

        lamda = 0.01
        sides = 7
        !do i=1,300
            call ising(self, lamda, sides)
            tmp = Real(self%e_vecs)
            M = self%Mop
            do i=1,size(M,dim=1)
                write(1,*) M(i,:)
            enddo
            write(2,*) tmp(1,:)
        !enddo
        close(1)
        close(2)

    end subroutine

end module diagstuff


Program test

    use diagstuff
    Implicit none

    integer :: sides, i
    double precision :: lamda
    type(store) :: val
    double precision :: a(300) = (/(i/100, i=1,300, 1)/)

    lamda = 1
    sides = 6
    print*, '-----------'
    call ising(val, lamda, sides)
    !call txt_stuff(val)
    !call txt_changelambda(val)

    !call tp(val, ons, pau_x)
    !print *, val%temp1
    !print *, '--------------'
    !print *, val%dim_1, val%dim_2

end program test