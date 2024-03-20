module diagstuff
    implicit none

    type store
        DOUBLE PRECISION, allocatable :: W(:)
        integer :: N
    end type store

contains
    subroutine rand_eigen(self, N, LDA)
        Implicit none

        type(store), INTENT(out) :: self
        INTEGER :: INFO, INFO2, LWORK, LDA, N, K, j, i, glob
        INTEGER :: ISEED(4)
        DOUBLE PRECISION, allocatable :: D(:), W(:), RWORK(:)
        DOUBLE PRECISION :: temp1, temp2
        COMPLEX*16, allocatable :: A(:,:), WORK(:), WORK2(:)
        self%N = N
        LWORK = -1
        allocate(A(LDA, N))
        allocate(WORK2(3000))
        allocate(D(N))
        allocate(WORK(2*N))
        allocate(RWORK(3*N-2))
        allocate(W(N))

        !can also unse an intrinsic function to create a complex hermitian matrix
        !K = 1
        !ISEED = (/122, 178, 1, 898/)
        !call ZLAGHE(N, K, D, A, LDA, ISEED, WORK, INFO)

        glob = 0 !1 if you want hermitian, for diagonal some other integer
        
            A = 0.0
            do i=1,N
                call random_number(temp1)
                temp1 = -1 + 2*temp1
                A(i,i) = complex(temp1,0)
            end do
            if (glob == 1) then
                do i = 2, N
                    do j = 1, i-1
                        call random_number(temp1)
                        call random_number(temp2)
                        temp1 = -1 + 2*temp1
                        temp2 = -1 + 2*temp2
                        A(i,j) = complex(temp1,temp2)
                        A(j,i) = conjg(A(i,j))
                    end do
                end do
            end if

        call zheev('N','L', N, A, LDA, W, WORK2, LWORK, RWORK, INFO2)

        LWORK = int(WORK2(1))
        deallocate(WORK2)
        allocate(WORK2(LWORK))
        call zheev('N','L', N, A, LDA, W, WORK2, LWORK, RWORK, INFO2)
        self%W = W
    end subroutine rand_eigen

    subroutine store_eigenvalues(self)
        Implicit none

        type(store), intent(in) :: self
        !character (len = *) :: filename
        open (1, file = 'eigenvalues.txt', status = 'old')
        write(1,*) self%W
        close(1)
    end subroutine store_eigenvalues

    subroutine store_normspacing(self)
        Implicit none

        type(store), intent(in) :: self
        integer :: i
        DOUBLE PRECISION :: S(self%N-1), avrg
        do i=1,self%N-1
            S(i) = (self%W(i+1)-self%W(i))
        end do
        avrg = sum(S)/(self%N-1)
        S = S/avrg
        open(1, file = 'normspacing.txt', status = 'old', position = 'append')
        write(1,*) S
        close(1)
    end subroutine store_normspacing
end module diagstuff

Program test
    use diagstuff
    Implicit none

    integer :: N, LDA, i
    type(store) :: val
    N = 1500
    LDA = 1500
    do i=1,50
        call rand_eigen(val, N, LDA)
        call store_eigenvalues(val)
        call store_normspacing(val)
    end do
end program test