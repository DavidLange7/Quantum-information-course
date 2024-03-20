module debugger
    implicit none
contains
    subroutine db_print_real(var, debug)
        real :: var
        integer :: debug
        if (debug == 1) then
            print *, var
        end if
    end subroutine

    subroutine db_print_imag(var, debug)
        double complex :: var
        integer :: debug
        if (debug == 1) then
            print *, var
        end if
    end subroutine

    subroutine db_print_int(var, debug, name)
        integer :: debug, var
        character (len = *) :: name
        if (debug == 1) then
            print *, 'Debugger level 1: variable ', name, ' =', var
        end if

        if (debug == 2) then
            print *, 'Debugger level 2: variable ', name, ' =', var
        end if

    end subroutine
end module

program exc2
    use debugger
    implicit none 
    INTEGER :: bin_glob, dlevel, n_mat1, i, j, k
    real :: timearray1(2), timearray2(2), timearray3(2), timediff1, timediff2, timediff3
    real, allocatable :: array2(:,:), array3(:,:), array4(:,:)

    dlevel = 0
	!Following is for user input, needed for all cases

    open(unit = 1, file = "settings.txt")
    read(1,*) n_mat1
    close(1)

    call db_print_int(bin_glob, dlevel, 'bin_glob')
    allocate(array2(n_mat1, n_mat1))
    allocate(array3(n_mat1, n_mat1))

    do i=1,n_mat1
        do j = 1, n_mat1
            array2(i, j) = i+j
        end do
    end do

    do i=1,n_mat1
        do j = 1, n_mat1
            array3(i, j) = i*j
        end do
    end do
	
	!FIRST LOOP ORDER #################################
    timediff1 = dtime( timearray1 )
    allocate(array4(n_mat1, n_mat1))
    array4(:,:) = 0

    do i = 1,n_mat1
            do k = 1,n_mat1
                do j = 1,n_mat1
                array4(i,j) = array4(i,j) + array2(i,k)*array3(k,j)
            end do
        end do
    end do
    
    timediff1 = dtime( timearray1 )
    open (2, file = 'one.txt', status = 'old', position='append')
    write(2,*) timediff1
    close(2)
    deallocate(array4)
	

	!SECOND LOOP ORDER #####################
    timediff2 = dtime( timearray2 )
    allocate(array4(n_mat1, n_mat1))
    array4(:,:) = 0

    do j = 1,n_mat1
            do k = 1,n_mat1
                do i = 1,n_mat1
                array4(i,j) = array4(i,j) + array2(i,k)*array3(k,j)
            end do
        end do
    end do
    timediff2 = dtime( timearray2 )

    open (1, file = 'two.txt', status = 'old', position='append')
    write(1,*) timediff2
    close(1)
    deallocate(array4)

	!INTRINSIC MATRIX MULTIPLICATION ########################
    timediff3 = dtime( timearray3 )
    allocate(array4(n_mat1, n_mat1))
    array4 = matmul(array2,array3)
    timediff3 = dtime( timearray3 )
    open (1, file = 'three.txt', status = 'old', position='append')
    write(1,*) timediff3
    close(1)
end program exc2