program matrix
	implicit none
	INTEGER*2 :: bin_glob, n_mat1, m_mat1, n_mat2, m_mat2, i, j, k
	real :: timearray1(2), timearray2(2), timearray3(2), timediff1, timediff2, timediff3
	real, allocatable :: array1(:)
	real, allocatable :: array2(:,:), array3(:,:), array4(:,:)
	
	!Following is for user input, needed for all cases
	print *, 'Type 1 If manual matrix otherwise 0:'
	read *, bin_glob
	
	print *, 'Enter dimensions of matrix A (n,m)'
	read *, n_mat1
	read *, m_mat1
	allocate(array2(n_mat1, m_mat1))

        print *, 'Enter dimensions of matrix B (n,m)'
        read *, n_mat2
        read *, m_mat2
        allocate(array3(n_mat2, m_mat2))

	!bin_glob is the variable which determines 2 things: bin_glob = 1 means you enter your matrix and also you get all results printed on screen, 
	!bin_glob = 0 means the matrices are initialized automatically with one being rowindex+columnindex and the other rowindex*columnindex
		!bin_glob = 0 also means the resulting computation times are written on three file, while the other results are "discarded", since not needed.
	if (bin_glob == 1) then
		do i = 1,n_mat1
			allocate(array1(m_mat1))
			print *, 'Enter your matrix row wise a,b,c,...:'
			read *, array1
			do j = 1,m_mat1
				array2(i,j) = array1(j)
			end do 
			deallocate(array1)
		end do
		print *, array2

        	do i = 1,n_mat2
                	allocate(array1(m_mat2))
                	print *, 'Enter your matrix row wise a,b,c,...:'
                	read *, array1
                	do j = 1,m_mat2
                        	array3(i,j) = array1(j)
                	end do
                	deallocate(array1)
        	end do
        	print *, array3
	else
		do i=1,n_mat1
      			do j = 1, m_mat1
         			array2(i, j) = i+j
      			end do
   		end do

                do i=1,n_mat2
                        do j = 1, m_mat2
                                array3(i, j) = i*j
                        end do
                end do

	end if
	! FROM HERE MATRIX MULTIPLICATION

	! First Raise Value Error if matrix dimensions dont fit multiplication
	if (m_mat1 /= n_mat2) then 
		print *, 'Wrong matrix dimensions, multiplication not possible!'
		stop
	end if
	
	!FIRST LOOP ORDER #################################
	timediff1 = dtime( timearray1 )
	allocate(array4(n_mat1, m_mat2))
	array4(:,:) = 0

	do i = 1,n_mat1
		do j = 1,m_mat2
			do k = 1,n_mat2
				array4(i,j) = array4(i,j) + array2(i,k)*array3(k,j)
			end do
		end do
	end do
	timediff1 = dtime( timearray1 )
	
	print *, '---------------------------------------'
	print *, 'End result of matrix multiplication by do do do (rowwise) loop is: (time for computation was', timediff1, 'seconds)'
	print *, '------------------------------------------------'

        if (bin_glob == 1) then
		print *,array4
        	print *, '---------------------------------------'
	end if
	
	if (bin_glob /= 1) then
		open (1, file = 'one.txt', status = 'old', position='append')
		write(1,*) timediff1
		close(1)
	end if 
	deallocate(array4)
	

	!SECOND LOOP ORDER #####################
        timediff3 = dtime( timearray3 )
        allocate(array4(n_mat1, m_mat2))
        array4(:,:) = 0

        do j = 1,m_mat2
                do i = 1,n_mat1
                        do k = 1,n_mat2
                                array4(i,j) = array4(i,j) + array2(i,k)*array3(k,j)
                        end do
                end do
        end do
        timediff3 = dtime( timearray3 )
	print *, 'End result of matrix multplication by do do do (column wise) loop is: (time for computation was', timediff3, 'seconds)'
        print *, '------------------------------------------------'

	if (bin_glob == 1) then
		print *,array4
                print *, '---------------------------------------'
        end if

        if (bin_glob /= 1) then
                open (1, file = 'two.txt', status = 'old', position='append')
                write(1,*) timediff3
                close(1)
        end if
        deallocate(array4)

	!INTRINSIC MATRIX MULTIPLICATION ########################
        timediff2 = dtime( timearray3 )
	allocate(array4(n_mat1, m_mat2))
	array4 = matmul(array2,array3)
        timediff2 = dtime( timearray2 )
	print *, 'End result of matrix multiplication by matmul intrinsic function is:(time for computation was', timediff2, 'seconds)'
        print *, '------------------------------------------------'

        if (bin_glob == 1) then
                print *,array4
                print *, '---------------------------------------'
        end if
        if (bin_glob /= 1) then
                open (1, file = 'three.txt', status = 'old', position='append')
                write(1,*) timediff2
                close(1)
        end if
end program matrix
