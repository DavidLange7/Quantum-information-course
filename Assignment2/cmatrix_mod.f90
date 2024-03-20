module cmatrix_mod
  implicit none

  type cmatrix
    integer :: n, m
    double complex, dimension(:,:), allocatable :: array, array_adj
    double complex, allocatable :: array_tr
  end type cmatrix

  interface init
    module procedure cmatrix_init
  end interface

  interface trace
    module procedure cmatrix_trace
  end interface
  
  interface adj
    module procedure cmatrix_adj
  end interface

contains

  subroutine cmatrix_init(self, bin_glob, sq)

    integer bin_glob, bin_temp, i, j, sq
    integer :: n,m 
    type(cmatrix), INTENT(INOUT)  :: self
    real, allocatable :: array_t_r(:), array_t_i2(:)
    double complex, allocatable :: array_t_i(:) 
    
    !First checkpoint of "bin_glob", which has to be given and is 0,1,2
    if (bin_glob /= 0 .and. bin_glob /= 1 .and. bin_glob /= 2) then 
      print *, 'Wrong input of bin_glob, please give 0,1 or 2'
      stop
    end if

    print *, 'Enter dimensions of matrix A (n,m)'
    read *, self%n
    read *, self%m
    n = self%n
    m = self%m

    if (m == n) then 
      sq = 1
    else 
      sq = 0
    end if

    if (allocated(self%array)) then
      deallocate(self%array)
    end if

    allocate(self%array(n, m))

    if (bin_glob == 0) then
      self%array = 1
    end if

    if (bin_glob == 1) then
      print *, 'Is your matrix real (enter 0) or complex (enter 1)?'
      read *, bin_temp
      do i = 1,n
          print *, 'Enter your matrix row wise a,b,c,...:'
          if (bin_temp == 0) then
            allocate(array_t_r(m))
            read *, array_t_r
            do j = 1,m
              self%array(i,j) = complex(array_t_r(j),0)
            end do
            deallocate(array_t_r)
          end if
          if (bin_temp == 1) then
            allocate(array_t_i(m))
            read *, array_t_i
            do j = 1,m
              self%array(i,j) = array_t_i(j)
            end do
            deallocate(array_t_i)
          end if
      end do
    end if

    if (bin_glob == 2) then
      print *, 'You want the random matrix to be real (enter 0) or complex (enter 1)?'
      read *, bin_temp
      do i = 1,n
          if (bin_temp == 0) then
            allocate(array_t_r(m))
            call random_number(array_t_r(m))
            do j = 1,m
              self%array(i,j) = complex(array_t_r(j),0)
            end do
            deallocate(array_t_r)
          end if
          if (bin_temp == 1) then
            allocate(array_t_r(m))
            allocate(array_t_i2(m))
            call random_number(array_t_r(m))
            call random_number(array_t_i2(m))
            do j = 1,m
              self%array(i,j) = complex(array_t_r(j), array_t_i2(j))
            end do
            deallocate(array_t_r)
            deallocate(array_t_i2)
          end if
      end do
    end if

    print *, self%array
  end subroutine cmatrix_init

  subroutine cmatrix_trace(self, sq)
    integer:: n, sq
    integer i
    type(cmatrix), INTENT(INOUT)  :: self
    double complex, allocatable :: array_t(:)
    n = self%n
    
    if (sq /= 1) then
      print *, 'Value Error: For Trace operation a square matrix has to be entered.'
      stop
    end if
    if (.not. allocated(self%array)) then 
      print *, 'Value Error: No matrix initialized'
      stop
    end if

    allocate(array_t(n))
    do i = 1,n
      array_t(i) = self%array(i,i)
    end do
    self%array_tr = sum(array_t)
    print *, '------- BELOW: Matrix Trace'
    print *, self%array_tr
  end subroutine cmatrix_trace

  subroutine cmatrix_adj(self)

    integer:: n, m, i, j
    type(cmatrix), INTENT(INOUT)  :: self
    double complex, allocatable :: array_test(:,:)

    if (.not. allocated(self%array)) then 
      print *, 'Value Error: No matrix initialized'
      stop
    end if
    n = self%n
    m = self%m
    allocate(self%array_adj(m,n))
    allocate(array_test(m,n))


    do i = 1,n
      do j = 1,m
        self%array_adj(j,i) = complex(real(self%array(i,j)), -1*imag(self%array(i,j)))
      end do
    end do

    array_test = conjg(transpose(self%array))

    print *, '---- BELOW: Explicit adjoint'
    print *, self%array_adj
    print *, '------ BELOW: Intrinsic adjoint'
    print *, array_test

  end subroutine cmatrix_adj

  subroutine cmatrix_writetxt(self)
    character(8)  :: date
    character(10) :: time
    integer :: i
    type(cmatrix), INTENT(INOUT)  :: self

    if (allocated(self%array)) then
      call date_and_time(DATE=date,TIME=time)
      open (1, file = 'cmatrix_mod_output.txt', status = 'old', position = 'append')
      write(1,*) '----------Date/time: ',date,'/', time, '----------'
      write (1,*) 'Matrix is given by:'
      do i = 1,self%n
        write (1,*) self%array(i,:)
      end do
      close(1)
    end if

    if (allocated(self%array_tr)) then
      open (1, file = 'cmatrix_mod_output.txt', status = 'old', position = 'append')
      write(1,*) 'The trace of the matrix is given by:'
      write(1,*) self%array_tr
      close(1)
    end if

    if (allocated(self%array_adj)) then
      open (1, file = 'cmatrix_mod_output.txt', status = 'old', position = 'append')
      write(1,*) 'The explicit adjoint calculation yields:'
      do i = 1,self%n
        write (1,*) self%array_adj(i,:)
      end do
      close(1)
      close(1)
    end if

  end subroutine
end module cmatrix_mod

 program exc

  use cmatrix_mod

  implicit none
  integer :: sq !placeholder variable to check if the matrix is square (needed)
  integer :: bin_glob !needed to specify what type of matrix is supplied
  type(cmatrix(:,:)) :: array1
  
  bin_glob = 2
  call cmatrix_init(array1, bin_glob, sq)
  call cmatrix_trace(array1, sq)
  call cmatrix_adj(array1)
  call cmatrix_writetxt(array1)

end program exc