module oscill    
    implicit none

    type store
        double precision, allocatable :: Yl(:), Yr(:), xmesh(:), Z(:,:), E(:)
        integer :: conv
    end type store

contains

    subroutine finitediff( self, n_mesh, xmax, xmin )
        implicit none

        type(store), INTENT(out) :: self
        integer :: n_mesh, i, LDZ, INFO
        double precision, allocatable :: xmesh(:), d(:), sd(:), Z(:,:), WORK(:)
        double precision :: xmax, xmin, h, hbar, m , omega

        !Energy in our scale is N+0.5 !!
        !self%E = E !we should first guess one, and then have something more sophisticated to find the right eigenvalues.. 
        self%conv = 0
        
        hbar = 1.0
        m = 1.0
        omega = 1.0
        LDZ = n_mesh
        allocate(xmesh(n_mesh))
        allocate(d(n_mesh))
        allocate(sd(n_mesh-1))
        allocate(WORK(max(1,2*n_mesh-2)))
        allocate(Z(LDZ, n_mesh))

        xmesh(1) = xmin
        h = (xmax-xmin)/DBLE(n_mesh)

        do i=1,n_mesh
            xmesh(i+1) =  xmesh(i) + h
        end do

        sd(:) = -hbar*hbar/(2*m*h*h)

        do i=1,n_mesh
            d(i) = hbar*hbar/(2.0*m*h*h)*2 + 1.0/2*m*omega*omega*xmesh(i)*xmesh(i)
        end do

        call dsteqr('I', n_mesh, d, sd, Z, LDZ, WORK, INFO)

        self%E = d
        self%Z = Z(:,:)

        open (1, file = 'Z.txt', status = 'old')
        write(1,*) Z(:,1)
        do i=2,30
            open (1, file = 'Z.txt', status = 'old', position = 'append')
            write(1,*) Z(:,i)
        end do
        close(1)
        open (2, file = 'mesh.txt', status = 'old')
        write(2,*) xmesh
        close(2)

        open (3, file = 'E.txt', status = 'old')
        write(3,*) d
        close(3)
    end subroutine finitediff

end module oscill

Program dostuff
    use oscill
    Implicit none

    integer :: n_mesh
    double precision :: xmax, xmin
    type(store):: val
    !integer :: n_mesh, i
    !double precision :: xmax, xmin, E
    n_mesh = 250
    xmax = 6.0
    xmin = -6.0

    call finitediff( val, n_mesh, xmax, xmin )
end program dostuff