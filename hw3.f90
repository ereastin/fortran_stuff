module funcs
    use, intrinsic :: iso_fortran_env, only: l => int64, d => real64
    implicit none
    integer, parameter :: nx = 31 + 1, ny = 21 + 1  ! extra nodes
    real, parameter :: wx = 1000.0e3, wy = 1500.0e3  ! m; x/y domain widths
    real(d), parameter :: eta = 10.0e21  ! Pa/s; viscosity
    real(d), parameter :: grav = 9.81  ! m/s2; gravity accel
    real(d), parameter :: rho_l = 3200.0, rho_r = 3300.0  ! kg/m3; density for left/right sides
    real(d), parameter :: dx = real(wx / (nx - 1)), dy = real(wy / (ny - 1))
    real(d), parameter :: Kscale = 2 * eta / (dx + dy)

contains
    !! have it accept boundary conditions :: need one version for free slip and another for no slip BCs
    subroutine solve_stokes(r)
        real(d), intent(inout) :: r(:,:)
        integer(l), intent(in) :: dt
        real(d), allocatable:: A(:,:), Ainv(:,:), coeffs_mass_cons(:), coeffs_xstokes(:), coeffs_ystokes(:)
        integer :: i, j, d, ivx, ivy, ip
        d = size(r, 1)

        ! allocate arrays
        allocate(A(d, d))
        allocate(Ainv(d, d))
        allocate(coeffs_mass_cons(d))
        allocate(coeffs_xstokes(d))
        allocate(coeffs_ystokes(d))

        ! initialize arrays
        A(:,:) = 0.0

        do i = 1, ny
            do j = 1, nx
            ! get index for vx, vy, and P at each node, use figure
            coeffs_mass_cons =  ! how to build this..
            A(i, :) = coeffs
            end do
        end do

        ! boundary conditions - FIX
        A(1,1) = 1.0
        A(d, d) = 1.0

        ! do computation
        Ainv = invert_matrix(A)
        u = matmul(Ainv, u)

        ! deallocate
        deallocate(A)
        deallocate(Ainv)
    end subroutine solve_stokes

    function invert_matrix(A) result(Ainv)
        real(d), intent(in) :: A(:, :)
        real(d) :: Ainv(size(A, 1), size(A, 2))
        real(d), allocatable :: work(:)
        integer(l), allocatable :: pivot(:)  ! apparently LAPACK routines expect long ints?
        integer(l) :: n, info

        external DGETRF
        external DGETRI

        Ainv = A
        n = size(A, 1)
        allocate(work(n))
        allocate(pivot(n))

        call DGETRF(n, n, Ainv, n, pivot, info)
        if (info /= 0) then
            stop 'matrix is numerically singular'
        end if

        call DGETRI(n, Ainv, n, pivot, work, n, info)
        if (info /= 0) then
            stop 'matrix inversion failed'
        end if

        deallocate(work)
        deallocate(pivot)
    end function invert_matrix

    subroutine init_uvec(u)
        real(d), intent(inout) :: u(:,:)
        u(:,:) = 0.0
    end subroutine init_uvec

    subroutine write_file(fname, arr)
        character(len=*), intent(in) :: fname
        real(d), intent(in) :: arr(:,:)
        integer :: i, funit

        open(newunit=funit, file=fname, status='replace', action='write')
        do i = 1, size(arr, 1)
            write(funit, *) dx * (i - 1), arr(i,1)
        end do
        close(funit)
    end subroutine write_file
end module funcs

program main
    use funcs
    implicit none
    real(d) :: u(3 * nx * ny, 1)  ! 3 unknowns; vx, vy, P at all nodes
    character(len=128) :: fname
    integer :: i

    call init_uvec(u)
    call analytic(u)
    write(fname, '(A)') 'analytic.txt'
    call write_file(trim(fname), u)
end program main
