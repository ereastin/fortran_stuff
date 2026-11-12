module funcs
    use, intrinsic :: iso_fortran_env, only: l => int64, d => real64
    implicit none
    real, parameter :: tmin = 0.0, tmax = 1350.0
    real(d), parameter :: kappa = 1.0e-6  ! m2/s, thermal diffusivity
    real(d), parameter :: dom_width = 150.0e3  ! m
    integer(l), parameter :: yr_to_sec = 3600 * 24 * 365
    integer(l), parameter :: time = 15e6 * yr_to_sec  ! 15 Myr
    integer, parameter :: space_res = 151
    real(d), parameter :: dx = real(dom_width / (space_res - 1))

contains
    subroutine analytic(u)
        real(d), intent(inout) :: u(:,:)
        integer :: i

        do i = 1, size(u, 1)
            u(i, 1) = tmax + (tmin - tmax) * erf((i - 1) * dx / (2 * sqrt(kappa * time)))
        end do
    end subroutine analytic

    subroutine forward_euler(u, dt)
        real(d), intent(inout) :: u(:,:)
        integer(l), intent(in) :: dt
        real(d), allocatable :: A(:,:)
        real(d) :: kernel(3), r
        integer :: i

        integer(l) :: N
        N = time / dt

        allocate(A(space_res, space_res))
        A(:, :) = 0.0
        r = kappa * dt / dx ** 2
        kernel = [r, 1-2*r, r]

        do i = 2, space_res - 1
            A(i, i-1:i+1) = kernel
        end do

        A(1, 1) = 1
        A(space_res, space_res) = 1

        do i = 1, N
            u = matmul(A, u)
        end do

        deallocate(A)
    end subroutine forward_euler

    subroutine backward_euler(u, dt)
        real(d), intent(inout) :: u(:,:)
        integer(l), intent(in) :: dt
        real(d), allocatable :: A(:,:), Ainv(:,:)
        real(d) :: kernel(3), r
        integer :: i

        integer(l) :: N
        N = time / dt

        allocate(A(space_res, space_res))
        allocate(Ainv(space_res, space_res))

        A(:,:) = 0.0
        r = kappa * dt / dx ** 2
        kernel = [-r, 1+2*r, -r]

        do i = 2, space_res - 1
            A(i, i-1:i+1) = kernel
        end do

        A(1, 1) = 1
        A(space_res, space_res) = 1
        Ainv = invert_matrix(A)

        do i = 1, N
            u = matmul(Ainv, u)
        end do

        deallocate(A)
        deallocate(Ainv)
    end subroutine backward_euler

    subroutine crank_nicolson(u, dt, vary_kappa)
        real(d), intent(inout) :: u(:,:)
        integer(l), intent(in) :: dt
        logical, intent(in), optional :: vary_kappa
        real(d), allocatable:: A(:,:), Ainv(:,:), B(:,:), C(:,:)
        real(d) :: kernel_A(3), kernel_B(3), kappa_scale(space_res), r, sr
        integer :: i
        integer(l) :: N

        N = time / dt
        kappa_scale(:) = 1.0
        if (present(vary_kappa)) then
            kappa_scale(1:30) = 10.0  ! increase first 30 km near tmax boundary
        end if

        ! allocate arrays
        allocate(A(space_res, space_res))
        allocate(Ainv(space_res, space_res))
        allocate(B(space_res, space_res))
        allocate(C(space_res, space_res))

        ! initialize arrays
        A(:,:) = 0.0
        B(:,:) = 0.0

        r = kappa * dt / dx ** 2
        !kernel_A = [-r/2, 1+r, -r/2]
        !kernel_B = [r/2, 1-r, r/2]

        do i = 2, space_res - 1
            sr = r * kappa_scale(i)
            kernel_A = [-sr/2, 1+sr, -sr/2]
            kernel_B = [sr/2, 1-sr, sr/2]
            A(i, i-1:i+1) = kernel_A
            B(i, i-1:i+1) = kernel_B
        end do

        A(1,1) = 1.0
        A(space_res, space_res) = 1.0
        B(1,1) = 1.0
        B(space_res, space_res) = 1.0

        ! do computation
        Ainv = invert_matrix(A)
        C = matmul(Ainv, B)

        do i = 1, N
            u = matmul(C, u)
        end do

        ! deallocate
        deallocate(A)
        deallocate(Ainv)
        deallocate(B)
        deallocate(C)
    end subroutine crank_nicolson

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
        u(2:space_res, 1) = tmin
        u(1, 1) = tmax
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
    real(d) :: u(space_res, 1)
    integer(l) :: dt(3)
    character(len=128) :: fname
    integer :: i
    i = 1
    dt = [4e11, 5e11, 6e11]

    call init_uvec(u)
    call analytic(u)
    write(fname, '(A)') 'analytic.txt'
    call write_file(trim(fname), u)

    do i = 1, 3
        call init_uvec(u)
        call forward_euler(u, dt(i))
        write(fname, '(A,I0,A)') 'forward_euler', i, '.txt'
        call write_file(trim(fname), u)
    end do

    do i = 1, 3
        call init_uvec(u)
        call backward_euler(u, dt(i))
        write(fname, '(A,I0,A)') 'backward_euler', i, '.txt'
        call write_file(fname, u)
    end do

    do i = 1, 3
        call init_uvec(u)
        call crank_nicolson(u, dt(i))
        write(fname, '(A,I0,A)') 'crank_nicolson', i, '.txt'
        call write_file(fname, u)
    end do

    call init_uvec(u)
    call crank_nicolson(u, dt(1), .true.)
    write(fname, '(A)') 'crank_nicolson_Kx.txt'
    call write_file(fname, u)
end program main
