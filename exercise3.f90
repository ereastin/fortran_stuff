module funcs
    implicit none
    real, parameter :: pi = 3.1415
    real, parameter :: grav_const = 6.672e-11
    real, parameter :: rho = 3000.0
    real, parameter :: dom_width = 1000.0e3
    integer, parameter :: res = 101

contains
    subroutine p1()
        real :: x(res), R(res), L(res, res), kernel(3), step
        integer :: i
        step = dom_width / res
        x = [((i - 1) * step, i = 1, res)]
        R(:) = 0
        R(2:res-1) = 4 * pi * grav_const * rho
        L(:,:) = 0
        do i = 2, res - 1
            L(i, i:i+size(kernel)-1) = kernel
        end do
        L(1,1) = 1
        L(res,res) = 1
    end subroutine p1
end module funcs

program main
    use funcs
    implicit none
    call p1()
end program main
