module funcs
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: RES = 50
    real, parameter :: xs = 0.0, xe = 2.0
    real :: dx
    dx = (xe - xs) / RES

contains
    function u(x) result(out)
        real, intent(in) :: x
        real :: out
        out = exp(x ** 2)
    end function u

    function du(x) result(out)
        real, intent(in) :: x
        real :: out
        out = 2 * x * exp(x ** 2)
    end function du

    function Dforward(x0, x1) result(out)
        real, intent(in) :: x0, x1
        out = (u(x0) - u(x1)) / dx
    end function Dforward

    function Dcentered(x0, x1) result(out)
        real, intent(in) :: x0, x1
        out = (u(x0) - u(x1)) / (2 * dx)
    end function Dcentered 

    subroutine plot()
        real :: df(RES), df_forward(RES - 1), df_centered(RES - 2), x
        x = [(i - 1) * dx, i = 1, RES]

        do i = 1, RES
            df(i) = du(x(i))
            write(*,*) x(i), df(i)
        end do
        do i = 1, RES - 1
            df_forward(i) = Dforward(x(i + 1), x(i))
            write(*,*) x(i), df_forward(i)
        end do
        do i = 2, RES - 1
            df_centered(i) = Dcentered(x(i + 1), x(i - 1))
            write(*,*) x(i), df_centered(i)
        end do
    end subroutine plot
end module funcs

program main
    call plot
end program main
