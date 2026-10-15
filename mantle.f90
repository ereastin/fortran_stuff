program mantle
    implicit none
    integer, parameter :: dp = kind(0.0d0)
    integer, parameter :: S = 31
    integer :: i, j
    real(dp), dimension(S, S) :: vx, vy, dvx_dx, dvx_dy, dvy_dy, dvy_dx, div
    real(dp) :: vx0, vy0, pi, H, W, stepx, stepy

    vx0 = 1e-9
    vy0 = 1e-9
    W = 1e6
    H = 1.5e6
    pi = 3.1415
    stepx = W / S
    stepy = H / S

    do i = 1, S
        do j =1, S
            vx(i, j) = 
            vy(i, j) = 
