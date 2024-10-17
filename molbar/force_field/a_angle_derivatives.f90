module a_angle_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_a_angle_derivatives
integer, parameter :: wp = selected_real_kind(15)

contains

subroutine get_a_angle_derivatives(xij,yij,zij,xkj,ykj,zkj,dadx,d2adxdy)
    real(8), intent(in) :: xij,yij,zij,xkj,ykj,zkj
    real(8), intent(out) :: dadx(9), d2adxdy(45)
    real(8) :: ddx, ddy, ddz
    integer :: i

    call get_agradient(xij,yij,zij,xkj,ykj,zkj,dadx)

    call get_ahessian(xij,yij,zij,xkj,ykj,zkj,d2adxdy)
    
end subroutine get_a_angle_derivatives


subroutine get_ahessian(xij,yij,zij,xkj,ykj,zkj,d2adxdy)
    real(8), intent(in) :: xij,yij,zij,xkj,ykj,zkj
    real(8), intent(out) :: d2adxdy(45)
    real(8) :: derivative
    integer :: ij

    d2adxdy = 0.0_wp

    call gidx(1,7, ij)
    d2adxdy(ij) = 1.0_wp

    call gidx(2,8, ij)
    d2adxdy(ij) = 1.0_wp

    call gidx(3,9, ij)
    d2adxdy(ij) = 1.0_wp

    call gidx(1,4, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(2,5, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(3,6, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(4,7, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(5,8, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(6,9, ij)
    d2adxdy(ij) = -1.0_wp

    call gidx(4,4, ij)
    d2adxdy(ij) = 2.0_wp

    call gidx(5,5, ij)
    d2adxdy(ij) = 2.0_wp

    call gidx(6,6, ij)
    d2adxdy(ij) = 2.0_wp


end subroutine get_ahessian


subroutine get_agradient(xij,yij,zij,xkj,ykj,zkj,dadx)
    real(8), intent(in) :: xij,yij,zij,xkj,ykj,zkj
    real(8), intent(out) :: dadx(9)
    real(8) :: ddx, ddy, ddz

    dadx = 0.0_wp

    call get_dadxi(xkj,ddx)
    call get_dadxi(ykj,ddy)
    call get_dadxi(zkj,ddz)

    dadx(1) = ddx
    dadx(2) = ddy
    dadx(3) = ddz 
    
    call get_dadxj(xij, xkj, ddx)
    call get_dadxj(yij, ykj, ddy)
    call get_dadxj(zij, zkj, ddz)

    dadx(4) = ddx
    dadx(5) = ddy
    dadx(6) = ddz 

    call get_dadxi(xij,ddx)
    call get_dadxi(yij,ddy)
    call get_dadxi(zij,ddz)

    dadx(7) = ddx
    dadx(8) = ddy
    dadx(9) = ddz

end subroutine get_agradient


subroutine get_dadxi(xkj,derivative)
real(8), intent(in) :: xkj
real(8), intent(out) :: derivative
derivative = xkj
end subroutine get_dadxi

subroutine get_dadxj(xij, xkj, derivative)
real(8), intent(in) :: xij, xkj
real(8), intent(out) :: derivative
derivative = -xij-xkj
end subroutine get_dadxj


end module a_angle_derivatives
