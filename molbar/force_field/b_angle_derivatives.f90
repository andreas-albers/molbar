module b_angle_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_b_angle_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_b_angle_derivatives(xij,yij,zij,bij,dbdx,d2bdxdy)
    real(8), intent(in) :: xij,yij,zij,bij
    real(8), intent(out) :: dbdx(9),d2bdxdy(45)

    call get_bgradient(xij,yij,zij,bij,dbdx)
    call get_bhessian(xij,yij,zij,bij,d2bdxdy)

end subroutine get_b_angle_derivatives


subroutine get_bhessian(xij,yij,zij,bij,d2bdxdy)
    real(8), intent(in) :: xij,yij,zij, bij
    real(8), intent(out) :: d2bdxdy(45)
    real(8) :: derivative
    integer :: ij

    d2bdxdy = 0.0_wp

    call get_d2bdxidxi(yij, zij, bij, derivative)
    call gidx(1, 1, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(zij, xij, bij, derivative)
    call gidx(2, 2, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(xij, yij, bij, derivative)
    call gidx(3, 3, ij)
    d2bdxdy(ij) = derivative

    call get_d2bdxidyi(xij, yij ,bij, derivative)
    call gidx(1, 2, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xij, zij ,bij, derivative)
    call gidx(1, 3, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(yij, zij ,bij, derivative)
    call gidx(2, 3, ij)
    d2bdxdy(ij) = derivative

    call get_d2bdxidxi(yij, zij, bij, derivative)
    call gidx(4, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(zij, xij, bij, derivative)
    call gidx(5, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(xij, yij, bij, derivative)
    call gidx(6, 6, ij)
    d2bdxdy(ij) = derivative

    call get_d2bdxidyi(xij, yij ,bij, derivative)
    call gidx(4, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xij, zij ,bij, derivative)
    call gidx(4, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(yij, zij ,bij, derivative)
    call gidx(5, 6, ij)
    d2bdxdy(ij) = derivative

    call get_d2bdxidxj(yij, zij,bij, derivative)
    call gidx(1, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxj(zij, xij,bij, derivative)
    call gidx(2, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxj(xij, yij,bij, derivative)
    call gidx(3, 6, ij)
    d2bdxdy(ij) = derivative

    call get_d2bdxidyj(xij, yij ,bij, derivative)
    call gidx(1, 5, ij)
    d2bdxdy(ij) = derivative
    call gidx(2, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(xij, zij ,bij, derivative)
    call gidx(1, 6, ij)
    d2bdxdy(ij) = derivative
    call gidx(3, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(yij, zij ,bij, derivative)
    call gidx(2, 6, ij)
    d2bdxdy(ij) = derivative
    call gidx(3, 5, ij)
    d2bdxdy(ij) = derivative

end subroutine get_bhessian


subroutine get_d2bdxidxi(yij, zij, bij, derivative)

real(8), intent(in) :: yij, zij, bij
real(8), intent(out) :: derivative
derivative = (yij**2+zij**2)/bij**3
end subroutine get_d2bdxidxi

subroutine get_d2bdxidyi(xij, yij ,bij, derivative)

real(8), intent(in) :: xij, yij, bij
real(8), intent(out) :: derivative
derivative = -(xij*yij)/bij**3
end subroutine get_d2bdxidyi

subroutine get_d2bdxidxj(yij, zij,bij, derivative)

real(8), intent(in) :: yij, zij, bij
real(8), intent(out) :: derivative
derivative = -(yij**2+zij**2)/bij**3
end subroutine get_d2bdxidxj

subroutine get_d2bdxidyj(xij, yij ,bij, derivative)

real(8), intent(in) :: xij, yij, bij
real(8), intent(out) :: derivative
derivative = (xij*yij)/bij**3
end subroutine get_d2bdxidyj


subroutine get_bgradient(xij,yij,zij,bij,dbdx)
    real(8), intent(in) :: xij,yij,zij, bij
    real(8), intent(out) :: dbdx(9)
    real(8) :: ddx, ddy, ddz

    dbdx = 0.0_wp

    call get_dbdxi(xij,bij,ddx)
    call get_dbdxi(yij,bij,ddy)
    call get_dbdxi(zij,bij,ddz)

    dbdx(1) = ddx
    dbdx(2) = ddy
    dbdx(3) = ddz 

    call get_dbdxi(xij,bij,ddx)
    call get_dbdxi(yij,bij,ddy)
    call get_dbdxi(zij,bij,ddz)

    dbdx(4) = -ddx
    dbdx(5) = -ddy
    dbdx(6) = -ddz 

end subroutine get_bgradient

subroutine get_dbdxi(xij,bij,derivative)
    real(8), intent(in) :: xij, bij
    real(8), intent(out) :: derivative
    derivative = xij/bij
end subroutine get_dbdxi

subroutine get_dbdxj(xij,bij,derivative)
    real(8), intent(in) :: xij, bij
    real(8), intent(out) :: derivative
    derivative = -xij/bij
end subroutine get_dbdxj
end module b_angle_derivatives
