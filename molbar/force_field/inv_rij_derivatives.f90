module inv_rij_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_inverse_distance_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_inverse_distance_derivatives(xij,yij,zij,rij,dqdx,d2qdxdy)
    real(8), intent(in) :: xij,yij,zij,rij
    real(8), intent(out) :: dqdx(6),d2qdxdy(21)

    call get_qgradient(xij,yij,zij,rij,dqdx)
    call get_qhessian(xij,yij,zij,rij,d2qdxdy)

end subroutine get_inverse_distance_derivatives


subroutine get_qhessian(xij,yij,zij,rij,d2qdxdy)
    real(8), intent(in) :: xij,yij,zij, rij
    real(8), intent(out) :: d2qdxdy(21)
    real(8) :: derivative
    integer :: ij

    d2qdxdy = 0.0_wp

    call get_d2qdxidxi(xij, yij, zij, rij, derivative)
    call gidx(1, 1, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxi(yij, zij, xij, rij, derivative)
    call gidx(2, 2, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxi(zij, xij, yij, rij, derivative)
    call gidx(3, 3, ij)
    d2qdxdy(ij) = derivative

    call get_d2qdxidyi(xij, yij, rij, derivative)
    call gidx(1, 2, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyi(xij, zij ,rij, derivative)
    call gidx(1, 3, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyi(yij, zij ,rij, derivative)
    call gidx(2, 3, ij)
    d2qdxdy(ij) = derivative

    call get_d2qdxidxi(xij, yij, zij, rij, derivative)
    call gidx(4, 4, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxi(yij, zij, xij, rij, derivative)
    call gidx(5, 5, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxi(zij, xij, yij, rij, derivative)
    call gidx(6, 6, ij)
    d2qdxdy(ij) = derivative

    call get_d2qdxidyi(xij, yij, rij, derivative)
    call gidx(4, 5, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyi(xij, zij ,rij, derivative)
    call gidx(4, 6, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyi(yij, zij ,rij, derivative)
    call gidx(5, 6, ij)
    d2qdxdy(ij) = derivative

    call get_d2qdxidxj(xij, yij, zij,rij, derivative)
    call gidx(1, 4, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxj(yij, zij, xij,rij, derivative)
    call gidx(2, 5, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidxj(zij, xij, yij,rij, derivative)
    call gidx(3, 6, ij)
    d2qdxdy(ij) = derivative

    call get_d2qdxidyj(xij, yij ,rij, derivative)
    call gidx(1, 5, ij)
    d2qdxdy(ij) = derivative
    call gidx(2, 4, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyj(xij, zij ,rij, derivative)
    call gidx(1, 6, ij)
    d2qdxdy(ij) = derivative
    call gidx(3, 4, ij)
    d2qdxdy(ij) = derivative
    call get_d2qdxidyj(yij, zij ,rij, derivative)
    call gidx(2, 6, ij)
    d2qdxdy(ij) = derivative
    call gidx(3, 5, ij)
    d2qdxdy(ij) = derivative

end subroutine get_qhessian

subroutine get_d2qdxidxi(xij, yij, zij, rij, derivative)
real(8), intent(in) :: xij,yij, zij,rij
real(8), intent(out) :: derivative
derivative = (2*xij**2-yij**2-zij**2)/rij**5
end subroutine get_d2qdxidxi

subroutine get_d2qdxidyi(xij, yij, rij, derivative)
real(8), intent(in) :: xij, yij, rij
real(8), intent(out) :: derivative
derivative = (3*xij*yij)/rij**5
end subroutine get_d2qdxidyi

subroutine get_d2qdxidxj(xij, yij, zij,rij, derivative)
real(8), intent(in) :: xij, yij, zij, rij
real(8), intent(out) :: derivative
derivative = (-2*xij**2+yij**2+zij**2)/rij**5
end subroutine get_d2qdxidxj

subroutine get_d2qdxidyj(xij, yij,rij, derivative)
real(8), intent(in) :: xij, yij, rij
real(8), intent(out) :: derivative
derivative = -(3*xij*yij)/rij**5
end subroutine get_d2qdxidyj

subroutine get_qgradient(xij,yij,zij,rij,dqdx)
    real(8), intent(in) :: xij,yij,zij, rij
    real(8), intent(out) :: dqdx(6)
    real(8) :: ddx, ddy, ddz

    dqdx = 0.0_wp

    call get_dqdxi(xij,rij,ddx)
    call get_dqdxi(yij,rij,ddy)
    call get_dqdxi(zij,rij,ddz)

    dqdx(1) = ddx
    dqdx(2) = ddy
    dqdx(3) = ddz 

    call get_dqdxi(xij,rij,ddx)
    call get_dqdxi(yij,rij,ddy)
    call get_dqdxi(zij,rij,ddz)

    dqdx(4) = -ddx
    dqdx(5) = -ddy
    dqdx(6) = -ddz 

end subroutine get_qgradient

subroutine get_dqdxi(xij,rij,derivative)
    real(8), intent(in) :: xij, rij
    real(8), intent(out) :: derivative
    derivative = -xij/rij**3
end subroutine get_dqdxi

end module inv_rij_derivatives
