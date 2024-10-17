module d_dihedral_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_d_dihedral_derivatives
integer, parameter :: wp = selected_real_kind(15)

contains

subroutine get_d_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,dddx,d2ddxdy)
     real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
     real(8), intent(out) :: dddx(12), d2ddxdy(78)
     real(8) :: ddx, ddy, ddz
    integer :: i

    call get_dgradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, dddx)

    call get_dhessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2ddxdy)
    
end subroutine get_d_dihedral_derivatives

subroutine get_dhessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2ddxdy)
     real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
     real(8), intent(out) :: d2ddxdy(78)
     real(8) :: derivative
    integer :: ij

    d2ddxdy = 0.0_wp
    ! i i = 0

    !j j
    call get_d2ddxjdxj(xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk, derivative)
    call gidx(4, 4, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdxj(yji, zji, xji, ykj, zkj, xkj, ylk, zlk, xlk,derivative)
    call gidx(5, 5, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdxj(zji, xji, yji, zkj, xkj, ykj, zlk, xlk, ylk, derivative)
    call gidx(6, 6, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyj(xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk, derivative)
    call gidx(4, 5, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyj(zji, xji, yji, zkj, xkj, ykj, zlk, xlk, ylk, derivative)
    call gidx(4, 6, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyj(yji, zji, xji, ykj, zkj, xkj, ylk, zlk,xlk, derivative)
    call gidx(5, 6, ij)
    d2ddxdy(ij) = derivative

    !i j
    call get_d2ddxidxj(xkj,ykj, zkj, ylk, zlk, derivative)
    call gidx(1, 4, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidxj(ykj,zkj, xkj, zlk, xlk, derivative) 
    call gidx(2, 5, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidxj(zkj,xkj,ykj, xlk, ylk, derivative) 
    call gidx(3, 6, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(xkj,ykj,zkj, xlk, ylk, zlk, derivative)
    call gidx(1, 5, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(xkj, zkj, ykj, xlk, zlk, ylk, derivative)
    call gidx(1, 6, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyj(ykj, zkj,xkj, ylk, zlk, xlk, derivative)
    call gidx(2, 6, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(ykj, xkj,zkj, ylk, xlk, zlk, derivative)
    call gidx(2, 4, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyj(zkj, xkj, ykj, zlk, xlk, ylk, derivative)
    call gidx(3, 4, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(zkj, ykj,xkj, zlk, ylk, xlk, derivative)
    call gidx(3, 5, ij)
    d2ddxdy(ij) = -derivative

    ! i k 

    call get_d2ddxidxk(xkj,ykj, zkj, ylk, zlk, derivative)
    call gidx(1, 7, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidxk(ykj,zkj, xkj, zlk, xlk, derivative)
    call gidx(2, 8, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidxk(zkj,xkj, ykj, xlk, ylk, derivative)
    call gidx(3, 9, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(xkj, ykj, zkj, ylk,zlk, derivative)
    call gidx(1, 8, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(xkj, zkj, ykj, zlk,ylk, derivative)
    call gidx(1, 9, ij)
    d2ddxdy(ij) = -derivative 
    call get_d2ddxidyk(ykj, zkj, xkj, zlk, xlk, derivative)
    call gidx(2, 9, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(ykj, xkj, zkj, xlk,zlk,derivative)
    call gidx(2, 7, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyk(zkj, xkj, ykj, xlk,ylk, derivative)
    call gidx(3, 7, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(zkj, ykj, xkj, ylk, xlk, derivative)
    call gidx(3, 8, ij)
    d2ddxdy(ij) = -derivative

    ! j k
    call get_d2ddxjdxk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
    call gidx(4, 7, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdxk(yji,xji,zji,ykj,xkj,zkj,ylk,xlk,zlk,derivative)
    call gidx(5, 8, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdxk(zji,xji,yji,zkj,xkj,ykj,zlk,xlk,ylk,derivative)
    call gidx(6, 9, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
    call gidx(4, 8, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyk(xji,zji,yji,xkj,zkj,ykj,xlk,zlk,ylk,derivative)
    call gidx(4, 9, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdyk(yji,zji,xji,ykj,zkj,xkj,ylk,zlk,xlk,derivative)
    call gidx(5, 9, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyk(yji,xji,zji,ykj,xkj,zkj,ylk,xlk,zlk,derivative)
    call gidx(5, 7, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdyk(zji,xji,yji,zkj,xkj,ykj,zlk,xlk,ylk,derivative)
    call gidx(6, 7, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxjdyk(zji,yji,xji,zkj,ykj,xkj,zlk,ylk,xlk,derivative)
    call gidx(6, 8, ij)
    d2ddxdy(ij) = -derivative

    ! k k 

    call get_d2ddxjdxj(xlk, ylk, zlk, xkj,ykj, zkj, xji, yji, zji, derivative)
    call gidx(7, 7, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdxj(ylk, zlk, xlk, ykj, zkj, xkj, yji, zji, xji, derivative)
    call gidx(8, 8, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdxj(zlk, xlk, ylk, zkj, xkj, ykj, zji, xji, yji, derivative)
    call gidx(9, 9, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdyj(xlk, ylk, zlk,xkj, ykj, zkj, xji, yji, zji, derivative)
    call gidx(7, 8, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdyj(zlk, xlk, ylk,zkj, xkj, ykj,zji, xji, yji,derivative)
    call gidx(7, 9, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxjdyj( ylk, zlk,xlk, ykj, zkj, xkj, yji, zji, xji, derivative)
    call gidx(8, 9, ij)
    d2ddxdy(ij) = -derivative

    ! i l

    call get_d2ddxidyl(xkj,ykj,zkj,derivative)
    call gidx(1, 11, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyl(xkj,zkj,ykj,derivative)
    call gidx(1, 12, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyl(ykj,zkj,xkj,derivative)
    call gidx(2, 12, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyl(ykj,xkj,zkj,derivative)
    call gidx(2, 10, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyl(zkj,xkj,ykj,derivative)
    call gidx(3, 10, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyl(zkj,ykj,xkj,derivative)
    call gidx(3, 11, ij)
    d2ddxdy(ij) = -derivative
    ! j l
    call get_d2ddxidxk(xkj,ykj, zkj, yji, zji, derivative)
    call gidx(4, 10, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidxk(ykj,zkj, xkj, zji, xji, derivative)
    call gidx(5, 11, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidxk(zkj,xkj, ykj, xji, yji, derivative)
    call gidx(6, 12, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyk(ykj, xkj, zkj, xji,zji, derivative)
    call gidx(4, 11, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(zkj, xkj, ykj, xji,yji, derivative)
    call gidx(4, 12, ij)
    d2ddxdy(ij) = -derivative 
    call get_d2ddxidyk(zkj, ykj, xkj, yji, xji, derivative)
    call gidx(5, 12, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(xkj, ykj, zkj, yji,zji,derivative)
    call gidx(5, 10, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyk(xkj, zkj, ykj, zji,yji, derivative)
    call gidx(6, 10, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyk(ykj, zkj, xkj, zji, xji, derivative)
    call gidx(6, 11, ij)
    d2ddxdy(ij) = -derivative

    ! k l
    call get_d2ddxidxj(xkj,ykj, zkj, yji, zji, derivative)
    call gidx(7, 10, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidxj(ykj,zkj, xkj, zji, xji, derivative)
    call gidx(8, 11, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidxj(zkj, xkj, ykj, xji, yji, derivative)
    call gidx(9, 12, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyj(ykj, xkj, zkj, yji, xji, zji, derivative)
    call gidx(7, 11, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(zkj,xkj,ykj, zji, xji, yji, derivative)
    call gidx(7, 12, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyj(zkj,ykj,xkj, zji, yji, xji, derivative)
    call gidx(8, 12, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(xkj,ykj, zkj, xji, yji, zji, derivative)
    call gidx(8, 10, ij)
    d2ddxdy(ij) = -derivative
    call get_d2ddxidyj(xkj,zkj,ykj, xji, zji, yji, derivative)
    call gidx(9, 10, ij)
    d2ddxdy(ij) = derivative
    call get_d2ddxidyj(ykj,zkj,xkj, yji, zji, xji, derivative)
    call gidx(9, 11, ij)
    d2ddxdy(ij) = -derivative 

    ! l l = 0            
    
end subroutine get_dhessian

subroutine get_d2ddxjdxj(xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk, derivative)
 real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk
 real(8), intent(out) :: derivative
 real(8) :: enumerator, denominator
enumerator = - 2*xkj*(xkj**2 + ykj**2 + zkj**2)*(yji*zlk + ykj*zlk - ylk*zji - ylk*zkj) &
            + (ykj**2 + zkj**2)*(xji*(ykj*zlk - ylk*zkj) - yji*(xkj*zlk - xlk*zkj) &
            + zji*(xkj*ylk - xlk*ykj))
denominator = (xkj**2 + ykj**2 + zkj**2)**1.5_wp
derivative = enumerator/denominator
end subroutine get_d2ddxjdxj

subroutine get_d2ddxjdyj(xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk, derivative)
 real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk
 real(8), intent(out) :: derivative
 real(8) :: enumerator, denominator
enumerator = - xkj*ykj*(xji*(ykj*zlk - ylk*zkj) - yji*(xkj*zlk - xlk*zkj) + &
             zji*(xkj*ylk - xlk*ykj)) - (- xkj*(xji*zlk + xkj*zlk - xlk*zji - xlk*zkj) &
            + ykj*(yji*zlk + ykj*zlk - ylk*zji - ylk*zkj))*(xkj**2 + ykj**2 + zkj**2)
denominator = (xkj**2 + ykj**2 + zkj**2)**1.5_wp
derivative = enumerator/denominator
end subroutine get_d2ddxjdyj

subroutine get_d2ddxidxj(xkj,ykj, zkj, ylk, zlk, derivative)
 real(8), intent(in) :: xkj, ykj, zkj, ylk, zlk
 real(8), intent(out) :: derivative
derivative = xkj*(ykj*zlk - ylk*zkj)/(sqrt(xkj**2 + ykj**2 + zkj**2))
end subroutine get_d2ddxidxj

subroutine get_d2ddxidyj(xkj, ykj,zkj, xlk, ylk, zlk, derivative)
intrinsic :: sqrt
 real(8), intent(in) :: xkj, ykj,zkj, xlk, ylk, zlk
 real(8), intent(out) :: derivative
derivative = (ykj*(ykj*zlk - ylk*zkj) + zlk*(xkj**2 + ykj**2 + zkj**2))/(sqrt(xkj**2 + ykj**2 + zkj**2))
end subroutine get_d2ddxidyj

subroutine get_d2ddxidxk(xkj,ykj, zkj, ylk, zlk, derivative)
intrinsic :: sqrt
 real(8), intent(in) :: xkj,ykj, zkj, ylk, zlk
 real(8), intent(out) :: derivative
derivative = xkj*(- ykj*zlk + ylk*zkj)/sqrt(xkj**2 + ykj**2 + zkj**2)
end subroutine get_d2ddxidxk

subroutine get_d2ddxidyk(xkj, ykj, zkj, ylk,zlk, derivative)
 real(8), intent(in) :: xkj, ykj, zkj, ylk,zlk
 real(8), intent(out) :: derivative
derivative = (- ykj*(ykj*zlk - ylk*zkj) - (zkj + zlk)*(xkj**2 + ykj**2 + zkj**2))/(sqrt(xkj**2 + ykj**2 + zkj**2))
end subroutine get_d2ddxidyk


subroutine get_d2ddxidyl(xkj,ykj,zkj,derivative)
 real(8), intent(in) :: xkj,ykj,zkj
 real(8), intent(out) :: derivative
derivative = zkj*sqrt(xkj**2+ykj**2+zkj**2)
end subroutine get_d2ddxidyl

subroutine get_d2ddxjdxk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
intrinsic :: sqrt
 real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
 real(8), intent(out) :: derivative
 real(8) :: enumerator, denominator
enumerator = xkj**2*(xji*(ykj*zlk - ylk*zkj) - yji*(xkj*zlk - xlk*zkj) + zji*(xkj*ylk - xlk*ykj)) &
            - (xkj**2 + ykj**2 + zkj**2)*(xji*(ykj*zlk - ylk*zkj) + xkj*(- yji*(zkj + zlk) + zji*(ykj + ylk)) & 
            - xkj*(yji*zlk + ykj*zlk - ylk*zji - ylk*zkj) - yji*(xkj*zlk - xlk*zkj) + zji*(xkj*ylk - xlk*ykj)) 
denominator = (xkj**2 + ykj**2 + zkj**2)**1.5_wp
derivative = enumerator/denominator
end subroutine get_d2ddxjdxk

subroutine get_d2ddxjdyk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
intrinsic :: sqrt
 real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
 real(8), intent(out) :: derivative
 real(8) :: enumerator, denominator
enumerator = xkj*ykj*(xji*(ykj*zlk - ylk*zkj) - yji*(xkj*zlk - xlk*zkj) + zji*(xkj*ylk - xlk*ykj)) &
             + (zji + zkj + zlk)*(xkj**2 + ykj**2 + zkj**2)**2 + (- xkj*(xji*(zkj + zlk) - zji*(xkj + xlk)) &
            + ykj*(yji*zlk + ykj*zlk - ylk*zji - ylk*zkj))*(xkj**2 + ykj**2 + zkj**2)
denominator = (xkj**2 + ykj**2 + zkj**2)**1.5_wp
derivative = enumerator/denominator
end subroutine get_d2ddxjdyk


subroutine get_dgradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, dddx)
     real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
     real(8), intent(out) :: dddx(12)
     real(8) :: ddx, ddy, ddz

    dddx = 0.0_wp

    call get_dddxi(xkj, ykj, zkj, ylk, zlk, ddx)
    call get_dddxi(ykj, zkj, xkj, zlk, xlk, ddy)
    call get_dddxi(zkj, xkj, ykj, xlk, ylk, ddz)

    dddx(1) = ddx
    dddx(2) = ddy
    dddx(3) = ddz 

    call get_dddxj(xkj, ykj, zkj, xlk,ylk, zlk, xji, yji, zji, ddx)
    call get_dddxj(ykj, zkj, xkj, ylk,zlk, xlk, yji, zji, xji, ddy)
    call get_dddxj(zkj, xkj, ykj, zlk, xlk, ylk, zji, xji, yji, ddz)

    dddx(4) = ddx
    dddx(5) = ddy
    dddx(6) = ddz 

    call get_dddxk(xkj, ykj, zkj, xlk, ylk, zlk, xji, yji, zji, ddx)
    call get_dddxk(ykj, zkj, xkj, ylk, zlk, xlk, yji, zji, xji, ddy)
    call get_dddxk(zkj, xkj, ykj, zlk, xlk, ylk, zji, xji, yji, ddz)

    dddx(7) = ddx
    dddx(8) = ddy
    dddx(9) = ddz 

    call get_dddxl(xkj, ykj, zkj, xji, yji, zji, ylk, zlk, ddx)
    call get_dddxl(ykj, zkj, xkj, yji, zji, xji, zlk, xlk, ddy)
    call get_dddxl(zkj, xkj, ykj, zji, xji, yji, xlk, ylk, ddz)

    dddx(10) = ddx
    dddx(11) = ddy
    dddx(12) = ddz 

end subroutine get_dgradient

subroutine get_dddxi(xkj, ykj, zkj, ylk, zlk, derivative)
    intrinsic :: sqrt
     real(8), intent(in) :: xkj, ykj, zkj, ylk, zlk
     real(8), intent(out) :: derivative

    derivative = -(ykj*zlk - ylk*zkj)*sqrt(xkj**2 + ykj**2 + zkj**2)

end subroutine get_dddxi

subroutine get_dddxj(xkj, ykj, zkj, xlk,ylk, zlk, xji, yji, zji, derivative)
    intrinsic :: sqrt
     real(8), intent(in) :: xkj, ykj, zkj, xlk, ylk, zlk, xji, yji, zji
     real(8), intent(out) :: derivative

    derivative = - (xkj * (xji * (ykj * zlk - ylk * zkj) + yji * (-xkj * zlk + xlk * zkj) + zji * (xkj * ylk - xlk * ykj))) / &
            sqrt(xkj**2 + ykj**2 + zkj**2) - (-yji * zlk + ylk * zji) * sqrt(xkj**2 + ykj**2 + zkj**2) +&
            (ykj * zlk - ylk * zkj) * sqrt(xkj**2 + ykj**2 + zkj**2)

end subroutine get_dddxj

subroutine get_dddxk(xkj, ykj, zkj, xlk, ylk, zlk, xji, yji, zji, derivative)
     real(8), intent(in) :: xkj, ykj, zkj, xlk, ylk, zlk, xji, yji, zji
     real(8), intent(out) :: derivative
    
    derivative = (xkj*(xji*(ykj*zlk - ylk*zkj) + yji*(-xkj*zlk + xlk*zkj) &
                + zji*(xkj*ylk - xlk*ykj)))/sqrt(xkj**2 + ykj**2 + zkj**2) - &
                (yji*zkj - ykj*zji)*sqrt(xkj**2 + ykj**2 + zkj**2) + (ylk*zji - yji*zlk)*sqrt(xkj**2 + ykj**2 + zkj**2)
end subroutine get_dddxk

subroutine get_dddxl(xkj, ykj, zkj, xji, yji, zji, ylk, zlk, derivative)
     real(8), intent(in) :: xkj, ykj, zkj, xji, yji, zji, ylk, zlk
     real(8), intent(out) :: derivative
    derivative = (yji * zkj - ykj * zji) * sqrt(xkj**2 + ykj**2 + zkj**2)
end subroutine get_dddxl

end module d_dihedral_derivatives