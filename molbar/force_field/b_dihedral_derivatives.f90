module b_dihedral_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_b_dihedral_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains


subroutine get_b_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, bijk, dbdx,d2bdxdy)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,bijk
    real(8), intent(out) :: dbdx(12),d2bdxdy(78)
    real(8) :: ddx, ddy, ddz

    call get_bgradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, bijk, dbdx)
    call get_bhessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2bdxdy)

end subroutine get_b_dihedral_derivatives


subroutine get_bhessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2bdxdy)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
    real(8), intent(out) :: d2bdxdy(78)
    real(8) :: derivative
    integer :: ij

    d2bdxdy = 0.0_wp

    ! i i
    call get_d2bdxidxi(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 1, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 2, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(3, 3, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 2, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xji, zji, yji, xkj, zkj, ykj, derivative)
    call gidx(1, 3, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 3, ij)
    d2bdxdy(ij) = derivative

    !i j
    call get_d2bdxidxj(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxj(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxj(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(3, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(xji, zji, yji, xkj, zkj, ykj, derivative)
    call gidx(1, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(yji, xji, zji, ykj, xkj, zkj, derivative)
    call gidx(2, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(3, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyj(zji, yji, xji, zkj, ykj, xkj, derivative)
    call gidx(3, 5, ij)
    d2bdxdy(ij) = derivative

    !j j
    call get_d2bdxjdxj(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(4, 4, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdxj(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(5, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdxj(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(6, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyj(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(4, 5, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyj(xji, zji, yji, xkj, zkj, ykj, derivative)
    call gidx(4, 6, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyj(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(5, 6, ij)
    d2bdxdy(ij) = derivative

    ! i k 

    call get_d2bdxidxk(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxk(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxk(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(3, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(1, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(xji, zji, yji, xkj, zkj, ykj, derivative)
    call gidx(1, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(2, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(yji, xji, zji, ykj, xkj, zkj, derivative)
    call gidx(2, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(3, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyk(zji, yji, xji, zkj, ykj, xkj, derivative)
    call gidx(3, 8, ij)
    d2bdxdy(ij) = derivative

    ! j k
    call get_d2bdxjdxk(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(4, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdxk(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(5, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdxk(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(6, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(xji, yji, zji, xkj, ykj, zkj, derivative)
    call gidx(4, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(xji, zji, yji, xkj, zkj, ykj, derivative)
    call gidx(4, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(yji, zji, xji, ykj, zkj, xkj, derivative)
    call gidx(5, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(yji, xji, zji, ykj, xkj, zkj, derivative)
    call gidx(5, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(zji, xji, yji, zkj, xkj, ykj, derivative)
    call gidx(6, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxjdyk(zji, yji, xji, zkj, ykj, xkj, derivative)
    call gidx(6, 8, ij)
    d2bdxdy(ij) = derivative

    ! k k 

    call get_d2bdxidxi(xkj, ykj, zkj,xji, yji, zji,  derivative)
    call gidx(7, 7, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(ykj, zkj, xkj, yji, zji, xji, derivative)
    call gidx(8, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidxi(zkj, xkj, ykj,zji, xji, yji,  derivative)
    call gidx(9, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xkj, ykj, zkj,xji, yji, zji,  derivative)
    call gidx(7, 8, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(xkj, zkj, ykj, xji, zji, yji,  derivative)
    call gidx(7, 9, ij)
    d2bdxdy(ij) = derivative
    call get_d2bdxidyi(ykj, zkj, xkj, yji, zji, xji,  derivative)
    call gidx(8, 9, ij)
    d2bdxdy(ij) = derivative

end subroutine get_bhessian


subroutine get_d2bdxidxi(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = ((ykj**2 + zkj**2)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 &
                + (yji*zkj - ykj*zji)**2) - (ykj*(xji*ykj - xkj*yji) + zkj*(xji*zkj - xkj*zji))**2)
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)

derivative = enumerator/denominator
end subroutine get_d2bdxidxi

subroutine get_d2bdxidyi(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = - xkj*ykj*((xji*ykj - xkj*yji)**2 + (xji * zkj - xkj * zji)**2 + (yji * zkj - ykj * zji)**2) &
                + (xkj* (xji * ykj - xkj * yji) - zkj * (yji * zkj - ykj * zji))*(ykj*(xji*ykj - xkj*yji) &
                + zkj*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)

derivative = enumerator/denominator
end subroutine get_d2bdxidyi

subroutine get_d2bdxjdxj(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = - ((yji + ykj)*(xji*ykj - xkj*yji) + (zji + zkj)*(xji*zkj - xkj*zji))**2 + ((yji + ykj)**2 &
                + (zji + zkj)**2)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)

denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)

derivative = enumerator/denominator
end subroutine get_d2bdxjdxj


subroutine get_d2bdxjdyj(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = - (xji + xkj)*(yji + ykj)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2) &
            + (- (xji + xkj)*(xji*ykj - xkj*yji) + (zji + zkj)*(yji*zkj - ykj*zji))*(- (yji + ykj)*(xji*ykj - xkj*yji) &
            - (zji + zkj)*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)

derivative = enumerator/denominator
end subroutine get_d2bdxjdyj

subroutine get_d2bdxidxj(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = - (ykj*(yji + ykj) + zkj*(zji + zkj))*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2) &
                + (ykj*(xji*ykj - xkj*yji) + zkj*(xji*zkj - xkj*zji))*((yji + ykj)*(xji*ykj - xkj*yji) + &
                 (zji + zkj)*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)

derivative = enumerator/denominator
end subroutine get_d2bdxidxj

subroutine get_d2bdxidyj(xji, yji, zji, xkj, ykj, zkj, derivative)

real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = (ykj*(xji*ykj - xkj*yji) + zkj*(xji*zkj - xkj*zji))*(- (xji + xkj)*(xji*ykj - xkj*yji) + &
             (zji + zkj)*(yji*zkj - ykj*zji)) - (- xji*ykj + xkj*yji - ykj*(xji + xkj))*((xji*ykj - xkj*yji)**2 &
              + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)
derivative = enumerator/denominator
end subroutine get_d2bdxidyj

subroutine get_d2bdxidxk(xji, yji, zji, xkj, ykj, zkj, derivative)
real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = (yji*ykj + zji*zkj)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2) &
            - (yji*(xji*ykj - xkj*yji) + zji*(xji*zkj - xkj*zji))*(ykj*(xji*ykj - xkj*yji) + zkj*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)
derivative = enumerator/denominator
end subroutine get_d2bdxidxk

subroutine get_d2bdxidyk(xji, yji, zji, xkj, ykj, zkj, derivative)
real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = (- 2*xji*ykj + xkj*yji)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2) &
            + (xji*(xji*ykj - xkj*yji) - zji*(yji*zkj - ykj*zji))*(ykj*(xji*ykj - xkj*yji) + zkj*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)
derivative = enumerator/denominator
end subroutine get_d2bdxidyk

subroutine get_d2bdxjdxk(xji, yji, zji, xkj, ykj, zkj, derivative)
real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = - (yji*(yji + ykj) + zji*(zji + zkj))*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 & 
            + (yji*zkj - ykj*zji)**2) + (yji*(xji*ykj - xkj*yji) + zji*(xji*zkj - xkj*zji))*((yji + ykj)*(xji*ykj - xkj*yji) &
            + (zji + zkj)*(xji*zkj - xkj*zji))
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)
derivative = enumerator/denominator
end subroutine get_d2bdxjdxk

subroutine get_d2bdxjdyk(xji, yji, zji, xkj, ykj, zkj, derivative)
real(8), intent(in) :: xji, yji, zji, xkj, ykj, zkj
real(8), intent(out) :: derivative
real(8) :: enumerator, denominator
enumerator = (xji*(xji*ykj - xkj*yji) - zji*(yji*zkj - ykj*zji))*(- (yji + ykj)*(xji*ykj - xkj*yji) - &
            (zji + zkj)*(xji*zkj - xkj*zji)) - (- xji*ykj - xji*(yji + ykj) &
            + xkj*yji)*((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)
denominator = (((xji*ykj - xkj*yji)**2 + (xji*zkj - xkj*zji)**2 + (yji*zkj - ykj*zji)**2)**1.5_wp)
derivative = enumerator/denominator
end subroutine get_d2bdxjdyk


subroutine get_bgradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, bijk, dbdx)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,bijk
    real(8), intent(out) :: dbdx(12)
    real(8) :: ddx, ddy, ddz

    dbdx = 0.0_wp

    call get_dbdxi(xji,yji,zji,xkj,ykj,zkj,bijk,ddx)
    call get_dbdxi(yji,zji,xji,ykj,zkj,xkj,bijk,ddy)
    call get_dbdxi(zji,xji,yji,zkj,xkj,ykj,bijk,ddz)

    dbdx(1) = ddx
    dbdx(2) = ddy
    dbdx(3) = ddz 


    call get_dbdxj(xji,yji,zji,xkj,ykj,zkj,bijk,ddx)
    call get_dbdxj(yji,zji,xji,ykj,zkj,xkj,bijk,ddy)
    call get_dbdxj(zji,xji,yji,zkj,xkj,ykj,bijk,ddz)

    dbdx(4) = ddx
    dbdx(5) = ddy
    dbdx(6) = ddz 

    call get_dbdxk(xji, yji, zji, xkj, ykj, zkj, bijk, ddx)
    call get_dbdxk(yji, zji, xji, ykj, zkj, xkj, bijk, ddy)
    call get_dbdxk(zji, xji, yji, zkj, xkj, ykj,bijk, ddz)

    dbdx(7) = ddx
    dbdx(8) = ddy
    dbdx(9) = ddz

end subroutine get_bgradient


subroutine get_dbdxi(xji,yji,zji,xkj,ykj,zkj,bijk,derivative)
    real(8), intent(in) :: xji,yji,xkj,ykj,zkj,zji,bijk
    real(8), intent(out) :: derivative
    real(8) :: numerator, denominator
    numerator = ykj*(xji*ykj - xkj*yji) - zkj*(-xji*zkj + xkj*zji)
    denominator = bijk
    derivative = -numerator / denominator
end subroutine get_dbdxi

subroutine get_dbdxj(xji,yji,zji,xkj,ykj,zkj,bijk,derivative)
    real(8), intent(in) :: xji,yji,xkj,ykj,zkj,zji,bijk
    real(8), intent(out) :: derivative
    real(8) :: num1, num2, den
    num1 = yji*(xji*ykj - xkj*yji) - zji*(-xji*zkj + xkj*zji)
    num2 = ykj*(xji*ykj - xkj*yji) - zkj*(-xji*zkj + xkj*zji)
    den = bijk
    derivative = (num1 + num2) / den
end subroutine get_dbdxj

subroutine get_dbdxk(xji, yji, zji, xkj, ykj, zkj, bijk, derivative)
    real(8), intent(in) :: xji, yji, xkj, ykj, zkj, zji, bijk
    real(8), intent(out) :: derivative
    derivative = (-yji*(xji*ykj - xkj*yji) - zji*(xji*zkj - xkj*zji)) / bijk
end subroutine get_dbdxk

end module b_dihedral_derivatives
