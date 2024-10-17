module a_dihedral_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_a_dihedral_derivatives
integer, parameter :: wp = selected_real_kind(15)

contains

subroutine get_a_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,dadx,d2adxdy)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
    real(8), intent(out) :: dadx(12), d2adxdy(78)
    real(8) :: ddx, ddy, ddz
    integer :: i

    call get_agradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,dadx)

    call get_ahessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2adxdy)
    
end subroutine get_a_dihedral_derivatives


subroutine get_ahessian(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,d2adxdy)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
    real(8), intent(out) :: d2adxdy(78)
    real(8) :: derivative
    integer :: ij

    d2adxdy = 0.0_wp
    ! i i = 0
    !j j
    call get_d2adxjdxj(yji, zji, ykj, zkj, ylk, zlk, derivative)
    call gidx(4, 4, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxj(xji, zji, xkj, zkj, xlk, zlk, derivative)
    call gidx(5, 5, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxj(xji, yji, xkj, ykj, xlk, ylk, derivative)
    call gidx(6, 6, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(xji, yji, xkj, ykj, xlk, ylk, derivative)
    call gidx(4, 5, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(xji, zji, xkj, zkj, xlk, zlk, derivative)
    call gidx(4, 6, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(yji, zji, ykj, zkj, ylk, zlk, derivative)
    call gidx(5, 6, ij)
    d2adxdy(ij) = derivative

    !i j
    call get_d2adxidxj(ykj, zkj, ylk, zlk, derivative)
    call gidx(1, 4, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxj(xkj, zkj, xlk, zlk, derivative)
    call gidx(2, 5, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxj(xkj, ykj, xlk, ylk, derivative)
    call gidx(3, 6, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xkj, ykj, xlk, ylk, derivative)
    call gidx(1, 5, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xkj, zkj, xlk, zlk, derivative)
    call gidx(1, 6, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(ykj, zkj, ylk, zlk, derivative)
    call gidx(2, 6, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xlk, ylk, xkj, ykj, derivative)
    call gidx(2, 4, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xlk, zlk, xkj, zkj, derivative)
    call gidx(3, 4, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(ylk, zlk, ykj, zkj, derivative)
    call gidx(3, 5, ij)
    d2adxdy(ij) = derivative

    ! i k 

    call get_d2adxidxk(ykj, zkj, ylk, zlk, derivative)
    call gidx(1, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxk(zkj, xkj, zlk, xlk, derivative)
    call gidx(2, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxk(xkj, ykj, xlk, ylk, derivative)
    call gidx(3, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(xkj, ykj, xlk, ylk,derivative)
    call gidx(1, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(xkj, zkj, xlk, zlk, derivative)
    call gidx(1, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(ykj, zkj, ylk, zlk, derivative)
    call gidx(2, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(ykj, xkj, ylk, xlk, derivative)
    call gidx(2, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(zkj, xkj, zlk, xlk, derivative)
    call gidx(3, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(zkj, ykj, zlk, ylk, derivative)
    call gidx(3, 8, ij)
    d2adxdy(ij) = derivative

    ! j k
    call get_d2adxjdxk(yji,zji,ykj,zkj,ylk,zlk,derivative)
    call gidx(4, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxk(zji,xji,zkj,xkj,zlk,xlk,derivative)
    call gidx(5, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxk(xji,yji,xkj,ykj,xlk,ylk,derivative)
    call gidx(6, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(xji,yji,xkj,ykj,xlk,ylk,derivative)
    call gidx(4, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(xji,zji,xkj,zkj,xlk,zlk,derivative)
    call gidx(4, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(yji,zji,ykj,zkj,ylk,zlk,derivative)
    call gidx(5, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(yji,xji,ykj,xkj,ylk,xlk,derivative)
    call gidx(5, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(zji,xji,zkj,xkj,zlk,xlk,derivative)
    call gidx(6, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyk(zji,yji,zkj,ykj,zlk,ylk,derivative)
    call gidx(6, 8, ij)
    d2adxdy(ij) = derivative

    ! k k 

    call get_d2adxjdxj(ylk, zlk, ykj, zkj, yji, zji, derivative)
    call gidx(7, 7, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxj(xlk, zlk, xkj, zkj, xji, zji, derivative)
    call gidx(8, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdxj(xlk, ylk, xkj, ykj, xji, yji, derivative)
    call gidx(9, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(xlk, ylk, xkj, ykj, xji, yji, derivative)
    call gidx(7, 8, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(xlk, zlk, xkj, zkj, xji, zji, derivative)
    call gidx(7, 9, ij)
    d2adxdy(ij) = derivative
    call get_d2adxjdyj(ylk, zlk, ykj, zkj, yji, zji, derivative)
    call gidx(8, 9, ij)
    d2adxdy(ij) = derivative

    ! i l

    call get_d2adxidxl(ykj,zkj,derivative)
    call gidx(1, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxl(zkj,xkj,derivative)
    call gidx(2, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxl(xkj,ykj,derivative)
    call gidx(3, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(xkj,ykj,derivative)
    call gidx(1, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(xkj,zkj,derivative)
    call gidx(1, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(zkj,ykj,derivative)
    call gidx(2, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(ykj,xkj,derivative)
    call gidx(2, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(zkj,xkj,derivative)
    call gidx(3, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyl(ykj,zkj,derivative)
    call gidx(3, 11, ij)
    d2adxdy(ij) = derivative
    ! j l
    call get_d2adxidxk(ykj, zkj, yji, zji, derivative)
    call gidx(4, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxk(zkj, xkj, zji, xji, derivative)
    call gidx(5, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxk(xkj, ykj, xji, yji, derivative)
    call gidx(6, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(ykj, xkj, yji, xji, derivative)
    call gidx(4, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(zkj, xkj, zji, xji, derivative)
    call gidx(4, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(zkj, ykj, zji, yji, derivative)
    call gidx(5, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(xkj, ykj, xji, yji,derivative)
    call gidx(5, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(xkj, zkj, xji, zji, derivative)
    call gidx(6, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyk(ykj, zkj, yji, zji, derivative)
    call gidx(6, 11, ij)
    d2adxdy(ij) = derivative

    ! k l
    call get_d2adxidxj(ykj, zkj, yji, zji, derivative)
    call gidx(7, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxj(xkj, zkj, xji, zji, derivative)
    call gidx(8, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidxj(xkj, ykj, xji, yji, derivative)
    call gidx(9, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xji, yji, xkj, ykj, derivative)
    call gidx(7, 11, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(xji, zji,xkj, zkj, derivative)
    call gidx(7, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(yji, zji, ykj, zkj, derivative)
    call gidx(8, 12, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(yji, xji, ykj, xkj,derivative)
    call gidx(8, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(zji, xji, zkj, xkj, derivative)
    call gidx(9, 10, ij)
    d2adxdy(ij) = derivative
    call get_d2adxidyj(zji, yji, zkj, ykj, derivative)
    call gidx(9, 11, ij)
    d2adxdy(ij) = derivative 

    ! l l = 0            
    
end subroutine get_ahessian


subroutine get_d2adxjdxj(yji, zji, ykj, zkj, ylk, zlk, derivative)
real(8), intent(in) :: yji, zji, ykj, zkj, ylk, zlk
real(8), intent(out) :: derivative
derivative = -2*ylk*(yji+ykj)-2*zlk*(zji+zkj)
end subroutine get_d2adxjdxj

subroutine get_d2adxjdyj(xji, yji, xkj, ykj, xlk, ylk, derivative)
real(8), intent(in) :: xji, yji, xkj, ykj, xlk, ylk
real(8), intent(out) :: derivative
derivative = xlk*(yji+ykj)+ylk*(xji+xkj)
end subroutine get_d2adxjdyj

subroutine get_d2adxidxj(ykj, zkj, ylk, zlk, derivative)
real(8), intent(in) :: ykj, zkj, ylk, zlk
real(8), intent(out) :: derivative
derivative = ykj*ylk+zkj*zlk
end subroutine get_d2adxidxj

subroutine get_d2adxidyj(xkj, ykj, xlk, ylk, derivative)
real(8), intent(in) :: xkj, ykj, xlk, ylk
real(8), intent(out) :: derivative
derivative = xkj*ylk-2*xlk*ykj
end subroutine get_d2adxidyj

subroutine get_d2adxidxk(ykj, zkj, ylk, zlk, derivative)
real(8), intent(in) :: ykj, zkj, ylk, zlk
real(8), intent(out) :: derivative
derivative = -ykj*(ykj+ylk)-zkj*(zkj+zlk)
end subroutine get_d2adxidxk

subroutine get_d2adxidyk(xkj, ykj, xlk, ylk,derivative)
real(8), intent(in) :: xkj, ykj, xlk, ylk
real(8), intent(out) :: derivative
derivative = -xkj*ylk+xlk*ykj+ykj*(xkj+xlk)
end subroutine get_d2adxidyk

subroutine get_d2adxidxl(ykj,zkj,derivative)
real(8), intent(in) :: ykj,zkj
real(8), intent(out) :: derivative
derivative = ykj**2+zkj**2
end subroutine get_d2adxidxl

subroutine get_d2adxidyl(xkj,ykj,derivative)
real(8), intent(in) :: xkj,ykj
real(8), intent(out) :: derivative
derivative = -xkj*ykj
end subroutine get_d2adxidyl

subroutine get_d2adxjdxk(yji,zji,ykj,zkj,ylk,zlk,derivative)
real(8), intent(in) :: yji,zji,ykj,zkj,ylk,zlk
real(8), intent(out) :: derivative
derivative = yji*ylk+zji*zlk+(yji+ykj)*(ykj+ylk)+(zji+zkj)*(zkj+zlk)
end subroutine get_d2adxjdxk

subroutine get_d2adxjdyk(xji,yji,xkj,ykj,xlk,ylk,derivative)
real(8), intent(in) :: xji,yji,xkj,ykj,xlk,ylk
real(8), intent(out) :: derivative
derivative = xji*ykj-xji*ylk-xkj*yji+xkj*ylk-xlk*ykj-(xkj+xlk)*(yji+ykj)
end subroutine get_d2adxjdyk


subroutine get_agradient(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,dadx)
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk
    real(8), intent(out) :: dadx(12)
    real(8) :: ddx, ddy, ddz

    dadx = 0.0_wp

    call get_dadxi(xkj,ykj,zkj,xlk,ylk,zlk,ddx)
    call get_dadxi(ykj,zkj,xkj,ylk,zlk,xlk,ddy)
    call get_dadxi(zkj,xkj,ykj,zlk,xlk,ylk,ddz)

    dadx(1) = ddx
    dadx(2) = ddy
    dadx(3) = ddz 
    
    call get_dadxj(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,ddx)
    call get_dadxj(yji,zji,xji,ykj,zkj,xkj,ylk,zlk,xlk,ddy)
    call get_dadxj(zji,xji,yji,zkj,xkj,ykj,zlk,xlk,ylk,ddz)


    dadx(4) = ddx
    dadx(5) = ddy
    dadx(6) = ddz 
    
    call get_dadxk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,ddx)
    call get_dadxk(yji,zji,xji,ykj,zkj,xkj,ylk,zlk,xlk,ddy)
    call get_dadxk(zji,xji,yji,zkj,xkj,ykj,zlk,xlk,ylk,ddz)

    dadx(7) = ddx
    dadx(8) = ddy
    dadx(9) = ddz

    call get_dadxl(xkj,ykj,zkj,xji,yji,zji,xlk,ylk,zlk,ddx)
    call get_dadxl(ykj,zkj,xkj,yji,zji,xji,ylk,zlk,xlk,ddy)
    call get_dadxl(zkj,xkj,ykj,zji,xji,yji,zlk,xlk,ylk,ddz)

    dadx(10) = ddx
    dadx(11) = ddy
    dadx(12) = ddz

end subroutine get_agradient


subroutine get_dadxi(xkj,ykj,zkj,xlk,ylk,zlk,derivative)
real(8), intent(in) :: xkj,ykj,zkj,xlk,ylk,zlk
real(8), intent(out) :: derivative
derivative = -ykj*(xkj*ylk - xlk*ykj) + zkj*(-xkj*zlk + xlk*zkj)
end subroutine get_dadxi

subroutine get_dadxj(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
real(8), intent(in) :: xji,yji,xkj,ykj,zkj,zji,xlk,ylk,zlk
real(8), intent(out) :: derivative
real(8) :: term1, term2, term3, term4, term5, term6
    term1= yji*(xkj*ylk-xlk*ykj)
    term2 = ykj*(xkj*ylk-xlk*ykj)
    term3 = ylk*(xji*ykj-xkj*yji)
    term4 = zji*(xkj*zlk-xlk*zkj)
    term5 = zkj*(xkj*zlk-xlk*zkj)
    term6 = zlk*(xji*zkj-xkj*zji)
    derivative = term1+term2-term3+term4+term5-term6
end subroutine get_dadxj

subroutine get_dadxk(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,derivative)
real(8), intent(in) :: xji,yji,xkj,ykj,zkj,zji,xlk,ylk,zlk
real(8), intent(out) :: derivative
real(8) :: term1, term2, term3, term4, term5, term6
    term1 = yji *(xkj *ylk - xlk *ykj)
    term2 = ykj *(xji* ykj - xkj* yji)
    term3 = ylk *(xji *ykj - xkj* yji)
    term4 = zji* (xkj* zlk - xlk* zkj)
    term5 = zkj* (xji* zkj - xkj* zji)
    term6 = zlk *(xji* zkj - xkj* zji)
    derivative = -term1+term2+term3-term4+term5+term6
end subroutine get_dadxk

subroutine get_dadxl(xkj,ykj,zkj,xji,yji,zji,xlk,ylk,zlk,derivative)
    real(8), intent(in) :: xkj,ykj,zkj,xji,yji,zji,xlk,ylk,zlk
    real(8), intent(out) :: derivative
    derivative = -ykj*(xji*ykj - xkj*yji) + zkj*(-xji*zkj + xkj*zji)
end subroutine get_dadxl


end module a_dihedral_derivatives
