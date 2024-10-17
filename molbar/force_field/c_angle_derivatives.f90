module c_angle_derivatives
use b_angle_derivatives, only: get_b_angle_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_c_angle_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_c_angle_derivatives(xkj,ykj,zkj,cjk, dcdx, d2cdxdy)
    intrinsic :: floor, sqrt
    real(8), intent(in) :: xkj,ykj,zkj,cjk
    real(8), intent(out) :: dcdx(9), d2cdxdy(45)
    real(8) :: dbdx(9), d2bdxdy(45)
    integer :: n, m, i, j, ij

    dcdx = 0.0_wp
    d2cdxdy = 0.0_wp

    call get_b_angle_derivatives(xkj,ykj,zkj,cjk,dbdx,d2bdxdy)
    
    do m = 1, 3
        dcdx(m+6) = dbdx(m)
    enddo
    do m = 4, 6
        dcdx(m) = dbdx(m)
    enddo

    do m = 1, 21
        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        if (i < 4) then
            i = i + 6
        endif
        if (j < 4) then
            j = j + 6
        endif
        call gidx(i, j, ij)
        d2cdxdy(ij) = d2bdxdy(m)
    enddo

end subroutine get_c_angle_derivatives

end module c_angle_derivatives
