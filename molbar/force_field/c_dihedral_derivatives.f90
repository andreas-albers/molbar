module c_dihedral_derivatives
use b_dihedral_derivatives, only: get_b_dihedral_derivatives
use fortran_helper, only: gidx
implicit none
public :: get_c_dihedral_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_c_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, cjkl, & 
                                dcdx, d2cdxdy)
    intrinsic :: floor, sqrt
    real(8), intent(in) :: xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,cjkl
    real(8), intent(out) :: dcdx(12), d2cdxdy(78)
    real(8) :: ddx, ddy, ddz, dbdx(12), d2bdxdy(78)
    integer :: n, m, i, j, ij

    dcdx = 0.0_wp
    d2cdxdy = 0.0_wp

    call get_b_dihedral_derivatives(xkj,ykj,zkj,xlk,ylk,zlk,xji,yji,zji, cjkl,& 
                                dbdx, d2bdxdy)
    
    do m = 1, 9
        dcdx(m+3) = dbdx(m)
    enddo

    do m = 1, 45
        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2
        i = i + 3
        j = j + 3
        call gidx(i, j, ij)
        d2cdxdy(ij) = d2bdxdy(m)
    enddo

end subroutine get_c_dihedral_derivatives

end module c_dihedral_derivatives
