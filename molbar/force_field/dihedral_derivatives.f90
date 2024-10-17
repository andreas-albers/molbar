module dihedral_derivatives
   
use fortran_helper, only :  get_aijkl, get_bijk, get_cjkl, get_dijkl
use a_dihedral_derivatives, only: get_a_dihedral_derivatives
use b_dihedral_derivatives, only: get_b_dihedral_derivatives
use c_dihedral_derivatives, only: get_c_dihedral_derivatives
use d_dihedral_derivatives, only: get_d_dihedral_derivatives
implicit none
public :: get_dihedral_derivatives, get_dihedral_gradient, get_dihedral_hessian
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_dihedral_derivatives(geometry, dihedrals, ideal_dihedrals, k_dihedrals, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_dihedrals(:)
    integer, intent(in) ::  dihedrals(:, :)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    integer :: n

    do n=1, size(ideal_dihedrals)

        call get_single_dihedral_derivative(geometry, dihedrals(:,n), ideal_dihedrals(n), k_dihedrals, gradient, hessian)
        
    enddo
end subroutine get_dihedral_derivatives

subroutine get_dihedral_gradient(geometry, dihedrals, ideal_dihedrals, k_dihedrals, gradient)
    real(8), intent(in) :: geometry(:, :), ideal_dihedrals(:)
    integer, intent(in) ::  dihedrals(:, :)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: gradient(:)
    integer :: n

    do n=1, size(ideal_dihedrals)

        call get_single_dihedral_gradient(geometry, dihedrals(:,n), ideal_dihedrals(n), k_dihedrals, gradient)
        
    enddo
end subroutine get_dihedral_gradient


subroutine get_dihedral_hessian(geometry, dihedrals, ideal_dihedrals, k_dihedrals, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_dihedrals(:)
    integer, intent(in) ::  dihedrals(:, :)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: hessian(:,:)
    integer :: n

    do n=1, size(ideal_dihedrals)

        call get_single_dihedral_hessian(geometry, dihedrals(:,n), ideal_dihedrals(n), k_dihedrals, hessian)
        
    enddo
end subroutine get_dihedral_hessian


subroutine get_single_dihedral_derivative(geometry, dihedral, phi0, k_dihedrals, gradient, hessian)
        
    real(8), intent(in) :: geometry(:, :), phi0
    integer, intent(in) ::  dihedral(4)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    real(8) :: cosphi0, sinphi0, cosphi, sinphi
    integer :: i,j,k,l, x_hat_i, y_hat_i, z_hat_i, x_hat_j, y_hat_j, z_hat_j
    integer :: xhat(12)
    real(8) :: gxi, gyi, gzi, gxj, gyj, gzj, gxk, gyk, gzk, gxl, gyl, gzl
    real(8) :: dcosdx(12), dsindx(12), d2cosdxdy(78), d2sindxdy(78), aijkl, bijk, cjkl, dijkl
    integer :: m

    i = dihedral(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = dihedral(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = dihedral(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    l = dihedral(4)

    xhat(10) = 1+3*(l-1)
    xhat(11) = 2+3*(l-1)
    xhat(12) = 3+3*(l-1)

    call get_phi_derivatives(geometry, i, j, k, l, dcosdx, d2cosdxdy,dsindx, d2sindxdy, aijkl, bijk, cjkl, dijkl)

    cosphi0 = cos(phi0)
    sinphi0 = sin(phi0)
    cosphi = aijkl/(bijk*cjkl)
    sinphi = dijkl/(bijk*cjkl)


    call build_dihedral_gradient(phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx, dsindx,k_dihedrals,xhat, gradient)

    call build_dihedral_hessian(phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx, d2cosdxdy, &
                                 dsindx, d2sindxdy, k_dihedrals, xhat, hessian)


end subroutine get_single_dihedral_derivative

subroutine get_single_dihedral_gradient(geometry, dihedral, phi0, k_dihedrals, gradient)
        
    real(8), intent(in) :: geometry(:, :), phi0
    integer, intent(in) ::  dihedral(4)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: gradient(:)
    real(8) :: cosphi0, sinphi0, cosphi, sinphi
    integer :: i,j,k,l, x_hat_i, y_hat_i, z_hat_i, x_hat_j, y_hat_j, z_hat_j
    integer :: xhat(12)
    real(8) :: gxi, gyi, gzi, gxj, gyj, gzj, gxk, gyk, gzk, gxl, gyl, gzl
    real(8) :: dcosdx(12), dsindx(12), d2cosdxdy(78), d2sindxdy(78), aijkl, bijk, cjkl, dijkl
    integer :: m

    i = dihedral(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = dihedral(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = dihedral(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    l = dihedral(4)

    xhat(10) = 1+3*(l-1)
    xhat(11) = 2+3*(l-1)
    xhat(12) = 3+3*(l-1)

    call get_phi_derivatives(geometry, i, j, k, l, dcosdx, d2cosdxdy,dsindx, d2sindxdy, aijkl, bijk, cjkl, dijkl)

    cosphi0 = cos(phi0)
    sinphi0 = sin(phi0)
    cosphi = aijkl/(bijk*cjkl)
    sinphi = dijkl/(bijk*cjkl)

    call build_dihedral_gradient(phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx, dsindx,k_dihedrals,xhat, gradient)

end subroutine get_single_dihedral_gradient

subroutine get_single_dihedral_hessian(geometry, dihedral, phi0, k_dihedrals, hessian)
        
    real(8), intent(in) :: geometry(:, :), phi0
    integer, intent(in) ::  dihedral(4)
    real(8), intent(in) :: k_dihedrals
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: cosphi0, sinphi0, cosphi, sinphi
    integer :: i,j,k,l, x_hat_i, y_hat_i, z_hat_i, x_hat_j, y_hat_j, z_hat_j
    integer :: xhat(12)
    real(8) :: gxi, gyi, gzi, gxj, gyj, gzj, gxk, gyk, gzk, gxl, gyl, gzl
    real(8) :: dcosdx(12), dsindx(12), d2cosdxdy(78), d2sindxdy(78), aijkl, bijk, cjkl, dijkl
    integer :: m

    i = dihedral(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = dihedral(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = dihedral(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    l = dihedral(4)

    xhat(10) = 1+3*(l-1)
    xhat(11) = 2+3*(l-1)
    xhat(12) = 3+3*(l-1)

    call get_phi_derivatives(geometry, i, j, k, l, dcosdx, d2cosdxdy,dsindx, d2sindxdy, aijkl, bijk, cjkl, dijkl)

    cosphi0 = cos(phi0)
    sinphi0 = sin(phi0)
    cosphi = aijkl/(bijk*cjkl)
    sinphi = dijkl/(bijk*cjkl)

    call build_dihedral_hessian(phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx, d2cosdxdy, &
                                 dsindx, d2sindxdy, k_dihedrals, xhat, hessian)


end subroutine get_single_dihedral_hessian

subroutine build_dihedral_gradient(phi0, cosphi0,sinphi0,sinphi, cosphi,dcosdx, dsindx, k_dihedrals,xhat, gradient)

    real(8), intent(in) :: phi0, cosphi0, sinphi0, sinphi, cosphi, dcosdx(12), dsindx(12), k_dihedrals
    integer, intent(in) :: xhat(12)
    real(8), intent(inout) :: gradient(:)
    real(8) :: g
    integer :: n

    do n=1, 12

        call get_dVdihedraldx(dcosdx(n), dsindx(n), k_dihedrals, cosphi0, sinphi0, cosphi, sinphi, g)
        gradient(xhat(n)) = gradient(xhat(n)) + g

    enddo
end subroutine build_dihedral_gradient

subroutine build_dihedral_hessian(phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx, d2cosdxdy, dsindx, d2sindxdy, k_dihedrals, xhat, & 
                                hessian)

    real(8), intent(in) :: phi0, cosphi0,sinphi0,sinphi, cosphi, dcosdx(12),d2cosdxdy(78), dsindx(12), d2sindxdy(78)
    real(8), intent(in) :: k_dihedrals
    integer, intent(in) :: xhat(12)
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: h
    integer :: m, i, j

    do m=1, 78

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        call get_d2Vdihedraldxdy(dcosdx(i), dcosdx(j), d2cosdxdy(m), dsindx(i), dsindx(j), d2sindxdy(m), k_dihedrals, &
                                cosphi0, sinphi0, cosphi, sinphi, h)

        hessian(xhat(i), xhat(j)) = hessian(xhat(i), xhat(j)) + h

        if (i /= j) then

            hessian(xhat(j), xhat(i)) = hessian(xhat(j), xhat(i)) + h

        endif
    enddo

end subroutine build_dihedral_hessian

subroutine get_phi_derivatives(geometry, i, j, k, l, dcosdx, d2cosdxdy, dsindx, d2sindxdy, aijkl, bijk, cjkl, dijkl)

    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k, l
    real(8), intent(out) :: dcosdx(12), d2cosdxdy(78), dsindx(12), d2sindxdy(78)
    real(8), intent(out) :: aijkl, bijk, cjkl, dijkl
    real(8) :: dadx(12),dbdx(12),dcdx(12),dddx(12)
    real(8) :: d2adxdy(78),d2bdxdy(78),d2cdxdy(78), d2ddxdy(78)
    real(8) :: xji, yji, zji, xkj, ykj, zkj, xlk, ylk, zlk
    integer :: m

    dadx = 0.0_wp
    dbdx = 0.0_wp
    dcdx = 0.0_wp
    dddx = 0.0_wp

    d2adxdy = 0.0_wp
    d2bdxdy = 0.0_wp
    d2cdxdy = 0.0_wp
    d2ddxdy = 0.0_wp

    xji = geometry(1,j) - geometry(1,i)
    yji = geometry(2,j) - geometry(2,i)
    zji = geometry(3,j) - geometry(3,i)

    xkj = geometry(1,k) - geometry(1,j)
    ykj = geometry(2,k) - geometry(2,j)
    zkj = geometry(3,k) - geometry(3,j)

    xlk = geometry(1,l) - geometry(1,k)
    ylk = geometry(2,l) - geometry(2,k)
    zlk = geometry(3,l) - geometry(3,k)


    call get_aijkl(geometry, i, j, k, l, aijkl)
    call get_bijk(geometry, i, j, k, bijk)
    call get_cjkl(geometry, j, k, l, cjkl)
    call get_dijkl(geometry, i, j, k, l, dijkl)

    call get_a_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,& 
                                dadx, d2adxdy)

    call get_b_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, bijk,& 
                                dbdx, d2bdxdy)

    call get_c_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk, cjkl,& 
                                dcdx, d2cdxdy)

    call get_d_dihedral_derivatives(xji,yji,zji,xkj,ykj,zkj,xlk,ylk,zlk,& 
                                dddx, d2ddxdy)

    call get_phi_gradient(aijkl, bijk, cjkl, dijkl, dadx, dbdx, dcdx, dddx, dcosdx, dsindx)

    call get_phi_hessian(aijkl, bijk, cjkl, dijkl, dadx, dbdx, dcdx, dddx, &
                            d2adxdy, d2bdxdy, d2cdxdy, d2ddxdy, d2cosdxdy, d2sindxdy)


end subroutine get_phi_derivatives

subroutine get_phi_gradient(aijkl, bijk, cjkl, dijkl, dadx, dbdx, dcdx, dddx, dcosdx, dsindx)

    real(8), intent(in) :: aijkl, bijk, cjkl, dijkl, dadx(12), dbdx(12), dcdx(12), dddx(12)
    real(8) :: cosx, sinx
    real(8), intent(out) :: dcosdx(12), dsindx(12)
    integer :: m
    
    do m=1, 12

        if (m < 4) then

        call get_ddxi(aijkl, bijk, cjkl, dadx(m), dbdx(m), cosx)
        call get_ddxi(dijkl, bijk, cjkl, dddx(m), dbdx(m), sinx)
        dcosdx(m) = cosx
        dsindx(m) = sinx

        else if (m > 9) then
        call get_ddxi(aijkl, cjkl, bijk, dadx(m), dcdx(m), cosx)
        call get_ddxi(dijkl, cjkl, bijk, dddx(m), dcdx(m), sinx)
        dcosdx(m) = cosx
        dsindx(m) = sinx

        else

        call get_ddxj(aijkl, bijk, cjkl, dadx(m), dbdx(m), dcdx(m), cosx)
        call get_ddxj(dijkl, bijk, cjkl, dddx(m), dbdx(m), dcdx(m), sinx)
        dcosdx(m) = cosx
        dsindx(m) = sinx
        
        endif

    enddo

end subroutine get_phi_gradient

subroutine get_phi_hessian(aijkl, bijk, cjkl, dijkl, dadx, dbdx, dcdx, dddx, &
                            d2adxdy, d2bdxdy, d2cdxdy, d2ddxdy, d2cosdxdy, d2sindxdy)
    real(8), intent(in) :: aijkl, bijk, cjkl, dijkl, dadx(12), dbdx(12), dcdx(12), dddx(12)
    real(8), intent(in) :: d2adxdy(78), d2bdxdy(78), d2cdxdy(78), d2ddxdy(78)
    real(8), intent(out) :: d2cosdxdy(78), d2sindxdy(78)
    real(8) :: cosxy, sinxy
    integer :: m, i, j
   
    do m=1, 78

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        call get_d2dxdy(aijkl, bijk, cjkl, dadx(i), dadx(j), d2adxdy(m), dbdx(i), dbdx(j), d2bdxdy(m), &
                         dcdx(i), dcdx(j), d2cdxdy(m), cosxy)

        call get_d2dxdy(dijkl, bijk, cjkl, dddx(i), dddx(j), d2ddxdy(m), dbdx(i), dbdx(j), d2bdxdy(m), &
                         dcdx(i), dcdx(j), d2cdxdy(m), sinxy)

        d2cosdxdy(m) = cosxy

        d2sindxdy(m) = sinxy
        
    enddo

end subroutine get_phi_hessian

subroutine get_d2dxdy(aijkl, bijk, cjkl, dadx, dady, d2adxdy, dbdx, dbdy, d2bdxdy, dcdx, dcdy, d2cdxdy, d2dxdy)

    real(8), intent(in) :: aijkl, bijk, cjkl, dadx, dady, d2adxdy, dbdx, dbdy, d2bdxdy, dcdx, dcdy, d2cdxdy
    real(8), intent(out) :: d2dxdy
    real(8) :: term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11

    term1 = aijkl*d2cdxdy/cjkl
    term2 = 2*aijkl*dcdx*dcdy/(cjkl**2)
    term3 = aijkl*d2bdxdy/bijk
    term4 = aijkl*dbdx*dcdy/(bijk*cjkl)
    term5 = aijkl*dbdy*dcdx/(bijk*cjkl)
    term6 = 2*aijkl*dbdx*dbdy/(bijk**2)
    term7 = d2adxdy
    term8 = dadx*dcdy/cjkl
    term9 = dady*dcdx/cjkl
    term10 = dadx*dbdy/bijk
    term11 = dady*dbdx/bijk
    d2dxdy = (-term1+term2-term3+term4+term5+term6+term7-term8-term9-term10-term11)/(bijk*cjkl)

end subroutine get_d2dxdy

subroutine get_dVdihedraldx(dcosdx, dsindx, k_dihedrals, cosphi0, sinphi0, cosphi, sinphi, dVdihedraldx)
real(8), intent(in) :: dcosdx, dsindx, k_dihedrals, cosphi0, sinphi0, cosphi, sinphi
real(8), intent(out) :: dVdihedraldx

dVdihedraldx = -2*k_dihedrals*((cosphi0-cosphi)*dcosdx+(sinphi0-sinphi)*dsindx)

end subroutine get_dVdihedraldx

subroutine get_d2Vdihedraldxdy(dcosdx, dcosdy, d2cosdxdy, dsindx, dsindy, d2sindxdy, k_dihedrals, cosphi0, sinphi0, &
                                cosphi, sinphi, d2Vdihedraldxdy)
real(8), intent(in) :: dcosdx, dcosdy, d2cosdxdy, dsindx, dsindy, d2sindxdy, k_dihedrals, cosphi0, sinphi0, cosphi, sinphi
real(8), intent(out) :: d2Vdihedraldxdy

d2Vdihedraldxdy = 2*k_dihedrals*dcosdx*dcosdy-2*k_dihedrals*(cosphi0 - cosphi)*d2cosdxdy &
                 + 2*k_dihedrals*dsindx*dsindy-2*k_dihedrals*(sinphi0 - sinphi)*d2sindxdy

end subroutine get_d2Vdihedraldxdy

subroutine get_ddxi(aijkl, bijk, cjkl, dadxi, dbdxi, ddxi)
    real(8), intent(in) :: aijkl, cjkl, bijk, dadxi, dbdxi
    real(8), intent(out) :: ddxi

    ddxi = - (aijkl*dbdxi)/(bijk**2*cjkl) + dadxi/(bijk*cjkl)

end subroutine get_ddxi

subroutine get_ddxj(aijkl, bijk, cjkl, dadxj, dbdxj, dcdxj, ddxj)
    real(8), intent(in) :: aijkl, cjkl, bijk, dadxj, dbdxj,dcdxj
    real(8), intent(out) :: ddxj

    ddxj = - (aijkl*dcdxj)/(bijk*cjkl**2)-(aijkl*dbdxj)/(bijk**2*cjkl)+dadxj/(bijk*cjkl)

end subroutine get_ddxj

end module dihedral_derivatives