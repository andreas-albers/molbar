module angle_derivatives

use fortran_helper, only :  get_aijk, get_bij
use a_angle_derivatives, only: get_a_angle_derivatives
use b_angle_derivatives, only: get_b_angle_derivatives
use c_angle_derivatives, only: get_c_angle_derivatives
implicit none
public :: get_angle_derivatives, get_angle_gradient, get_angle_hessian
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_angle_derivatives(geometry, angles, ideal_angles, k_angles, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_angles(:)
    integer, intent(in) ::  angles(:, :)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    integer :: n

    do n=1, size(ideal_angles)
        call get_single_angle_derivative(geometry, angles(:,n), ideal_angles(n), k_angles, gradient, hessian)
    enddo

end subroutine get_angle_derivatives


subroutine get_angle_gradient(geometry, angles, ideal_angles, k_angles, gradient)
    real(8), intent(in) :: geometry(:, :), ideal_angles(:)
    integer, intent(in) ::  angles(:, :)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) :: gradient(:)
    integer :: n

    do n=1, size(ideal_angles)
        call get_single_angle_gradient(geometry, angles(:,n), ideal_angles(n), k_angles, gradient)
    enddo

end subroutine get_angle_gradient

subroutine get_angle_hessian(geometry, angles, ideal_angles, k_angles, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_angles(:)
    integer, intent(in) ::  angles(:, :)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) :: hessian(:,:)
    integer :: n

    do n=1, size(ideal_angles)
        call get_single_angle_hessian(geometry, angles(:,n), ideal_angles(n), k_angles, hessian)
    enddo

end subroutine get_angle_hessian

subroutine get_single_angle_derivative(geometry, angle, theta0, k_angles, gradient, hessian)
        
    real(8), intent(in) :: geometry(:, :), theta0
    integer, intent(in) ::  angle(3)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    real(8) :: costheta0, costheta
    integer :: i,j,k,l,m
    integer :: xhat(9)
    real(8) :: dcosdx(9), d2cosdxdy(45), aijk , bij, cjk

    i = angle(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = angle(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = angle(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    call get_theta_derivatives(geometry, i, j, k, dcosdx, d2cosdxdy, aijk, bij, cjk)

    costheta0 = cos(theta0)

    costheta = aijk/(bij*cjk)

    call build_angle_gradient(costheta0, costheta, dcosdx, k_angles, xhat, gradient)

    call build_angle_hessian(costheta0, costheta, dcosdx, d2cosdxdy, k_angles, xhat, hessian)

end subroutine get_single_angle_derivative

subroutine get_single_angle_gradient(geometry, angle, theta0, k_angles, gradient)
        
    real(8), intent(in) :: geometry(:, :), theta0
    integer, intent(in) ::  angle(3)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) :: gradient(:)
    real(8) :: costheta0, costheta
    integer :: i,j,k,l,m
    integer :: xhat(9)
    real(8) :: dcosdx(9), d2cosdxdy(45), aijk , bij, cjk

    i = angle(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = angle(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = angle(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    call get_theta_derivatives(geometry, i, j, k, dcosdx, d2cosdxdy, aijk, bij, cjk)

    costheta0 = cos(theta0)

    costheta = aijk/(bij*cjk)

    call build_angle_gradient(costheta0, costheta, dcosdx, k_angles, xhat, gradient)

end subroutine get_single_angle_gradient

subroutine get_single_angle_hessian(geometry, angle, theta0, k_angles, hessian)
        
    real(8), intent(in) :: geometry(:, :), theta0
    integer, intent(in) ::  angle(3)
    real(8), intent(in) :: k_angles
    real(8), intent(inout) ::  hessian(:,:)
    real(8) :: costheta0, costheta
    integer :: i,j,k,l,m
    integer :: xhat(9)
    real(8) :: dcosdx(9), d2cosdxdy(45), aijk , bij, cjk

    i = angle(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = angle(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    k = angle(3)

    xhat(7) = 1+3*(k-1)
    xhat(8) = 2+3*(k-1)
    xhat(9)= 3+3*(k-1)

    call get_theta_derivatives(geometry, i, j, k, dcosdx, d2cosdxdy, aijk, bij, cjk)

    costheta0 = cos(theta0)

    costheta = aijk/(bij*cjk)

    call build_angle_hessian(costheta0, costheta, dcosdx, d2cosdxdy, k_angles, xhat, hessian)
    
end subroutine get_single_angle_hessian


subroutine build_angle_gradient(costheta0, costheta, dcosdx, k_angles, xhat, gradient)

    real(8), intent(in) :: costheta0, costheta, dcosdx(9), k_angles
    integer, intent(in) :: xhat(9)
    real(8), intent(inout) :: gradient(:)
    real(8) :: g
    integer :: n

    do n=1, 9

        call get_dVangledx(dcosdx(n), k_angles, costheta0, costheta, g)
        gradient(xhat(n)) = gradient(xhat(n)) + g

    enddo
end subroutine build_angle_gradient

subroutine build_angle_hessian(costheta0, costheta, dcosdx, d2cosdxdy, k_angles, xhat, hessian)

    real(8), intent(in) :: costheta0, costheta, dcosdx(9), d2cosdxdy(45)
    real(8), intent(in) :: k_angles
    integer, intent(in) :: xhat(9)
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: h
    integer :: m, i, j

    do m=1, 45

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        call get_d2Vangledxdy(dcosdx(i), dcosdx(j), d2cosdxdy(m), k_angles, costheta0, costheta, h)

        hessian(xhat(i), xhat(j)) = hessian(xhat(i), xhat(j)) + h

        if (i /= j) then

            hessian(xhat(j), xhat(i)) = hessian(xhat(j), xhat(i)) + h

        endif
    enddo

end subroutine build_angle_hessian

subroutine get_theta_derivatives(geometry, i, j, k, dcosdx, d2cosdxdy, aijk, bij, cjk)

    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k
    real(8), intent(out) :: dcosdx(9), d2cosdxdy(45)
    real(8), intent(out) :: aijk, bij, cjk
    real(8) :: dadx(9),dbdx(9),dcdx(9)
    real(8) :: d2adxdy(45), d2bdxdy(45), d2cdxdy(45)
    real(8) :: xij, yij, zij, xkj, ykj, zkj
    integer :: m

    xij = geometry(1,i) - geometry(1,j)
    yij = geometry(2,i) - geometry(2,j)
    zij = geometry(3,i) - geometry(3,j)

    xkj = geometry(1,k) - geometry(1,j)
    ykj = geometry(2,k) - geometry(2,j)
    zkj = geometry(3,k) - geometry(3,j)

    call get_aijk(geometry, i, j, k, aijk)
    call get_bij(geometry, i, j, bij)
    call get_bij(geometry, j, k, cjk)

    call get_a_angle_derivatives(xij,yij,zij,xkj,ykj,zkj,dadx,d2adxdy)

    call get_b_angle_derivatives(xij,yij,zij,bij,dbdx,d2bdxdy)

    call get_c_angle_derivatives(xkj,ykj,zkj,cjk, dcdx, d2cdxdy)

    call get_theta_gradient(aijk, bij, cjk, dadx, dbdx, dcdx, dcosdx)

    call get_theta_hessian(aijk, bij, cjk, dadx, dbdx, dcdx, &
                            d2adxdy, d2bdxdy, d2cdxdy, d2cosdxdy)

end subroutine get_theta_derivatives

subroutine get_theta_gradient(aijk, bij, cjk, dadx, dbdx, dcdx, dcosdx)

    real(8), intent(in) :: aijk, bij, cjk, dadx(9), dbdx(9), dcdx(9)
    real(8) :: cosx
    real(8), intent(out) :: dcosdx(9)
    integer :: m
    
    dcosdx = 0.0_wp

    do m=1, 9

        if (m < 4) then

        call get_ddxi(aijk, bij, cjk, dadx(m), dbdx(m), cosx)
        dcosdx(m) = cosx

        else if (m > 6) then
        call get_ddxi(aijk, cjk, bij, dadx(m), dcdx(m), cosx)
        dcosdx(m) = cosx

        else
        call get_ddxj(aijk, bij, cjk, dadx(m), dbdx(m), dcdx(m), cosx)
        dcosdx(m) = cosx
        
        endif

    enddo

end subroutine get_theta_gradient

subroutine get_theta_hessian(aijk, bij, cjk, dadx, dbdx, dcdx, d2adxdy, d2bdxdy, d2cdxdy, d2cosdxdy)
  
    real(8), intent(in) :: aijk, bij, cjk, dadx(9), dbdx(9), dcdx(9) 
    real(8), intent(in) :: d2adxdy(45), d2bdxdy(45), d2cdxdy(45)
    real(8), intent(out) :: d2cosdxdy(45)
    real(8) :: cosxy
    integer :: m, i, j
   
    do m=1, 45

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        call get_d2dxdy(aijk, bij, cjk, dadx(i), dadx(j), d2adxdy(m), dbdx(i), dbdx(j), d2bdxdy(m), &
                         dcdx(i), dcdx(j), d2cdxdy(m), cosxy)

        d2cosdxdy(m) = cosxy  
    enddo

end subroutine get_theta_hessian

subroutine get_d2dxdy(aijk, bij, cjk, dadx, dady, d2adxdy, dbdx, dbdy, d2bdxdy, dcdx, dcdy, d2cdxdy, d2dxdy)

    real(8), intent(in) :: aijk, bij, cjk, dadx, dady, d2adxdy, dbdx, dbdy, d2bdxdy, dcdx, dcdy, d2cdxdy
    real(8), intent(out) :: d2dxdy
    real(8) :: term1, term2, term3, term4, term5, term6, term7, term8, term9, term10, term11

    term1 = aijk*d2cdxdy/cjk
    term2 = 2*aijk*dcdx*dcdy/(cjk**2)
    term3 = aijk*d2bdxdy/bij
    term4 = aijk*dbdx*dcdy/(bij*cjk)
    term5 = aijk*dbdy*dcdx/(bij*cjk)
    term6 = 2*aijk*dbdx*dbdy/(bij**2)
    term7 = d2adxdy
    term8 = dadx*dcdy/cjk
    term9 = dady*dcdx/cjk
    term10 = dadx*dbdy/bij
    term11 = dady*dbdx/bij
    d2dxdy = (-term1+term2-term3+term4+term5+term6+term7-term8-term9-term10-term11)/(bij*cjk)

end subroutine get_d2dxdy

subroutine get_dVangledx(dcosdx, k_angles, costheta0, costheta, dVangledx)

real(8), intent(in) :: dcosdx, k_angles, costheta0, costheta
real(8), intent(out) :: dVangledx

dVangledx = -2*k_angles*(dacos(costheta0)-dacos(costheta))*(-dcosdx/(sqrt(1-costheta**2)))

end subroutine get_dVangledx

subroutine get_d2Vangledxdy(dcosdx, dcosdy, d2cosdxdy, k_angles, costheta0, costheta, d2Vangledxdy)
real(8), intent(in) :: dcosdx, dcosdy, d2cosdxdy, k_angles, costheta0, costheta
real(8), intent(out) :: d2Vangledxdy
real(8) :: dthetadcos, d2thetadcos2

dthetadcos = -1/sqrt(1-costheta**2)

d2thetadcos2 = - costheta/sqrt(1-costheta**2)**1.5_wp

d2Vangledxdy = 2*k_angles*(-(dacos(costheta0)-dacos(costheta))*dthetadcos*d2cosdxdy &
            -(dacos(costheta0)-dacos(costheta))*d2thetadcos2*dcosdx*dcosdy+dthetadcos**2*dcosdx*dcosdy)

end subroutine get_d2Vangledxdy

subroutine get_ddxi(aijk, bij, cjk, dadxi, dbdxi, ddxi)

    real(8), intent(in) :: aijk, bij, cjk, dadxi, dbdxi
    real(8), intent(out) :: ddxi

    ddxi = - (aijk*dbdxi)/(bij**2*cjk) + dadxi/(bij*cjk)

end subroutine get_ddxi

subroutine get_ddxj(aijk, bij, cjk, dadxj, dbdxj, dcdxj, ddxj)

    real(8), intent(in) :: aijk, cjk, bij, dadxj, dbdxj,dcdxj
    real(8), intent(out) :: ddxj

    ddxj = - (aijk*dcdxj)/(bij*cjk**2)-(aijk*dbdxj)/(bij**2*cjk)+dadxj/(bij*cjk)

end subroutine get_ddxj

end module angle_derivatives