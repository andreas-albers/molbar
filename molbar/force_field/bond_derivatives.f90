module bond_derivatives

use fortran_helper, only : get_bij
use b_angle_derivatives, only: get_b_angle_derivatives
implicit none
public :: get_bond_derivatives, get_bond_gradient, get_bond_hessian
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_bond_derivatives(geometry, bonds, ideal_bonds, k_bonds, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_bonds(:)
    integer, intent(in) ::  bonds(:, :)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    integer :: n

    do n=1, size(ideal_bonds)
        call get_single_bond_derivative(geometry, bonds(:,n), ideal_bonds(n), k_bonds, gradient, hessian)
    enddo

end subroutine get_bond_derivatives

subroutine get_bond_gradient(geometry, bonds, ideal_bonds, k_bonds, gradient)
    real(8), intent(in) :: geometry(:, :), ideal_bonds(:)
    integer, intent(in) ::  bonds(:, :)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: gradient(:)
    integer :: n

    do n=1, size(ideal_bonds)
        call get_single_bond_gradient(geometry, bonds(:,n), ideal_bonds(n), k_bonds, gradient)
    enddo

end subroutine get_bond_gradient

subroutine get_bond_hessian(geometry, bonds, ideal_bonds, k_bonds, hessian)
    real(8), intent(in) :: geometry(:, :), ideal_bonds(:)
    integer, intent(in) ::  bonds(:, :)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: hessian(:,:)
    integer :: n

    do n=1, size(ideal_bonds)
        call get_single_bond_hessian(geometry, bonds(:,n), ideal_bonds(n), k_bonds, hessian)
    enddo

end subroutine get_bond_hessian

subroutine get_single_bond_gradient(geometry, bond, r0, k_bonds, gradient)
    real(8), intent(in) :: geometry(:, :), r0
    integer, intent(in) ::  bond(2)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: gradient(:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: drdx(6), d2rdxdy(21)

    i = bond(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = bond(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_r_derivatives(geometry, i, j, drdx, d2rdxdy, rij)

    call build_bond_gradient(r0, rij, drdx, k_bonds, xhat, gradient)

end subroutine get_single_bond_gradient

subroutine get_single_bond_hessian(geometry, bond, r0, k_bonds, hessian)
    real(8), intent(in) :: geometry(:, :), r0
    integer, intent(in) ::  bond(2)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: drdx(6), d2rdxdy(21)

    i = bond(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = bond(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_r_derivatives(geometry, i, j, drdx, d2rdxdy, rij)

    call build_bond_hessian(r0, rij, drdx, d2rdxdy, k_bonds, xhat, hessian)

end subroutine get_single_bond_hessian

subroutine get_single_bond_derivative(geometry, bond, r0, k_bonds, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), r0
    integer, intent(in) ::  bond(2)
    real(8), intent(in) :: k_bonds
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: drdx(6), d2rdxdy(21)

    i = bond(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = bond(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_r_derivatives(geometry, i, j, drdx, d2rdxdy, rij)

    call build_bond_gradient(r0, rij, drdx, k_bonds, xhat, gradient)

    call build_bond_hessian(r0, rij, drdx, d2rdxdy, k_bonds, xhat, hessian)

end subroutine get_single_bond_derivative

subroutine build_bond_gradient(r0, rij, drdx, k_bonds, xhat, gradient)

    real(8), intent(in) :: r0, rij, drdx(6), k_bonds
    integer, intent(in) :: xhat(6)
    real(8), intent(inout) :: gradient(:)
    real(8) :: g
    integer :: n

    do n=1, 6
        call get_dVbonddx(drdx(n), k_bonds, r0, rij, g)
        gradient(xhat(n)) = gradient(xhat(n)) + g
    enddo

end subroutine build_bond_gradient

subroutine build_bond_hessian(r0, rij, drdx, d2rdxdy, k_bonds, xhat, hessian)

    real(8), intent(in) :: r0, rij, drdx(6), d2rdxdy(21)
    real(8), intent(in) :: k_bonds
    integer, intent(in) :: xhat(6)
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: h
    integer :: m, i, j

    do m=1, 21

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        call get_d2Vbonddxdy(drdx(i), drdx(j), d2rdxdy(m), k_bonds, r0, rij, h)

        hessian(xhat(i), xhat(j)) = hessian(xhat(i), xhat(j)) + h

        if (i /= j) then

            hessian(xhat(j), xhat(i)) = hessian(xhat(j), xhat(i)) + h

        endif
    enddo

end subroutine build_bond_hessian

subroutine get_r_derivatives(geometry, i, j, drdx, d2rdxdy, rij)

    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j
    real(8), intent(out) :: drdx(6), d2rdxdy(21)
    real(8), intent(out) :: rij
    real(8) :: xij, yij, zij, dbdx(9),d2bdxdy(45)
    integer :: m

    xij = geometry(1,i) - geometry(1,j)
    yij = geometry(2,i) - geometry(2,j)
    zij = geometry(3,i) - geometry(3,j)

    call get_bij(geometry, i, j, rij)

    call get_b_angle_derivatives(xij,yij,zij,rij,dbdx,d2bdxdy)

    drdx = 0.0_wp
    do m=1, 6
        drdx(m) = dbdx(m)
    enddo

    d2rdxdy = 0.0_wp
    do m=1, 21
        d2rdxdy(m) = d2bdxdy(m)
    enddo

end subroutine get_r_derivatives


subroutine get_dVbonddx(drdx, k_bonds, r0, rij, dVbonddx)
real(8), intent(in) ::drdx, k_bonds, r0, rij
real(8), intent(out) :: dVbonddx

dVbonddx = -2*k_bonds*(r0-rij)*drdx

end subroutine get_dVbonddx

subroutine get_d2Vbonddxdy(drdx, drdy, d2rdxdy, k_bonds, r0, rij, d2Vbonddxdy)
real(8), intent(in) :: drdx, drdy, d2rdxdy, k_bonds, r0, rij
real(8), intent(out) :: d2Vbonddxdy

d2Vbonddxdy = 2*k_bonds*drdx*drdy-2*k_bonds*(r0 - rij)*d2rdxdy
end subroutine get_d2Vbonddxdy

end module bond_derivatives