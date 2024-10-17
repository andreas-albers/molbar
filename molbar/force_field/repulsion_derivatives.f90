module repulsion_derivatives

use fortran_helper, only : get_bij
use inv_rij_derivatives, only: get_inverse_distance_derivatives
implicit none
public :: get_repulsion_derivatives, get_repulsion_gradient, get_repulsion_hessian
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_repulsion_derivatives(geometry, repulsive_partner, charges, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), charges(:)
    integer, intent(in) ::  repulsive_partner(:, :)
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    integer :: n

    do n=1, size(charges)
        call get_single_repulsion_derivative(geometry, repulsive_partner(:,n), charges(n), gradient, hessian)
    enddo

end subroutine get_repulsion_derivatives

subroutine get_repulsion_gradient(geometry, repulsive_partner, charges, gradient)
    real(8), intent(in) :: geometry(:, :), charges(:)
    integer, intent(in) ::  repulsive_partner(:, :)
    real(8), intent(inout) :: gradient(:)
    integer :: n

    do n=1, size(charges)
        call get_single_repulsion_gradient(geometry, repulsive_partner(:,n), charges(n), gradient)
    enddo

end subroutine get_repulsion_gradient

subroutine get_repulsion_hessian(geometry, repulsive_partner, charges, hessian)
    real(8), intent(in) :: geometry(:, :), charges(:)
    integer, intent(in) ::  repulsive_partner(:, :)
    real(8), intent(inout) :: hessian(:,:)
    integer :: n

    do n=1, size(charges)
        call get_single_repulsion_hessian(geometry, repulsive_partner(:,n), charges(n), hessian(:,:))
    enddo

end subroutine get_repulsion_hessian

subroutine get_single_repulsion_gradient(geometry, partner, charge, gradient)
    real(8), intent(in) :: geometry(:, :), charge
    integer, intent(in) ::  partner(2)
    real(8), intent(inout) :: gradient(:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: dqdx(6), d2qdxdy(21)

    i = partner(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = partner(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_inv_dist_derivatives(geometry, i, j, dqdx, d2qdxdy, rij)

    call build_repulsion_gradient(rij, dqdx, charge, xhat, gradient)

end subroutine get_single_repulsion_gradient

subroutine get_single_repulsion_hessian(geometry, partner, charge, hessian)
    real(8), intent(in) :: geometry(:, :), charge
    integer, intent(in) ::  partner(2)
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: dqdx(6), d2qdxdy(21)

    i = partner(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = partner(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_inv_dist_derivatives(geometry, i, j, dqdx, d2qdxdy, rij)

    call build_repulsion_hessian(rij, d2qdxdy, charge, xhat, hessian)

end subroutine get_single_repulsion_hessian

subroutine get_single_repulsion_derivative(geometry, partner, charge, gradient, hessian)
    real(8), intent(in) :: geometry(:, :), charge
    integer, intent(in) ::  partner(2)
    real(8), intent(inout) :: gradient(:), hessian(:,:)
    real(8) :: rij
    integer :: i,j
    integer :: xhat(6)
    real(8) :: dqdx(6), d2qdxdy(21)

    i = partner(1)

    xhat(1) = 1+3*(i-1)
    xhat(2) = 2+3*(i-1)
    xhat(3) = 3+3*(i-1)

    j = partner(2)

    xhat(4) = 1+3*(j-1)
    xhat(5) = 2+3*(j-1)
    xhat(6) = 3+3*(j-1)

    call get_inv_dist_derivatives(geometry, i, j, dqdx, d2qdxdy, rij)

    call build_repulsion_gradient(rij, dqdx, charge, xhat, gradient)

    call build_repulsion_hessian(rij, d2qdxdy, charge, xhat, hessian)

end subroutine get_single_repulsion_derivative

subroutine build_repulsion_gradient(rij, dqdx, charge, xhat, gradient)

    real(8), intent(in) :: rij, dqdx(6), charge
    integer, intent(in) :: xhat(6)
    real(8), intent(inout) :: gradient(:)
    real(8) :: g
    integer :: n

    do n=1, 6
        gradient(xhat(n)) = gradient(xhat(n)) + charge*dqdx(n)
    enddo

end subroutine build_repulsion_gradient

subroutine build_repulsion_hessian(rij, d2qdxdy, charge, xhat, hessian)

    real(8), intent(in) :: rij,d2qdxdy(21)
    real(8), intent(in) :: charge
    integer, intent(in) :: xhat(6)
    real(8), intent(inout) :: hessian(:,:)
    real(8) :: h
    integer :: m, i, j

    do m=1, 21

        j = floor((sqrt(8*real(m) - 7) + 1) / 2)
        i = m - (j-1)*j/2

        hessian(xhat(i), xhat(j)) = hessian(xhat(i), xhat(j)) + charge*d2qdxdy(m)

        if (i /= j) then

            hessian(xhat(j), xhat(i)) = hessian(xhat(j), xhat(i)) + charge*d2qdxdy(m)
        endif
    enddo

end subroutine build_repulsion_hessian

subroutine get_inv_dist_derivatives(geometry, i, j, dqdx, d2qdxdy, rij)
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j
    real(8), intent(out) :: dqdx(6), d2qdxdy(21)
    real(8), intent(out) :: rij
    real(8) :: xij, yij, zij, dbdx(9)
    integer :: m

    xij = geometry(1,i) - geometry(1,j)
    yij = geometry(2,i) - geometry(2,j)
    zij = geometry(3,i) - geometry(3,j)

    call get_bij(geometry, i, j, rij)

    call get_inverse_distance_derivatives(xij,yij,zij,rij,dqdx,d2qdxdy)

end subroutine get_inv_dist_derivatives


end module repulsion_derivatives