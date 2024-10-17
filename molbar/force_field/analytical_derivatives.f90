module derivatives

use bond_derivatives, only: get_bond_derivatives, get_bond_gradient, get_bond_hessian
use angle_derivatives, only: get_angle_derivatives, get_angle_gradient, get_angle_hessian
use dihedral_derivatives, only: get_dihedral_derivatives, get_dihedral_gradient, get_dihedral_hessian
use repulsion_derivatives, only: get_repulsion_derivatives, get_repulsion_gradient, get_repulsion_hessian
implicit none
public :: get_derivatives
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_derivatives(n_atoms, geometry, bonds, ideal_bonds, angles, ideal_angles, &
    dihedrals, ideal_dihedrals, repulsive_partner, charges, k_bonds, k_angles, k_dihedrals, gradient, hessian)
        real(8), intent(in) :: geometry(:, :), ideal_bonds(:), ideal_angles(:), ideal_dihedrals(:), charges(:)
        integer, intent(in) ::  n_atoms, bonds(:, :), angles(:, :), dihedrals(:,:), repulsive_partner(:,:)
        real(8), intent(in) :: k_bonds, k_angles, k_dihedrals
        real(8), intent(out) :: gradient(n_atoms*3), hessian(n_atoms*3, n_atoms*3)

        gradient = 0.0_wp
        hessian = 0.0_wp

        call get_bond_derivatives(geometry, bonds, ideal_bonds, k_bonds, gradient, hessian)

        call get_angle_derivatives(geometry, angles, ideal_angles, k_angles, gradient, hessian)

        call get_dihedral_derivatives(geometry, dihedrals, ideal_dihedrals, k_dihedrals, gradient, hessian)

        call get_repulsion_derivatives(geometry, repulsive_partner, charges, gradient, hessian)

end subroutine


subroutine get_gradient(n_atoms, geometry, bonds, ideal_bonds, angles, ideal_angles, &
    dihedrals, ideal_dihedrals, repulsive_partner, charges, k_bonds, k_angles, k_dihedrals, gradient)
        real(8), intent(in) :: geometry(:, :), ideal_bonds(:), ideal_angles(:), ideal_dihedrals(:), charges(:)
        integer, intent(in) ::  n_atoms, bonds(:, :), angles(:, :), dihedrals(:,:), repulsive_partner(:,:)
        real(8), intent(in) :: k_bonds, k_angles, k_dihedrals
        real(8), intent(out) :: gradient(n_atoms*3)

        gradient = 0.0_wp

        call get_bond_gradient(geometry, bonds, ideal_bonds, k_bonds, gradient)

        call get_angle_gradient(geometry, angles, ideal_angles, k_angles, gradient)

        call get_dihedral_gradient(geometry, dihedrals, ideal_dihedrals, k_dihedrals, gradient)

        call get_repulsion_gradient(geometry, repulsive_partner, charges, gradient)

end subroutine

subroutine get_hessian(n_atoms, geometry, bonds, ideal_bonds, angles, ideal_angles, &
    dihedrals, ideal_dihedrals, repulsive_partner, charges, k_bonds, k_angles, k_dihedrals, hessian)
        real(8), intent(in) :: geometry(:, :), ideal_bonds(:), ideal_angles(:), ideal_dihedrals(:), charges(:)
        integer, intent(in) ::  n_atoms, bonds(:, :), angles(:, :), dihedrals(:,:), repulsive_partner(:,:)
        real(8), intent(in) :: k_bonds, k_angles, k_dihedrals
        real(8), intent(out) :: hessian(n_atoms*3, n_atoms*3)

        hessian = 0.0_wp

        call get_bond_hessian(geometry, bonds, ideal_bonds, k_bonds, hessian)

        call get_angle_hessian(geometry, angles, ideal_angles, k_angles, hessian)

        call get_dihedral_hessian(geometry, dihedrals, ideal_dihedrals, k_dihedrals, hessian)

        call get_repulsion_hessian(geometry, repulsive_partner, charges, hessian)

end subroutine

end module derivatives


