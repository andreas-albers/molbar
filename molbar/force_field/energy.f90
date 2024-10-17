module force_field_energy

use fortran_helper, only :  get_R_ij, get_theta_ijk, get_cosphi_ijkl, get_sinphi_ijkl

implicit none
public :: get_energy

contains

subroutine get_energy(n_atoms, geometry, bonds, ideal_bonds, angles, ideal_angles, &
    dihedrals, ideal_dihedrals, repulsion, charges, k_bonds, k_angles, k_dihedrals, energy)

    intrinsic :: selected_real_kind
    integer, parameter :: wp = selected_real_kind(15)
    real(8), intent(in) :: geometry(:, :), ideal_bonds(:), ideal_angles(:), ideal_dihedrals(:), charges(:)
    integer, intent(in) ::  n_atoms, bonds(:, :), angles(:, :), dihedrals(:,:), repulsion(:,:)
    real(8), intent(in) :: k_bonds, k_angles, k_dihedrals
    real(8), intent(out) :: energy
    real(8) :: Rij, thetaijk, sinphi_ijkl, cosphi_ijkl, E_bond, E_angle, E_dihedral, E_rep
    integer :: n, i, j, k,l

    energy = 0.0_wp
    E_bond = 0.0_wp
    E_angle = 0.0_wp
    E_dihedral = 0.0_wp
    E_rep = 0.0_wp

    do n=1, size(ideal_bonds)

        i = bonds(1, n)
        j = bonds(2, n)
        call get_R_ij(geometry, i, j, Rij)

        E_bond = E_bond + k_bonds*(ideal_bonds(n)-Rij)**2

    enddo

    do n=1, size(ideal_angles)

        i = angles(1, n)
        j = angles(2, n)
        k = angles(3, n)

        call get_theta_ijk(geometry, i, j, k, thetaijk)

        E_angle = E_angle + k_angles*(ideal_angles(n)-thetaijk)**2

    enddo

    do n=1, size(ideal_dihedrals)

        i = dihedrals(1, n)
        j = dihedrals(2, n)
        k = dihedrals(3, n)
        l = dihedrals(4, n)

        call get_cosphi_ijkl(geometry, i, j, k,l, cosphi_ijkl)
        call get_sinphi_ijkl(geometry, i, j, k,l, sinphi_ijkl)

        E_dihedral = E_dihedral + k_dihedrals*((cos(ideal_dihedrals(n))-cosphi_ijkl)**2 &
                + (sin(ideal_dihedrals(n))-sinphi_ijkl)**2)
    enddo

    do n=1, size(charges)

        i = repulsion(1, n)
        j = repulsion(2, n)

        call get_R_ij(geometry, i, j, Rij)

        E_rep = E_rep + charges(n)/(Rij)

    enddo

    energy = E_bond + E_angle + E_dihedral + E_rep

end subroutine get_energy

end module force_field_energy