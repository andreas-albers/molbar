module unsaturated_edges
implicit none
public :: get_unsaturated_atoms, get_unsaturated_bond_matrix, get_unsaturated_bonds
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_unsaturated_atoms(n_atoms, valences, current_valences, unsaturated_atoms, degree_unsaturation)
    integer :: n_atoms
    integer, intent(in) :: valences(:)
    integer, intent(in) :: current_valences(:)
    integer :: i, j, k, l, atomic_valence, current_valence
    integer, intent(out) :: unsaturated_atoms(n_atoms), degree_unsaturation(n_atoms)

    unsaturated_atoms = 0
    degree_unsaturation = 0

    do i=1, n_atoms

        atomic_valence = valences(i)

        current_valence = current_valences(i)

        if (atomic_valence > current_valence) then
                unsaturated_atoms(i) = 1
                degree_unsaturation(i) = atomic_valence - current_valence
        endif
    enddo

end subroutine get_unsaturated_atoms

subroutine get_unsaturated_bond_matrix(n_atoms, cn_matrix, unsaturated_atoms, unsaturated_bond_matrix, n_bonds)
    integer :: n_atoms
    integer, intent(in) :: cn_matrix(:,:)
    integer, intent(in) :: unsaturated_atoms(:)
    integer :: i, j, k, l
    integer, intent(out):: unsaturated_bond_matrix(n_atoms,n_atoms)
    integer, intent(out):: n_bonds

    unsaturated_bond_matrix = 0
    n_bonds = 0
    do i=1, n_atoms
        do j=1, n_atoms
            if (j > i .and. cn_matrix(i,j) > 0 .and. unsaturated_atoms(i) == 1 .and. unsaturated_atoms(j) == 1) then
                n_bonds = n_bonds + 1
                unsaturated_bond_matrix(i,j) = 1
                unsaturated_bond_matrix(j,i) = 1
            endif
        enddo
    enddo

end subroutine get_unsaturated_bond_matrix

subroutine get_unsaturated_bonds(unsaturated_bond_matrix, unsaturated_bonds)
    integer, intent(in) :: unsaturated_bond_matrix(:,:)
    integer, intent(out) :: unsaturated_bonds(:,:)
    integer :: i, j, k, l, counter
    unsaturated_bonds = 0
    k = 1
    do i=1, size(unsaturated_bond_matrix, dim=1)
        do j=1, size(unsaturated_bond_matrix, dim=2)
            if (j > i .and. unsaturated_bond_matrix(i,j) > 0) then
                unsaturated_bonds(1,k) = i 
                unsaturated_bonds(2,k) = j
                k = k + 1
            endif
        enddo
    enddo
end subroutine get_unsaturated_bonds

end module unsaturated_edges