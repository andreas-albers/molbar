module bond_order_detection

use unsaturated_edges, only: get_unsaturated_atoms, get_unsaturated_bond_matrix, get_unsaturated_bonds
use valence_conditions, only: get_total_charge, check_valence
use unsaturated_graph, only: get_combinations, get_matching
implicit none
public :: get_best_bo_matrix
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_best_bo_matrix(n_atoms, cn_matrix, cn, atoms, all_possible_valences, are_core, best_bo_matrix)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: cn_matrix(:,:)
    integer, intent(in) :: cn(:)
    integer, intent(in) :: all_possible_valences(:,:)
    integer, intent(in) :: are_core(:)
    integer, intent(in) :: atoms(:)
    integer, intent(out) :: best_bo_matrix(n_atoms,n_atoms)
    integer :: bo_matrix(n_atoms,n_atoms)
    integer :: bo_charge, energy
    integer :: best_bo_charge, best_energy
    integer :: i, j, k, x, y

    best_bo_charge = 0
    best_energy = 0

    do i=1, size(all_possible_valences, dim=2)
        !write(*,*) all_possible_valences(:,i)
        call get_bo_matrix_for_valence(n_atoms,  atoms, cn_matrix, cn, & 
                                        all_possible_valences(:,i), are_core, bo_matrix, bo_charge, energy)

        !write(*,*) "new_best", energy, best_energy, bo_charge, best_bo_charge
        if (i == 1) then
            best_bo_matrix = bo_matrix
            best_bo_charge = bo_charge
            best_energy = energy
        else
            if (energy > best_energy .and. bo_charge <= best_bo_charge) then
                best_bo_charge = bo_charge
                best_energy = energy
                best_bo_matrix = bo_matrix
            else if (energy <= best_energy .and. bo_charge < best_bo_charge) then
                best_bo_charge = bo_charge
                best_energy = energy
                best_bo_matrix = bo_matrix
            else if (energy == best_energy .and. bo_charge == best_bo_charge) then
                do j=1, n_atoms
                    do k=1, n_atoms
                        if (k > j .and. bo_matrix(j,k) > best_bo_matrix(j,k)) then
                            best_bo_matrix(j,k) = bo_matrix(j,k)
                            best_bo_matrix(k,j) = bo_matrix(k,j)
                        endif
                    enddo
                enddo
            endif
        endif
    enddo

end subroutine get_best_bo_matrix

subroutine get_bo_matrix_for_valence(n_atoms, atoms, cn_matrix, cn, valences,&
                                         are_core, best_bo_matrix, best_bo_charge, best_energy)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: cn_matrix(:,:)
    integer, intent(in) :: cn(:)
    integer, intent(in) :: valences(:)
    integer, intent(in) :: atoms(:)
    integer, intent(in) :: are_core(:)
    integer, intent(out) :: best_bo_matrix(:,:)
    integer, intent(out) :: best_bo_charge
    integer, intent(out) :: best_energy
    integer :: bo_matrix(n_atoms,n_atoms)
    integer :: unsaturated_atoms(n_atoms), degree_unsaturation(n_atoms), unsaturated_bond_matrix(n_atoms,n_atoms)
    integer :: i, j, k, l, m, x,y,n_bonds
    integer, allocatable :: unsaturated_bonds(:,:)
    integer, allocatable :: combinations(:,:,:)
    integer, allocatable :: pairs(:,:)
    integer :: energy, sum_abs_charge, abs_sum_charge, bo_charge, summe
    logical :: healthy

    best_bo_matrix = cn_matrix

    call get_total_charge(n_atoms, atoms, valences, are_core, cn_matrix, energy, sum_abs_charge, abs_sum_charge)

    best_energy = energy

    best_bo_charge = sum_abs_charge

    call get_unsaturated_atoms(n_atoms, valences, cn, unsaturated_atoms, degree_unsaturation)

    call get_unsaturated_bond_matrix(n_atoms, cn_matrix, unsaturated_atoms, unsaturated_bond_matrix, n_bonds)

    if (n_bonds == 0) then
        return
    endif
    
    allocate(unsaturated_bonds(2,n_bonds))

    call get_unsaturated_bonds(unsaturated_bond_matrix, unsaturated_bonds)
    
    k = count(unsaturated_atoms == 1)/2

    if (k == 0) then
        return
    endif 

    allocate(combinations(2, k, n_bonds))

    call get_combinations(n_atoms, unsaturated_bonds, n_bonds, k, combinations)

    allocate(pairs(2,k))

    !write(*,*) "new_valence"
    do i=1, n_bonds

        pairs = combinations(:,:,i)

        !write(*,*) pairs
        
        call get_bo_matrix(n_atoms, cn_matrix, unsaturated_atoms, degree_unsaturation, valences, & 
                            pairs, bo_matrix)
        
        call get_total_charge(n_atoms, atoms, valences, are_core, bo_matrix, energy, sum_abs_charge, abs_sum_charge)

        call check_valence(n_atoms, bo_matrix, valences, healthy)

        !write(*,*) energy, sum_abs_charge, healthy

        if (energy > best_energy .and. sum_abs_charge < best_bo_charge .and. healthy) then

            best_bo_charge = sum_abs_charge
            best_energy = energy
            best_bo_matrix = bo_matrix
            !write(*,*) "best", best_energy, best_bo_charge

        else if (energy == best_energy .and. sum_abs_charge < best_bo_charge .and. healthy) then

            best_bo_charge = sum_abs_charge
            best_energy = energy
            best_bo_matrix = bo_matrix
            !write(*,*) "best", best_energy, best_bo_charge

        else if (energy > best_energy .and. sum_abs_charge == best_bo_charge .and. healthy) then

            best_bo_charge = sum_abs_charge
            best_energy = energy
            best_bo_matrix = bo_matrix
            !write(*,*) "best", best_energy, best_bo_charge

        else if (energy == best_energy .and. sum_abs_charge == best_bo_charge .and. healthy) then

            do j=1, n_atoms
                do k=1, n_atoms
                    if (k > j .and. bo_matrix(j,k) > best_bo_matrix(j,k)) then
                        best_bo_matrix(j,k) = bo_matrix(j,k)
                        best_bo_matrix(k,j) = bo_matrix(k,j)
                    endif
                enddo
            enddo
        endif
    enddo

end subroutine get_bo_matrix_for_valence

subroutine get_bo_matrix(n_atoms, cn_matrix, given_unsaturated_atoms, given_degree_unsaturation, valences, & 
                            starting_pairs, bo_matrix)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: cn_matrix(:,:)
    integer, intent(in) :: given_unsaturated_atoms(:)
    integer, intent(in) :: given_degree_unsaturation(:)
    integer, intent(in) :: valences(:)
    integer, intent(in) :: starting_pairs(:,:)
    integer, intent(out) :: bo_matrix(:,:)
    integer :: i, j, k, l, n_bonds, n_combi
    integer, allocatable :: combinations(:,:,:)
    integer, allocatable :: unsaturated_bonds(:,:)
    integer :: old_degree_unsaturation(n_atoms)
    integer :: current_valences(n_atoms)
    integer :: unsaturated_bond_matrix(n_atoms,n_atoms)
    integer :: current_valence
    integer :: flag
    integer, allocatable :: pairs(:,:), unsaturated_atoms(:), degree_unsaturation(:)

    bo_matrix = cn_matrix

    old_degree_unsaturation = 0

    flag = 0

    allocate(pairs(2,size(starting_pairs, dim=2)))

    pairs = starting_pairs

    allocate(unsaturated_atoms(n_atoms))
    
    unsaturated_atoms = given_unsaturated_atoms

    allocate(degree_unsaturation(n_atoms))

    degree_unsaturation = given_degree_unsaturation

    do while (flag==0)

        flag = 1  

        do j = 1, n_atoms
            if (old_degree_unsaturation(j) /= degree_unsaturation(j)) then
                flag = 0
                exit
            endif
        enddo

        if (flag == 1) then
            exit
        endif

        do i = 1, size(pairs, dim=2)
            if (pairs(1,i) /= 0 .and. pairs(2,i) /= 0) then
                bo_matrix(pairs(1,i), pairs(2,i)) = bo_matrix(pairs(1,i), pairs(2,i)) + 1
                bo_matrix(pairs(2,i), pairs(1,i)) = bo_matrix(pairs(2,i), pairs(1,i)) + 1
            endif 
        enddo

        deallocate(pairs)

        current_valences = 0

        do i=1, size(bo_matrix,dim=2)
            current_valences(i) = sum(bo_matrix(:,i))
        enddo
           
        old_degree_unsaturation = degree_unsaturation

        call get_unsaturated_atoms(n_atoms, valences, current_valences, unsaturated_atoms, degree_unsaturation)
        call get_unsaturated_bond_matrix(n_atoms, cn_matrix, unsaturated_atoms, unsaturated_bond_matrix, n_bonds)

        if (n_bonds == 0) then
            return
        endif

        allocate(unsaturated_bonds(2,n_bonds))
        call get_unsaturated_bonds(unsaturated_bond_matrix, unsaturated_bonds)

        k = count(unsaturated_atoms == 1)/2

        if (k == 0) then
            return
        endif

        allocate(pairs(2,k))

        call get_matching(n_atoms, 1, unsaturated_bonds, k, pairs)

        deallocate(unsaturated_bonds)
    enddo

end subroutine get_bo_matrix

end module bond_order_detection

