module valence_conditions
implicit none
public :: get_total_charge, get_atomic_charge, check_bo_matrix,check_valence
integer, parameter :: wp = selected_real_kind(15)

integer, parameter :: valence_electrons(118) = (/ 1, 2, &
                    1, 2, 3, 4, 5, 6, 7, 8, &
                    1, 2, 3, 4, 5, 6, 7, 8, &
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, &
                    3, 4, 5, 6, 7, 8, &
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, &
                    3, 4, 5, 6, 7, 8, &
                    1, 2, &
                    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, &
                    4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8, & 
                    1, 2, &
                    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, &
                    4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 5, 6, 7, 8 /)
contains


subroutine get_total_charge(n_atoms, atoms, ideal_valence, are_core, bo_matrix, energy, sum_abs_charge, abs_sum_charge)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: bo_matrix(:,:)
    integer, intent(in) :: ideal_valence(:)
    integer, intent(in) :: are_core(:)
    integer, intent(in) :: atoms(:)
    integer, intent(out) :: energy, sum_abs_charge, abs_sum_charge
    integer :: current_valences(n_atoms)
    integer :: current_valence,sum_charge
    integer :: i, j, k, l, q
    !write(*,*) "are_core", are_core
    current_valences = 0
    energy = 0 
    do k = 1, n_atoms
        if (are_core(k) == 0) cycle
        current_valence = 0
        do l = 1, n_atoms
            current_valence = current_valence + bo_matrix(k,l)
            energy = energy + bo_matrix(k,l)
        enddo
        current_valences(k) = current_valence
    enddo

    sum_abs_charge = 0
    sum_charge = 0

    do i=1, n_atoms
        if (are_core(i) == 0) cycle
        call get_atomic_charge(atoms(i), current_valences(i), q)
        !write(*,*) i, atoms(i), current_valences(i), q
        !q = ideal_valence(i) - current_valences(i)
        sum_abs_charge = sum_abs_charge + abs(q)
        sum_charge = sum_charge + q
    enddo
    
    abs_sum_charge = abs(sum_charge)

end subroutine get_total_charge

subroutine get_atomic_charge(atom, current_valence, charge)
    integer, intent(in) :: atom
    integer, intent(in) :: current_valence
    integer, intent(out) :: charge
    integer :: i

    if (atom == 1) then
        charge = 1 - current_valence
    else if (atom == 5) then
        charge = 3 - current_valence
    else if (atom == 7 .and. current_valence==3) then
        charge = 0
    else if (atom == 7 .and. current_valence==4) then
        charge = 1
    else if (atom == 7 .and. current_valence==2) then
        charge = -1
    else if (atom == 7 .and. current_valence==1) then
        charge = -2
    else if (atom == 8 .and. current_valence==2) then
        charge = 0
    else if (atom == 8 .and. current_valence==1) then
        charge = -1
    else if (atom == 15 .and. current_valence==5) then
        charge = 0
    else if (atom == 16 .and. current_valence==6) then
        charge = 0
    else if (current_valence /= 0) then
        charge = valence_electrons(atom) - current_valence
    else
        charge = 0
    endif

end subroutine get_atomic_charge

subroutine check_bo_matrix(n_atoms, cn_matrix, bo_matrix, valence, degree_unsaturation, healthy)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: cn_matrix(:,:), bo_matrix(:,:)
    integer, intent(in) :: valence(:)
    integer, intent(in) :: degree_unsaturation(:)
    logical, intent(out) :: healthy
    logical :: checked_valence
    integer :: degree, i, j, sum_of_additional_bo

    call check_valence(n_atoms, bo_matrix, valence, checked_valence)

    if (checked_valence .eqv. .false.) then
        healthy = .false.
        return
    endif

    degree  = sum(degree_unsaturation)
    sum_of_additional_bo = 0
    do i=1, size(cn_matrix, dim=1)
        do j=1, size(cn_matrix, dim=2)
            sum_of_additional_bo = sum_of_additional_bo + bo_matrix(i,j) - cn_matrix(i,j)
        enddo
    enddo

    if (sum_of_additional_bo /= degree) then
        healthy = .false.
        return
    endif

    healthy = .true.

end subroutine check_bo_matrix

subroutine check_valence(n_atoms, bo_matrix, valence, healthy)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: bo_matrix(:,:)
    integer, intent(in) :: valence(:)
    logical, intent(out) :: healthy
    integer :: i
    integer :: actual_valence_per_atom(n_atoms)

    do i=1, n_atoms
        actual_valence_per_atom(i) = sum(bo_matrix(:,i))
    enddo
    do i=1, n_atoms
        if (actual_valence_per_atom(i) > valence(i)) then
            healthy = .false.
            return
        endif
    enddo

    healthy = .true.

end subroutine check_valence

end module valence_conditions