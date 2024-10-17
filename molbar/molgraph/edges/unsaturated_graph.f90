module unsaturated_graph
implicit none
public :: get_combinations, get_matching
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine get_combinations(n_atoms, bonds, n_bonds, k, combinations)
    integer, intent(in) :: n_atoms
    integer, intent(in) :: n_bonds
    integer, intent(in) :: k
    integer, intent(in) :: bonds(:,:)
    integer, intent(out) :: combinations(2,k,n_bonds)
    integer :: i
    integer :: matching(2,k)

    do i=1, n_bonds
        call get_matching(n_atoms, i, bonds, k, matching)
        !write(*,*) "Final matching", matching
        combinations(:,:,i) = matching
    enddo

end subroutine get_combinations

subroutine get_matching(n_atoms, starting_id, bonds, k, matching)
    implicit none
    integer, intent(in), dimension(:, :) :: bonds
    integer, intent(in) :: n_atoms, starting_id, k
    integer, intent(out) :: matching(:,:)

    integer :: edge_id, i, j, n_degenerate_atoms, lowest_degree
    integer, dimension(:), allocatable :: degrees
    integer :: next_edge(2)
    integer, dimension(:), allocatable :: lowest_degree_atoms
    integer, dimension(:,:), allocatable :: next_edges
    logical :: success
    integer, dimension(:), allocatable :: visited

    allocate(visited(n_atoms))

    allocate(degrees(n_atoms))

    matching = 0

    matching(:, 1) = bonds(:, starting_id)

    edge_id = 2

    visited = 0
    visited(matching(1, 1)) = 1
    visited(matching(2, 1)) = 1

    do i = 1, k

        !write(*,*) i, "matching", matching
        !write(*,*) i, "visited", visited

        call get_degrees(bonds, visited, degrees)

        !write(*,*) i, "th degrees", degrees

        call get_lowest_degree(degrees, visited, n_degenerate_atoms, lowest_degree, success)

        if (.not. success) then
            return
        end if

        allocate(lowest_degree_atoms(n_degenerate_atoms))

        call get_lowest_degree_atoms(degrees, visited, lowest_degree, lowest_degree_atoms)
        
        !write(*,*) i, "th lowest degree atoms", lowest_degree_atoms

        !allocate(next_edges(2, n_degenerate_atoms))

        do j=1, size(lowest_degree_atoms)

            call get_next_edge_for_atom(lowest_degree_atoms(j), bonds, visited, degrees, next_edge, success)

            if (.not. success) then

                cycle

            else
                matching(:, edge_id) = next_edge
                visited(next_edge(1)) = 1
                visited(next_edge(2)) = 1
                edge_id = edge_id + 1
                EXIT
            end if
        end do

        deallocate(lowest_degree_atoms)

    end do

    deallocate(visited)
end subroutine get_matching

!subroutine get_matching(n_atoms, starting_edge_id, bonds, k, matching)

!> Calculates the degree of a given atom in a molecule.
!! Inputs:
!!   atom: integer, the index of the atom whose degree is to be calculated.
!!   bonds: integer, 2D array containing the indices of the atoms that are bonded.
!!   visited: integer, 1D array containing the indices of the atoms that have already been visited.
!! Outputs:
!!   degree: integer, the degree of the given atom.
subroutine get_degree(atom, bonds, visited, degree)
    implicit none
    integer, intent(in) :: atom
    integer, intent(in), dimension(:,:) :: bonds
    integer, intent(in), dimension(:) :: visited
    integer, intent(out) :: degree

    integer :: i, trial_atom

    if (visited(atom) == 1) then
        degree = 0
        return
    end if

    degree = 0
    do i = 1, size(bonds,2)
        if (atom == bonds(1,i) .or. atom == bonds(2,i)) then
            if (atom == bonds(1,i)) then
                trial_atom = bonds(2,i)
            else
                trial_atom = bonds(1,i)
            end if

            if (visited(trial_atom) == 0) then
                degree = degree + 1
            end if
        end if
    end do

end subroutine get_degree

!> Calculates the degree of each vertex in a graph represented by a bonds.
!! Arguments:
!! bonds: integer, intent(in), dimension(:,:) :: bonds
!! visited: integer, intent(in), dimension(:) :: array of visited vertices
!! degrees: integer, intent(out), dimension(:) :: array of degrees of each vertex
subroutine get_degrees(bonds, visited, degrees)
    implicit none
    integer, intent(in), dimension(:,:) :: bonds
    integer, intent(in), dimension(:) :: visited
    integer, intent(out), dimension(:) :: degrees

    integer :: i
    degrees = 0
    do i = 1, size(visited)
        call get_degree(i, bonds, visited, degrees(i))
    end do

end subroutine get_degrees

subroutine get_lowest_degree(degrees, visited, n_degenerate_atoms, lowest_degree, success)

    implicit none
    integer, intent(in), dimension(:) :: degrees
    integer, intent(in), dimension(:) :: visited
    integer, intent(out) :: n_degenerate_atoms, lowest_degree
    logical, intent(out) :: success

    integer :: i

    lowest_degree = HUGE(1)
    n_degenerate_atoms = 0
    success = .true.

    do i = 1, size(degrees)
        if (degrees(i) < lowest_degree .and. degrees(i) /= 0 .and. visited(i) == 0) then
            lowest_degree = degrees(i)
            n_degenerate_atoms = 1
        elseif (degrees(i) == lowest_degree .and. visited(i) == 0) then
            n_degenerate_atoms = n_degenerate_atoms + 1
        end if
    end do
    if (lowest_degree > 1000) then
        !write(*,*) "lowest_degree", lowest_degree, "n_degenerate_atoms", n_degenerate_atoms
        success = .false.
        return
    end if

end subroutine get_lowest_degree

subroutine get_lowest_degree_atoms(degrees, visited, lowest_degree, lowest_degree_atoms)
    implicit none
    integer, intent(in), dimension(:) :: degrees
    integer, intent(in), dimension(:) :: visited
    integer, intent(in) :: lowest_degree
    integer, intent(out), dimension(:) :: lowest_degree_atoms
    integer :: i, idx

    idx = 1

    do i = 1, size(degrees)
        if (degrees(i) == lowest_degree .and. visited(i) == 0) then
            lowest_degree_atoms(idx) = i
            idx = idx + 1
        end if
    end do

end subroutine get_lowest_degree_atoms

subroutine get_next_edge_for_atom(atom, bonds, visited, degrees, next_edge, success)
    implicit none
    integer, intent(in) :: atom
    integer, intent(in), dimension(:, :) :: bonds
    integer, intent(in), dimension(:) :: visited
    integer, intent(in), dimension(:) :: degrees
    integer, intent(out), dimension(2) :: next_edge
    logical, intent(out) :: success

    integer :: i, n_bonds, j
    integer, dimension(:), allocatable :: trial_bonds, trial_degrees
    integer :: min_degree, next_atom

    n_bonds = 0

    next_atom = 0

    next_edge = 0
    ! Count the number of bonds for the given atom
    do i = 1, size(bonds, 2)
        if (bonds(1,i) == atom .and. visited(bonds(2, i)) == 0) then
            n_bonds = n_bonds + 1
        elseif (bonds(2,i) == atom .and. visited(bonds(1,i)) == 0) then
            n_bonds = n_bonds + 1
        end if
    end do

    if (n_bonds == 0) then
        success = .false.
        return
    else
        allocate(trial_bonds(n_bonds))
        allocate(trial_degrees(n_bonds))
    end if

    j = 1

    ! Populate trial_bonds array
    do i = 1, size(bonds, 2)
        if (bonds(1,i) == atom .and. visited(bonds(2,i)) == 0) then
            trial_bonds(j) = bonds(2,i)
            j = j + 1
        elseif (bonds(2,i) == atom .and. visited(bonds(1,i)) == 0) then
            trial_bonds(j) = bonds(1,i)
            j = j + 1
        end if
    end do

    ! Populate trial_degrees array
    do i = 1, n_bonds
        trial_degrees(i) = degrees(trial_bonds(i))
    end do

    min_degree = HUGE(1)

    ! Find the minimum degree and corresponding atom
    do i = 1, n_bonds
        if (trial_degrees(i) < min_degree) then
            min_degree = trial_degrees(i)
            next_atom = trial_bonds(i)
        end if
    end do

    ! Set the output variables
    next_edge(1) = atom
    next_edge(2) = next_atom
    success = .true.
    ! Deallocate the temporary arrays
    deallocate(trial_bonds)
    deallocate(trial_degrees)

end subroutine get_next_edge_for_atom

subroutine get_next_edges_for_atoms(atoms, bonds, visited, degrees, next_edges, final_success)
    implicit none
    integer, intent(in), dimension(:) :: atoms
    integer, intent(in), dimension(:, :) :: bonds
    integer, intent(inout), dimension(:) :: visited
    integer, intent(in), dimension(:) :: degrees
    integer, intent(out), dimension(:,:) :: next_edges
    logical, intent(out) :: final_success
    integer, dimension(:), allocatable :: updated_degrees

    integer :: i
    integer:: next_edge(2)
    logical :: success

    next_edges = 0

    final_success = .false.

    allocate(updated_degrees(size(degrees)))

    updated_degrees = degrees

    do i = 1, size(atoms)
        call get_next_edge_for_atom(atoms(i), bonds, visited, updated_degrees, next_edge, success)
        !write(*,*) "next_edge", next_edge
        if (success .and. visited(next_edge(1)) == 0 .and. visited(next_edge(2)) == 0) then
            final_success = .true.
            next_edges(:,i) = next_edge
            visited(next_edge(1)) = 1
            visited(next_edge(2)) = 1
            call get_degrees(bonds, visited, updated_degrees)
        end if
    end do

end subroutine get_next_edges_for_atoms
end module unsaturated_graph
