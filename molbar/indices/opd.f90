module opd
implicit none
public :: get_index
integer, parameter :: wp = selected_real_kind(15)

contains

subroutine get_index(n_atoms, geometry, priorities, G_os)

    real(8), intent(in) :: geometry(:, :), priorities(:)
    integer, intent(in) ::  n_atoms
    real(8), intent(out) :: G_os
    integer :: i,j,k,l,n,m
    real(8) :: G_ijkl, G, n_atoms_float

    n_atoms_float = n_atoms

    n = 2

    m = 1

    G = 0.0_wp

    do i = 1, n_atoms
        do j = 1, n_atoms
            do k = 1, n_atoms
                do l = 1, n_atoms
                    if (i /= j .and. j /= k .and. k /= l .and. i /= l) then
                        call get_G_ijkl(i,j,k,l, geometry, priorities, G_ijkl)
                        G = G + G_ijkl
                    end if
                end do
            end do
        end do
    end do

    G_os = 24.0_wp/(n_atoms_float**4)*(1.0_wp/3.0_wp)*G

end subroutine

subroutine get_G_ijkl(i,j,k,l, geometry, priorities, G_ijkl)

    real(8), intent(in) :: geometry(:, :), priorities(:)
    integer, intent(in) ::  i, j, k, l
    real(8), intent(out) :: G_ijkl
    real(8) :: r_ij(3), r_kl(3), r_il(3), r_jk(3),c_ijkl(3)
    real(8) :: r_ij_mag, r_kl_mag, r_il_mag, r_jk_mag, mw, d_ijjk, d_jkkl, d_ijklil
    real(8) :: G_ijkl_up, G_ijkl_down

    call get_vector(geometry, i, j, r_ij, r_ij_mag)

    call get_vector(geometry, k, l, r_kl, r_kl_mag)

    call get_vector(geometry, i, l, r_il, r_il_mag)

    call get_vector(geometry, j, k, r_jk, r_jk_mag)

    mw = priorities(i)*priorities(j)*priorities(k)*priorities(l)

    call get_cross(r_ij, r_kl, c_ijkl)

    call get_dot(r_ij, r_jk, d_ijjk)

    call get_dot(r_jk, r_kl, d_jkkl)

    call get_dot(c_ijkl, r_il, d_ijklil)

    G_ijkl_up = d_ijklil*d_ijjk*d_jkkl

    G_ijkl_down = ((r_ij_mag*r_jk_mag*r_kl_mag)**2)*((r_il_mag)**1)

    G_ijkl = mw*G_ijkl_up/G_ijkl_down


end subroutine


subroutine get_vector(geometry, i, j, vec, mag)

    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) ::  i,j
    real(8), intent(out) :: vec(3), mag

    integer :: n

    do n=1, 3

        vec(n) = geometry(n,i)-geometry(n,j)

    enddo

    mag = 0.0_wp

    do n=1, 3

        mag = mag + vec(n)**2

    enddo

    mag = sqrt(mag)

end subroutine

subroutine get_cross(a, b, c)

    real(8), intent(in) :: a(3), b(3)
    real(8), intent(out) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
end subroutine

subroutine get_dot(a, b, c)

    real(8), intent(in) :: a(3), b(3)
    real(8), intent(out) :: c

    c = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)

end subroutine

end module

