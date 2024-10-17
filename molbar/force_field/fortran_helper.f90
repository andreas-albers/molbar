module fortran_helper

implicit none
public :: gidx, get_R_ij, get_theta_ijk, get_cosphi_ijkl, get_sinphi_ijkl, get_aijkl, get_bijk, get_cjkl, get_dijkl
integer, parameter :: wp = selected_real_kind(15)
contains

subroutine gidx(i, j, ij)
    integer, intent(in) :: i, j
    integer, intent(out) :: ij
    if (i < j) then
        ij = (j - 1) * j  / 2 + i 
    else
        ij = (i - 1)* i / 2 + j 
    endif
end subroutine gidx

subroutine get_R_ij(geometry, i, j, R_ij)

    intrinsic :: sqrt   

    real(8), intent(in) :: geometry(:, :)

    integer, intent(in) :: i, j

    real(8), intent(out) :: R_ij

    real(8) :: x_ij, y_ij, z_ij

    x_ij = geometry(1,i)-geometry(1,j)
    y_ij = geometry(2,i)-geometry(2,j)
    z_ij = geometry(3,i)-geometry(3,j)

    R_ij = sqrt(x_ij**2 + y_ij**2 + z_ij**2)

end subroutine get_R_ij

subroutine get_theta_ijk(geometry, i, j, k, theta)

    intrinsic :: dot_product, acos

    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k
    real(8), intent(out) :: theta
    real(8) :: r_ij(3), r_kj(3)
    real(8) :: numerator,  denominator, norm_r_ij, norm_r_kj

    r_ij = (/ geometry(1,i)-geometry(1,j), geometry(2,i)-geometry(2,j), geometry(3,i)-geometry(3,j) /)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j) /)

    call get_l2norm(r_ij, norm_r_ij)

    call get_l2norm(r_kj, norm_r_kj)

    numerator = dot_product(r_ij, r_kj)

    denominator = norm_r_ij*norm_r_kj

    theta = acos(numerator/denominator)

end subroutine get_theta_ijk

subroutine get_aijk(geometry, i, j, k, aijk)

    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k
    real(8), intent(out) :: aijk
    real(8) :: r_ij(3), r_kj(3)

    r_ij = (/ geometry(1,i)-geometry(1,j), geometry(2,i)-geometry(2,j), geometry(3,i)-geometry(3,j)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)

    aijk = dot_product(r_ij, r_kj)

end subroutine get_aijk

subroutine get_bij(geometry, i, j, bij)

    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j
    real(8), intent(out) :: bij
    real(8) :: r_ij(3)

    r_ij = (/ geometry(1,i)-geometry(1,j), geometry(2,i)-geometry(2,j), geometry(3,i)-geometry(3,j)/)

    call get_l2norm(r_ij, bij)

end subroutine get_bij

subroutine get_aijkl(geometry, i, j, k, l, aijkl)
    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k, l
    real(8), intent(out) :: aijkl
    real(8) :: r_ji(3), r_kj(3), r_lk(3), cross_jikj(3), cross_kjlk(3)

    r_ji = (/ geometry(1,j)-geometry(1,i), geometry(2,j)-geometry(2,i), geometry(3,j)-geometry(3,i)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)
    
    r_lk = (/ geometry(1,l)-geometry(1,k), geometry(2,l)-geometry(2,k), geometry(3,l)-geometry(3,k)/)

    call cross_product(r_ji, r_kj, cross_jikj)

    call cross_product(r_kj, r_lk, cross_kjlk)

    aijkl = dot_product(cross_jikj, cross_kjlk)

end subroutine get_aijkl

subroutine get_bijk(geometry, i, j, k, bijk)
    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k
    real(8), intent(out) :: bijk
    real(8) :: r_ji(3), r_kj(3), cross_jikj(3)

    r_ji = (/ geometry(1,j)-geometry(1,i), geometry(2,j)-geometry(2,i), geometry(3,j)-geometry(3,i)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)
    
    call cross_product(r_ji, r_kj, cross_jikj)

    call get_l2norm(cross_jikj, bijk)

end subroutine get_bijk

subroutine get_cjkl(geometry, j, k,l, cjkl)
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: j,k,l
    real(8), intent(out) :: cjkl
    real(8) ::r_kj(3), r_lk(3), cross_kjlk(3)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)

    r_lk = (/ geometry(1,l)-geometry(1,k), geometry(2,l)-geometry(2,k), geometry(3,l)-geometry(3,k)/)
    
    call cross_product(r_kj, r_lk, cross_kjlk)

    call get_l2norm(cross_kjlk, cjkl)

end subroutine get_cjkl

subroutine get_dijkl(geometry, i, j, k, l, dijkl)
    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k, l
    real(8), intent(out) :: dijkl
    real(8) :: r_ji(3), r_kj(3), r_lk(3), cross_jikj(3), cross_kjlk(3)
    real(8) :: norm_r_kj,dot_jikjlk

    r_ji = (/ geometry(1,j)-geometry(1,i), geometry(2,j)-geometry(2,i), geometry(3,j)-geometry(3,i)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)
    
    r_lk = (/ geometry(1,l)-geometry(1,k), geometry(2,l)-geometry(2,k), geometry(3,l)-geometry(3,k)/)

    call cross_product(r_ji, r_kj, cross_jikj)

    call cross_product(r_kj, r_lk, cross_kjlk)

    call get_l2norm(r_kj, norm_r_kj)

    dot_jikjlk = dot_product(r_ji, cross_kjlk)

    dijkl = norm_r_kj*dot_jikjlk

end subroutine get_dijkl



subroutine get_cosphi_ijkl(geometry, i, j, k, l, cosphi)
    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k, l
    real(8), intent(out) :: cosphi
    real(8) :: r_ji(3), r_kj(3), r_lk(3), cross_jikj(3), cross_kjlk(3)
    real(8) :: numerator, denominator, norm_cross_jikj, norm_cross_kjlk

    r_ji = (/ geometry(1,j)-geometry(1,i), geometry(2,j)-geometry(2,i), geometry(3,j)-geometry(3,i)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)
    
    r_lk = (/ geometry(1,l)-geometry(1,k), geometry(2,l)-geometry(2,k), geometry(3,l)-geometry(3,k)/)

    call cross_product(r_ji, r_kj, cross_jikj)

    call cross_product(r_kj, r_lk, cross_kjlk)

    numerator = dot_product(cross_jikj, cross_kjlk)

    call get_l2norm(cross_jikj, norm_cross_jikj)

    call get_l2norm(cross_kjlk, norm_cross_kjlk)

    denominator = norm_cross_jikj * norm_cross_kjlk

    cosphi = numerator/denominator

end subroutine get_cosphi_ijkl

subroutine get_sinphi_ijkl(geometry, i, j, k, l, sinphi)
    intrinsic :: dot_product
    real(8), intent(in) :: geometry(:, :)
    integer, intent(in) :: i, j, k, l
    real(8), intent(out) :: sinphi
    real(8) :: r_ji(3), r_kj(3), r_lk(3), cross_jikj(3), cross_kjlk(3)
    real(8) :: numerator, denominator, norm_cross_jikj, norm_cross_kjlk, norm_r_kj,dot_jikjlk

    r_ji = (/ geometry(1,j)-geometry(1,i), geometry(2,j)-geometry(2,i), geometry(3,j)-geometry(3,i)/)

    r_kj = (/ geometry(1,k)-geometry(1,j), geometry(2,k)-geometry(2,j), geometry(3,k)-geometry(3,j)/)
    
    r_lk = (/ geometry(1,l)-geometry(1,k), geometry(2,l)-geometry(2,k), geometry(3,l)-geometry(3,k)/)

    call cross_product(r_ji, r_kj, cross_jikj)

    call cross_product(r_kj, r_lk, cross_kjlk)

    call get_l2norm(r_kj, norm_r_kj)

    dot_jikjlk = dot_product(r_ji, cross_kjlk)

    numerator = norm_r_kj*dot_jikjlk

    call get_l2norm(cross_jikj, norm_cross_jikj)

    call get_l2norm(cross_kjlk, norm_cross_kjlk)

    denominator = norm_cross_jikj * norm_cross_kjlk

    sinphi = numerator/denominator

end subroutine get_sinphi_ijkl

subroutine cross_product(vec1, vec2, vec3)

    real(8), intent(in) :: vec1(3), vec2(3)
    real(8), intent(out) :: vec3(3)

    vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
    
end subroutine cross_product

subroutine get_l2norm(vec1, norm)

    intrinsic :: sqrt, selected_real_kind
    integer, parameter :: wp = selected_real_kind(15)
    real(8), intent(in) :: vec1(:)
    real(8), intent(out) :: norm
    integer :: i

    norm = 0.0_wp

    do i=1, size(vec1)
        norm = norm + vec1(i)**2
    enddo

    norm = sqrt(norm)

end subroutine

end module fortran_helper
