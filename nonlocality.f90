!-----------------------------------------------------------------------------------------------------------------------------------
!                                            Nonlocality (quantifiers and related functions)
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function CHSH_2qb(rho)  ! Returns the CHSH parameter (its maximazation)
! Ref: R. Horodecki, P. Horodecki e M. Horodecki, “Violating Bell inequality by mixed spin-1 states: Necessary and sufficient
!      condition”, Phys. Lett. A 200, 3402 (1995).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: mK(1:3,1:3)  ! Auxiliary matrix
real(8) :: W(1:3)  ! For the eigenvalues of K

call corrmat_gellmann_unopt(2, 2, rho, corr)
mK = matmul(transpose(corr),corr) ;
call lapack_dsyevd('N', 3, mK, W)
 CHSH_2qb = 2.d0*sqrt( W(2) + W(3) )

end
!-----------------------------------------------------------------------------------------------------------------------------------
