!-----------------------------------------------------------------------------------------------------------------------------------
!                                                        Coherence
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function coh_l1n(d, rho)  ! Returns the l1-norm quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014).
! C_l1 = sum_j/=k |rho(j,k)|
implicit none
integer :: d ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
integer :: j, k  ! Auxiliary variables for counters

 coh_l1n = 0.d0
 do j = 1, d-1 ;   do k = j+1, d
   coh_l1n = coh_l1n + sqrt( (dble(rho(j,k)))**2.d0 + (aimag(rho(j,k)))**2.d0 )
 enddo ;   enddo
 coh_l1n = 2.d0*coh_l1n

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function coh_re(d, rho)  ! Returns the relative entropy of quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014).
! C_re = H(rho(j,j)) - S(rho)
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
real(8) :: shannon, neumann  ! For the Shannon and von Neumman entropy functions
real(8) :: pv(1:d)  ! Auxiliary probability vector
integer :: j  ! Auxiliary variable for counters

forall(j=1:d) pv(j) = dble(rho(j,j)) ;   coh_re = shannon(d, pv) - neumann(d, rho)

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function asy_wyd(d,s,rho,L)
! Retorns the Wigner-Yanase-Dyson asymmetry of rho in relation to the generator L: tr(ρL^2 − tr(ρ^s Lρ^1−s L)
! I. Marvian and R. W. Spekkens, Nature Communications 5, 3821 (2014)
integer :: d ! system Dimension
real(8) :: s
complex(8) :: rho(d,d), L(d,d), As(d,d), A1ms(d,d)

asy_wyd = trace(d, matmul(rho,matmul(L,L)))
call mat_pows(d, rho, s, As);  call mat_pows(d, rho, 1.d0-s, A1ms)
asy_wyd = asy_wyd - trace(d, matmul(As,matmul(L,matmul(A1ms,L))))

end
!-----------------------------------------------------------------------------------------------------------------------------------
