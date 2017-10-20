!-----------------------------------------------------------------------------------------------------------------------------------
!                                                 Some popular Quantum STATES
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine psi_qubit(theta, phi, psi)  ! General one-qubit pure state
implicit none
real(8) :: theta, phi
complex(8) :: psi(4) 

psi(1) = cos(theta/2.d0) ;  psi(2) = sin(theta/2.d0)*(cos(phi) + sin(phi)*(0.d0,1.d0))

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_qubit(r1, r2, r3, rho)  ! General one-qubit mixed state
implicit none
real(8) :: r1, r2, r3  ! Components of Bloch's vector
complex(8) :: rho(2,2), s0(2,2), s1(2,2), s2(2,2), s3(2,2)

call pauli_group(s0, s1, s2, s3);  rho = 0.5d0*(s0 + r1*s1 + r2*s2 + r3*s3)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine bell_basis(psi_p, psi_m, phi_p, phi_m)  ! Defines the BELL BASIS states
implicit none
complex(8) :: psi_p(4), psi_m(4), phi_p(4), phi_m(4) 

psi_p(1) = 0.d0 ;  psi_p(2) = 1.d0/sqrt(2.d0) ;   psi_p(3) = 1.d0/sqrt(2.d0) ;   psi_p(4) = 0.d0
psi_m(1) = 0.d0 ;  psi_m(2) = 1.d0/sqrt(2.d0) ;   psi_m(3) = -1.d0/sqrt(2.d0) ;   psi_m(4) = 0.d0
phi_p(1) = 1.d0/sqrt(2.d0) ;  phi_p(2) = 0.d0 ;   phi_p(3) = 0.d0 ;   phi_p(4) = 1.d0/sqrt(2.d0)
phi_m(1) = 1.d0/sqrt(2.d0) ;  phi_m(2) = 0.d0 ;   phi_m(3) = 0.d0 ;   phi_m(4) = -1.d0/sqrt(2.d0)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_bds(c1, c2, c3, rho)  ! Two-qubit Bell-diagonal state (BDS)
implicit none
real(8) :: c1, c2, c3  ! Components of the correlation vector
complex(8) :: rho(1:4,1:4)  ! For the Bell-diagobal density matrix
complex(8), allocatable :: kp(:,:)  ! For the Kronecker product between elements of the Pauli group
complex(8), allocatable :: sigma_0(:,:), sigma_1(:,:), sigma_2(:,:), sigma_3(:,:)

allocate( kp(1:4,1:4), sigma_0(1:2,1:2), sigma_1(1:2,1:2), sigma_2(1:2,1:2), sigma_3(1:2,1:2) )
call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)

rho = 0.d0
call kronecker_product_c(sigma_0, 2, 2, sigma_0, 2, 2, kp) ;   rho = rho + kp  ! Identity term
call kronecker_product_c(sigma_1, 2, 2, sigma_1, 2, 2, kp) ;   rho = rho + c1*kp  ! sigma_x term
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, kp) ;   rho = rho + c2*kp  ! sigma_y term
call kronecker_product_c(sigma_3, 2, 2, sigma_3, 2, 2, kp) ;   rho = rho + c3*kp  ! sigma_z term
rho = (1.d0/4.d0)*rho

deallocate( kp, sigma_0, sigma_1, sigma_2, sigma_3 )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_x(c11, c22, c33, a3, b3, rho)  ! Two-qubit X state, in the standar form
implicit none
real(8) :: c11, c22, c33, a3, b3  ! Components of the correlation vector and polarizations
complex(8) :: rho(1:4,1:4)  ! For the Bell-diagobal density matrix
complex(8), allocatable :: kp(:,:)  ! For the Kronecker product between elements of the Pauli group
complex(8), allocatable :: sigma_0(:,:), sigma_1(:,:), sigma_2(:,:), sigma_3(:,:)

allocate( kp(1:4,1:4), sigma_0(1:2,1:2), sigma_1(1:2,1:2), sigma_2(1:2,1:2), sigma_3(1:2,1:2) )
call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)

rho = 0.d0
call kronecker_product_c(sigma_0, 2, 2, sigma_0, 2, 2, kp) ;   rho = rho + kp  ! Identity term
call kronecker_product_c(sigma_1, 2, 2, sigma_1, 2, 2, kp) ;   rho = rho + c11*kp  ! sigma_x term
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, kp) ;   rho = rho + c22*kp  ! sigma_y term
call kronecker_product_c(sigma_3, 2, 2, sigma_3, 2, 2, kp) ;   rho = rho + c33*kp  ! sigma_z term
call kronecker_product_c(sigma_3, 2, 2, sigma_0, 2, 2, kp) ;   rho = rho + a3*kp  ! system A z magnetization
call kronecker_product_c(sigma_0, 2, 2, sigma_3, 2, 2, kp) ;   rho = rho + b3*kp  ! system B z magnetization
rho = (1.d0/4.d0)*rho

deallocate( kp, sigma_0, sigma_1, sigma_2, sigma_3 )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_werner(da, x, rho)  ! Returns the two-qudit WERNER's state, as written in the reference below:
! Ref: S. J. Akhtarshenas, H. Mohammadi, S. Karimi, and Z. Azmi, Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
integer :: da ! Dimension of sub-system a (db = da and d = da*db = da^2)
complex(8) :: rho(1:da**2,1:da**2)  ! The density matrix (output)
real(8) :: x   ! For the "mixedness" parameter x (x \in [-1,1])
real(8) :: y  ! Auxiliary variable for computing the matrix elements
integer :: k, l  ! Auxiliary variables for counters

rho = 0.d0
! 'Identity' term
y = (dble(da)-x)/dble(da*(da**2-1))
do k = 1, da**2 ;   rho(k,k) = y ;   enddo
! Off-diagonal terms
y = (dble(da)*x - 1.d0)/dble(da*(da**2-1))
do k = 1, da ;   do l = 1, da
  if ( k == l ) then ;   rho((k-1)*da+l,(l-1)*da+k) = rho((k-1)*da+l,(l-1)*da+k) + y
  else if ( k /= l ) then ;   rho((k-1)*da+l,(l-1)*da+k) = y
  endif
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_isotropic(da, x, rho)  ! Returns the two-qudit ISOTROPIC state
! Ref: S. J. Akhtarshenas, H. Mohammadi, S. Karimi, and Z. Azmi, Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
integer :: da ! Dimension of sub-system a (db = da and d = da*db = da^2)
complex(8) :: rho(1:da**2,1:da**2)  ! The density matrix (output)
real(8) :: x   ! For the "mixedness" parameter x (x \in [-1,1])
real(8) :: y  ! Auxiliary variable for computing the matrix elements
integer :: k, l  ! Auxiliary variables for counters

rho = 0.d0
! 'Identity' term
y = (1.d0 - x)/dble(da**2-1)
do k = 1, da**2 ;   rho(k,k) = y ;   enddo
! Off-diagonal terms
y = (x*dble(da**2) - 1.d0)/dble(da**2-1)
do k = 1, da ;   do l = 1, da
  if ( k == l ) then ;   rho((k-1)*da+k,(l-1)*da+l) = rho((k-1)*da+k,(l-1)*da+l) + y/dble(da)
  else if ( k /= l ) then ;   rho((k-1)*da+k,(l-1)*da+l) = y/dble(da)
  endif
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_GHZ(nqb, GHZ)  ! Returns the ket for the GHZ state
implicit none
integer :: nqb  ! No. of qubits
integer :: d  ! Dimension
complex(8) :: GHZ(1:2**nqb)  ! Ket for the GHZ state (of a number nqb of qubits)
real(8) :: x  ! Auxiliary variable

x = 1.d0/dsqrt(2.d0) ;   d = 2**nqb ;   GHZ = 0.d0 ;   GHZ(1) = x ;   GHZ(d) = x

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_W(nqb, W)  ! Returns the ket for the W state
implicit none
integer :: nqb  ! No. of qubits
integer :: d  ! Dimension of the ket
complex(8) :: W(1:2**nqb)  ! Ket for the W state (of a number nqb of qubits)
real(8) :: x  ! Auxiliary variable
integer :: j  ! Auxiliary variable for counters
 
d = 2**nqb
if ( nqb == 1 ) then ;  x = 1.d0/sqrt(dble(d)) ;   W = x ;   return ;   endif  ! For one qubit
! For two or more qubits
x = 1.d0/sqrt(dble(nqb)) ;   W = 0.d0 ;   W(2) = x ;   W(2**(nqb-1)+1) = x 
if ( nqb > 2 ) then ;   do j = nqb-1, 2, -1 ;   W(2**(nqb-j)+1) = x ;   enddo ;   endif

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine state_thermal(kBT, H, d, rhoT)  ! Returns the Gibbs THERMAL state 
implicit none
integer :: d  ! Dimension
complex(8) :: H(1:d,1:d), rhoT(1:d,1:d)  ! For the Hamiltonian and the associated thermal state
complex(8), allocatable :: A(:,:), proj(:,:) ! Auxiliary matrices
real(8), allocatable :: W(:)  ! For the Hamiltonian eigenvalues
real(8) :: kBT  ! For the temperature, multiplied by the Boltzmann constant
real(8) :: ZZ  ! For the partition function
integer :: j, k  ! Auxiliary variables for counters

allocate(A(1:d,1:d), proj(1:d,1:d), W(1:d))
if ( kBT < 1.d-15) kBT = 1.d-15 ;   A = H ;   call lapack_zheevd('V', d, A, W)                            
forall(j = 1:d) W(j) = W(j) - W(1) ! Changes the energy scale: E >= 0
ZZ = 0.d0 ;   do j = 1, d ;   ZZ = ZZ + dexp(-W(j)/kBT) ;   enddo ;   if ( ZZ < 1.d-15) kBT = 1.d-15 
rhoT = 0.d0 ;   do j = 1, d ;   call projector(A(:,j), d, proj) ;  rhoT = rhoT + (exp(-W(j)/kBT)/ZZ)*proj ;   enddo
deallocate(A, proj, W)

end 
!-----------------------------------------------------------------------------------------------------------------------------------