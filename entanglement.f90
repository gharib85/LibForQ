!-----------------------------------------------------------------------------------------------------------------------------------
!                                                       Entanglement
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function eentropy(psi, d, da, db, ssys)  ! Returns the ENTANGLEMENT ENTROPY of a bipartite pure state
implicit none
integer :: d, da, db  ! Dimension of whole state space and of the marginal spaces
complex(8) :: psi(1:d)  ! The bipartite state vector
complex(8), allocatable :: rho(:,:)  ! The density operator corresponding to psi, i.e., rho = |psi><psi|
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! Reduced density operators
character(1) :: ssys  ! The sub-system whose entropy is to be computed (ssys is 'a' or 'b')
real(8) :: neumann  ! For the von Neumann's entropy function

allocate( rho(1:d,1:d) ) ;   call projector(psi, d, rho)
if ( ssys == 'a' ) then
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a) ;   deallocate( rho )
  eentropy = neumann(da, rho_a) ;   deallocate( rho_a )
else if ( ssys == 'b' ) then
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b) ;   deallocate( rho )
  eentropy = neumann(db, rho_b) ;   deallocate( rho_b )
endif

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine schmidt_coefficients(psi, d, da, db, schcoeff, eigvec_a, eigvec_b)  ! Returns the Schmidt coeff. for a bipartite pure state
! The matrix evecs has dimension ds x ds and contains in its columns the eigenvectors of rho_s
integer :: d, da, db  ! Dimension of whole state space and of the marginal spaces
complex(8) :: psi(1:d)  ! The bipartite state vector
complex(8), allocatable :: rho(:,:)  ! The density operator corresponding to psi, i.e., rho = |psi><psi|
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! Reduced density operators
real(8), allocatable :: Wa(:), Wb(:)  ! For the eigenvalues of rho_a and rho_b
real(8) :: schcoeff(1:min(da,db))  ! For the Schmidt coefficients
complex(8) :: eigvec_a(1:da,1:da), eigvec_b(1:db,1:db)  ! For the eigenvectors of the reduced density matrices

allocate( rho(1:d,1:d) ) ;   call projector(psi, d, rho)
allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
allocate( rho_b(1:da,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
eigvec_a = rho_a ;   eigvec_b = rho_b ;   deallocate( rho, rho_a, rho_b )
allocate( Wa(1:da), Wb(1:db) ) ;   call lapack_zheevd('V', da, eigvec_a, Wa) ;   call lapack_zheevd('V', db, eigvec_b, Wb)
if ( da <= db ) then !;   allocate( schcoeff(1:da) ) ;
  forall(j=1:da) schcoeff(j) = sqrt(Wa(j))
else if ( da > db ) then !;  allocate( schcoeff(1:db) ) ;
  forall(j=1:db) schcoeff(j) = sqrt(Wb(j))
endif
deallocate( Wa, Wb )

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function concurrence_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(4,4)  ! Density matrix we want to compute the concurrence
complex(8) :: R(4,4), rho_tilde(4,4), s2_kp_s2(4,4)  ! Auxiliary matrices
complex(8) :: egv(4) ! Eigenvalues of R = rho*rho^tilde
real(8) :: egv_max  ! The greater eigenvalue of R
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, s2_kp_s2) ;   rho_tilde = matmul( matmul(s2_kp_s2,conjg(rho)) , s2_kp_s2 )
R = matmul(rho,rho_tilde) ;   call lapack_zgeev('N', 4, R, egv)

egv_max = max( real(egv(1)), real(egv(2)), real(egv(3)), real(egv(4)) )
 concurrence_2qb = max( 0.d0, (2.d0*sqrt(egv_max)-sqrt(real(egv(1)))-sqrt(real(egv(2)))-sqrt(real(egv(3)))-sqrt(real(egv(4)))))

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function EoF_2qb(rho)  ! Returns the entanglement of formation, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(1:4,1:4)  ! Density matrix we want to compute the concurrence
real(8) :: concurrence_2qb  ! For the concurrence function
real(8) :: pv(1:2), shannon  ! Probability vector and Shannon's entropy

pv(1) = (1.d0 + sqrt(1.d0 - concurrence_2qb(rho)**2.d0))/2.d0 ;   pv(2) = 1.d0 - pv(1) ;   EoF_2qb = shannon(2, pv)

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function negativity(d, rho_pt)
  ! Returns the entanglement negativity of a "bipartite" system
! This is an simplified version of the subroutine below.
! Here only the partial transposed matrix is given as input.
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement,
! Phys. Rev. A 65, 032314 (2002).
implicit none
integer :: d ! Dimension of the state space
complex(8) :: rho_pt(1:d,1:d)  ! Partial transposed of a state
real(8) :: norm_tr  ! For the trace norm function

negativity = 0.5d0*(norm_tr(d, rho_pt) - 1.d0)
! It's equal to the sum of the negative eigenvalues of rho_pt

!--------------------------------------------------
!real(8) function negativity(da, db, ssys, rho)  ! Returns the entanglement negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
!implicit none
!character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
!integer :: da, db ! Dimensions of the subsystems
!complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: norm_tr  ! For the trace norm function
!if (ssys == 'a') then
!  call partial_transpose_a(da, db, rho, rho_pt)
!else if (ssys == 'b') then
!  call partial_transpose_b(da, db, rho, rho_pt)
!endif
!negativity = 0.5d0*(norm_tr(da*db, rho_pt) - 1.d0)
!end
!--------------------------------------------------
end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function log_negativity(d, rho_pt)  ! Returns the entanglement logaritmic negativity of a "bipartite" system
! This is an simplified version of the code below. Here only the partial transposed matrix is given as input.
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
implicit none
integer :: d ! Dimension of the state space
complex(8) :: rho_pt(1:d,1:d)  ! Partial transposed of a state
real(8) :: norm_tr  ! For the trace norm function
real(8) :: log2  ! For the log base two
real(8) :: negativity  ! For the negativity of entanglement

!log_negativity = log2(2.d0*negativity(d, rho_pt)+1.d0)
log_negativity = log2( norm_tr(d, rho_pt) )

!--------------------------------------------------
!real(8) function log_negativity(da, db, ssys, rho)  ! Returns the entanglement logaritmic negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
!implicit none
!character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
!integer :: da, db ! Dimensions of the subsystems
!complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: norm_tr  ! For the trace norm function
!real(8) :: log2  ! For the log base two
!if (ssys == 'a') then
!  call partial_transpose_a(da, db, rho, rho_pt)
!else if (ssys == 'b') then
!  call partial_transpose_b(da, db, rho, rho_pt)
!endif
!log_negativity = log2( norm_tr(da*db, rho_pt) )
!end
!--------------------------------------------------
end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine entanglement_hs(d, rho_pt, Ehs, css)  ! Returns the Hilbert-Schmidt entanglement of two qudits
! Ref: J. Maziero, Computing partial transposes and related entanglement measures, Braz. J. Phys. 46, 605 (2016),  arXiv:1609.00323
implicit none
integer :: d, dm, dp, dpp  ! For the dimensions (d is the whole system dimension)
complex(8) :: rho_pt(1:d,1:d)  ! The partial transpose of the state under analysis (input). On exit, if css = 'y' and
                               ! Ehs > 0 then the closest separable state (CSS) is returned via this variable
real(8) :: Ehs  ! For the Hilbert-Schmidt entanglement
character(1) :: css  ! If css = 'y' the CSS is computed and returned in rho_pt, if css = 'n' and/or Ehs = 0 the CSS is not
                     ! computed and rho_pt is not modified
complex(8), allocatable :: A(:,:)  ! Auxiliary variable for sending the PT to Lapack (it may returns the eigenvectors on exit)
real(8), allocatable :: W(:), Wd(:)  ! For the eigenvalues of the PT, in ascending and descending order, respectively
real(8) :: sw, sn1, sn2, sp1, sp2, xi  ! Auxiliary variable for the sums of eigenvalues and for xi
integer :: j  ! Auxiliary variable for counters
complex(8), allocatable :: proj(:,:)  ! Auxiliary variable for projectors

allocate( A(1:d,1:d), W(1:d) ) ;   A = rho_pt  ! Computes the eigenvalues and eigenvectors of the PT
if ( css == 'n' ) then ; call lapack_zheevd('N', d, A, W) ; else if ( css == 'y' ) then ; call lapack_zheevd('V', d, A, W) ; endif

j = 0 ;   dm = 0 ;   do ;   j = j + 1 ;  if ( W(j) >= 0.d0 ) exit ;   dm = j ;   enddo  ! Computes d-

if ( dm == 0 ) then  ! In this case Ehs is null
  Ehs = 0.d0
else if ( dm > 0 ) then  ! Computes the auxiliary dimension and Ehs
  sn1 = 0.d0 ;   sn2 = 0.d0 ;   do j = 1, dm ;   sn1 = sn1 + dabs(W(j)) ;   sn2 = sn2 + (W(j))**2.d0 ;   enddo
  sp1 = 0.d0 ;   sp2 = 0.d0 ;  dp = d - dm ;   allocate( Wd(1:d) ) ;   forall ( j = 1:d ) Wd(j) = W(d-j+1)
  sw = 0.d0 ;   j = 0 ;   do ;   j = j + 1 ;   sw = sw + Wd(j)   ;  if ( (sw > 1.d0) .or. (j >= dp) ) exit ; enddo ;   dpp = j
  if ( dp > dpp ) then ;   do j = dpp+1 , dp ;   sp1 = sp1 + Wd(j) ;   sp2 = sp2 + (Wd(j))**2.d0 ;   enddo ;    endif
  Ehs = dsqrt( (sn1 - sp1)**2.d0 + sp2 + sn2 )
  if ( css == 'y' ) then  ! For computing the closest separable state
    allocate( proj(1:d,1:d) ) ;   rho_pt = 0.d0 ;   xi = 0.d0
    do j = 1, dpp-1 ;   call projector(A(:,d-j+1), d, proj) ;   rho_pt = rho_pt + Wd(j)*proj ;   xi = xi + Wd(j) ;   enddo
    j = dpp ;   xi = 1.d0 - xi ;   call projector(A(:,d-j+1), d, proj) ;   rho_pt = rho_pt + xi*proj ;   deallocate( proj )
  endif ;   deallocate( Wd )
endif

deallocate( A, W )

end
!-----------------------------------------------------------------------------------------------------------------------------------
