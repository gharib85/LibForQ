!###################################################################################################################################
!                                            Nonlocality (quantifiers and related functions)
!###################################################################################################################################
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

call stokes_parameters_2qb(rho, ma, mb, corr) ;   mK = matmul(transpose(corr),corr) ;   call lapack_dsyevd('N', 3, mK, W)
 CHSH_2qb = 2.d0*sqrt( W(2) + W(3) ) 

end
!-----------------------------------------------------------------------------------------------------------------------------------
! needs more tests
real(8) function min_hs(ssys, da, db, rho)  ! Returns the MEASUREMENT-INDUCED NONLOCALITY (with Hilbert-Schmidt distance)
! Refs: Phys. Rev. Lett. 106, 120401 (2011)
implicit none
character(1), intent(in) :: ssys  ! Tells if sub-system a or b is measured, in the minimization
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! The reduced states of the subsystems
real(8), allocatable :: bv(:)  ! For the Bloch vector, of rho_a or of rho_b
real(8), allocatable :: corrmat(:,:)  ! For the correlation matrix, as defined in arXiv:1603.05284
real(8), allocatable :: Xi(:,:)  ! For the product of the correlation matrix and its transpose, as defined in PRL 106, 120401 (2011)
real(8), allocatable :: W(:)  ! For the eigenvalues of Xi
integer :: dda, ddb  ! Auxiliary variable for the dimensions
real(8) :: norm_r  ! For the norm of a real vector
real(8) :: trace_re  ! For the trace of a real matrix
real(8) :: matrix_average_sy  ! For the average of a real symmetric matrix, computed using a real vector: <x|A|x>

dda = da**2 - 1 ;   ddb = db**2 - 1
allocate(corrmat(1:dda,1:ddb)) ;   call corrmat_gellmann(da, db, rho, corrmat)  ! Computes the correlation matrix C

if ( ssys == 'a' ) then ! measurements over a
  allocate( Xi(1:dda,1:dda) ) ;   Xi = (4.d0/dble(da*da*db*db))*matmul(corrmat,transpose(corrmat)) ! Xi = TT^t = (4/da^2*db^2)*C*C^t
  deallocate(corrmat) ;   allocate(W(1:dda)) ;   call lapack_dsyevd('N', dda, Xi, W)
  if ( da == 2 ) then  ! exact expression
    allocate(rho_a(1:da,1:da)) ;   call partial_trace_b_he(rho, da, db, rho_a)
    allocate(bv(1:dda)) ;   call bloch_vector_gellmann(da, rho_a, bv) ;   deallocate(rho_a)  ! Computes the Bloch vector
    if ( norm_r(dda, bv) > 1.d-15 ) then 
      min_hs = trace_re(dda, Xi) - matrix_average_sy(dda, bv, Xi)/(norm_r(dda, bv)**2.d0)
    else
      min_hs = trace_re(dda, Xi) - min(W(1),W(2),W(3)) ;   deallocate(W)
    endif
  else if ( da > 2 ) then  ! lower bound
    min_hs = sum(W(da:dda)) ;   deallocate(W)
  endif 
  deallocate(Xi)
  
else if ( ssys == 'b' ) then  ! measurements over b
  allocate( Xi(1:ddb,1:ddb) ) ;   Xi = (4.d0/dble(da*da*db*db))*matmul(transpose(corrmat),corrmat) ! Xi = T^t*T = (4/da^2*db^2)*C^t*C
  deallocate(corrmat) ;   allocate(W(1:ddb)) ;   call lapack_dsyevd('N', ddb, Xi, W)
  if ( db == 2 ) then  ! exact expression
    allocate(rho_b(1:db,1:db)) ;   call partial_trace_a_he(rho, da, db, rho_b)
    allocate(bv(1:ddb)) ;   call bloch_vector_gellmann(db, rho_b, bv) ;   deallocate(rho_b)  ! Computes the Bloch vector
    if ( norm_r(ddb, bv) > 1.d-15 ) then 
      min_hs = trace_re(ddb, Xi) - matrix_average_sy(ddb, bv, Xi)/(norm_r(ddb, bv)**2.d0)
    else
      min_hs = trace_re(ddb, Xi) - min(W(1),W(2),W(3)) ;   deallocate(W)
    endif
  else if ( db > 2 ) then  ! lower bound
    min_hs = sum(W(db:ddb)) ;   deallocate(W)
  endif 
  deallocate(Xi)
endif

end
!###################################################################################################################################