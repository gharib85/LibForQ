!-----------------------------------------------------------------------------------------------------------------------------------
!                                   TWO-QUBIT evolved states for some quantum channels
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_ad(p, rho, rhop)  ! Returns the two-qubit evolved state for LOCAL AMPLITUDE DAMPING channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading",  Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2), K1d(1:2,1:2)  ! For the Kraus operators and their adjoints
complex(8) :: kp1(1:4,1:4), kp2(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = 1.d0 ;   K0(2,2) = dsqrt(1.d0-p)
K1 = 0.d0 ;   K1(1,2) = dsqrt(p) ;   K1d = 0.d0 ;   K1d(2,1) = K1(1,2)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp1) ;   rhop = rhop + matmul(matmul(kp1,rho),kp1)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp1) ;   call kronecker_product_c(K0, 2, 2, K1d, 2, 2, kp2) 
rhop = rhop + matmul(matmul(kp1,rho),kp2)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp1) ;   call kronecker_product_c(K1d, 2, 2, K0, 2, 2, kp2) 
rhop = rhop + matmul(matmul(kp1,rho),kp2)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp1) ;   call kronecker_product_c(K1d, 2, 2, K1d, 2, 2, kp2) 
rhop = rhop + matmul(matmul(kp1,rho),kp2)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_pd(p, rho, rhop)  ! Returns the two-qubit evolved state for local PHASE DAMPING channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading", Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2), K2(1:2,1:2)  ! For the Kraus operators
complex(8) :: kp(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = dsqrt(1.d0-p) ;   K0(2,2) = K0(1,1)
K1 = 0.d0 ;   K1(1,1) = dsqrt(p) 
K2 = 0.d0 ;   K2(2,2) = K1(1,1)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_bf(p, rho, rhop)  ! Returns the two-qubit evolved state for local BIT FLIP channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading",  Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2)  ! For the Kraus operators
complex(8) :: kp(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = dsqrt(1.d0-p) ;   K0(2,2) = K0(1,1)
K1 = 0.d0 ;   K1(1,2) = dsqrt(p) ;   K1(2,1) = K1(1,2)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_pf(p, rho, rhop)  ! Returns the two-qubit evolved state for local PHASE FLIP channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading", Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2)  ! For the Kraus operators
complex(8) :: kp(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = dsqrt(1.d0-p) ;   K0(2,2) = K0(1,1)
K1 = 0.d0 ;   K1(1,1) = dsqrt(p) ;   K1(2,2) = -K1(1,1)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_bpf(p, rho, rhop)  ! Returns the two-qubit evolved state for local BIT-PHASE FLIP channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading", Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2)  ! For the Kraus operators
complex(8) :: kp(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = dsqrt(1.d0-p) ;   K0(2,2) = K0(1,1)
K1 = 0.d0 ;   K1(2,1) = (0.d0,1.d0)*dsqrt(p) ;   K1(1,2) = -K1(2,1)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_d(p, rho, rhop)  ! Returns the two-qubit evolved state for local DEPOLARIZING channels
! Ref: M. B. Pozzobom and J. Maziero, "Environment-induced quantum coherence spreading",  Ann. Phys. 377, 243 (2017), arXiv:1605.04746
implicit none
real(8) :: p  ! Parametrized time
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)  ! Initial and evolved states
complex(8) :: K0(1:2,1:2), K1(1:2,1:2), K2(1:2,1:2), K3(1:2,1:2)  ! For the Kraus operators
complex(8) :: kp(1:4,1:4)  ! For the Kronecker product

K0 = 0.d0 ;   K0(1,1) = dsqrt(1.d0-3.d0*p/4.d0) ;   K0(2,2) = K0(1,1)
K1 = 0.d0 ;   K1(1,2) = dsqrt(p/4.d0) ;   K1(2,1) = K1(1,2)
K2 = 0.d0 ;   K2(2,1) = (0.d0,1.d0)*dsqrt(p/4.d0) ;   K2(1,2) = -K2(2,1)
K3 = 0.d0 ;   K3(1,1) = dsqrt(p/4.d0) ;   K3(2,2) = -K3(1,1)

rhop = 0.d0
call kronecker_product_c(K0, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K0, 2, 2, K3, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K1, 2, 2, K3, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K2, 2, 2, K3, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K3, 2, 2, K0, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K3, 2, 2, K1, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K3, 2, 2, K2, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)
call kronecker_product_c(K3, 2, 2, K3, 2, 2, kp) ;   rhop = rhop + matmul(matmul(kp,rho),kp)

end
!-----------------------------------------------------------------------------------------------------------------------------------