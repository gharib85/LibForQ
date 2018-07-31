!-----------------------------------------------------------------------------------------------------------------------------------
!                                           Hamiltonians
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine Ising_1D(ns, h, ham) ! Returns the Hamiltonian for the 1D ISING chain in a transverse field (J := 1):
! H = -0.5*J*\sum_{j}\sigma_{j}^{z}\sigma_{j+1}^{z} - h*\sum_{j}\sigma_{j}^{z}
implicit none
real(8) :: h  ! Transverse field
integer :: ns  ! Number of spins
integer :: d  ! For the Hamiltonian's dimension (d=2^ns)
complex(8) :: ham(1:2**ns,1:2**ns)  ! The Hamiltonian to be returned
complex(8), allocatable :: ham_J(:,:)  ! Hamiltonian terms involving exchange interactions in the z direction
complex(8), allocatable :: ham_h(:,:)  ! Hamiltonian terms involving the transverse field in the x direction
integer, allocatable :: ord_pm(:)  ! Order of the Pauli matrices for the tensor products
complex(8), allocatable :: kp_pauli_mat(:,:)  ! For the Kronecker product among several Pauli matrices
integer :: j, k ! Auxiliary variables for counters

d = 2**ns
allocate ( ham_J(1:d,1:d), ham_h(1:d,1:d), kp_pauli_mat(1:d,1:d), ord_pm(1:ns) )

! Summing up the hamiltonian terms -J*s_{j}^{z}s_{k}^{z} (J := 1, exchange interaction in the z=3 direction)
ham_J = 0.d0 ;   ord_pm = 0
do j = 1, ns-1
  if ( j == 1 ) ord_pm(j) = 3 ;   ord_pm(j+1) = 3
  call kron_prod_pauli_mat(ord_pm, ns, kp_pauli_mat) ;   forall ( k = 1:d ) ham_J(k,k) = ham_J(k,k) + kp_pauli_mat(k,k)
  ord_pm(j) = 0
enddo
ham_J = -0.5d0*ham_J

! Summing up the Hamiltonian terms -h*s_{j}^{x} (magnetic field in the x=1 direction)
ham_h = 0.d0 ;   ord_pm = 0
do j = 1, ns
  ord_pm(j) = 1 ;   call kron_prod_pauli_mat(ord_pm, ns, kp_pauli_mat) ;   ham_h = ham_h + kp_pauli_mat
  ord_pm(j) = 0
enddo
ham_h = -h*ham_h

ham = ham_J + ham_h

deallocate ( ham_J, ham_h, ord_pm, kp_pauli_mat )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine XY_1D(ns, h, g, ham) ! Returns the Hamiltonian for the 1D XY chain (the exchange energy is set to unity: J=1):
!H = -0.5*J(1+g)*\sum_{j}\sigma_{j}^{x}\sigma_{j+1}^{x} -0.5*J(1-g)*\sum_{j}\sigma_{j}^{y}\sigma_{j+1}^{y} -h*\sum_{j}\sigma_{j}^{z}
implicit none
real(8) :: h  ! Transverse field
real(8) :: g  ! Anisotropy
integer :: ns  ! Number of spins
integer :: d  ! For the Hamiltonian's dimension (d=2^ns)
complex(8) :: ham(1:2**ns,1:2**ns)  ! The Hamiltonian to be returned
complex(8), allocatable :: ham_Jp(:,:), ham_Jm(:,:)  ! Hamiltonian terms involving exchange interactions in the x and y directions
complex(8), allocatable :: ham_h(:,:)  ! Hamiltonian terms involving the transverse field in the z direction
integer, allocatable :: ord_pm(:)  ! Order of the Pauli matrices for the tensor products
complex(8), allocatable :: kp_pauli_mat(:,:)  ! For the Kronecker product among several Pauli matrices
integer :: j, k ! Auxiliary variables for counters

allocate( ham_Jp(1:2**ns,1:2**ns), ham_Jm(1:2**ns,1:2**ns), ham_h(1:2**ns,1:2**ns), kp_pauli_mat(1:2**ns,1:2**ns), ord_pm(1:ns) )

! Summing up the hamiltonian terms -J*(1+g)*s_{j}^{x}s_{k}^{x}
ham_Jp = 0.d0 ;   ord_pm = 0
do j = 1, ns-1
  if ( j == 1 ) ord_pm(j) = 1 ;   ord_pm(j+1) = 1
  call kron_prod_pauli_mat(ord_pm, ns, kp_pauli_mat) ;   ham_Jp = ham_Jp + kp_pauli_mat
  ord_pm(j) = 0
enddo
ham_Jp = -0.5d0*(1.d0+g)*ham_Jp

! Summing up the hamiltonian terms -J*(1+g)*s_{j}^{x}s_{k}^{x}
ham_Jm = 0.d0 ;   ord_pm = 0
do j = 1, ns-1
  if ( j == 1 ) ord_pm(j) = 2 ;   ord_pm(j+1) = 2
  call kron_prod_pauli_mat(ord_pm, ns, kp_pauli_mat) ;   ham_Jm = ham_Jm + kp_pauli_mat
  ord_pm(j) = 0
enddo
ham_Jm = -0.5d0*(1.d0-g)*ham_Jm

! Summing up the Hamiltonian terms -h*s_{j}^{x} (magnetic field in the x=1 direction)
ham_h = 0.d0 ;   ord_pm = 0
do j = 1, ns
  ord_pm(j) = 3
  call kron_prod_pauli_mat(ord_pm, ns, kp_pauli_mat) ;   ham_h = ham_h + kp_pauli_mat
  ord_pm(j) = 0
enddo
ham_h = -h*ham_h

ham = ham_Jp + ham_Jm + ham_h

deallocate( ham_Jp, ham_Jm, ham_h, ord_pm, kp_pauli_mat )

end
!-----------------------------------------------------------------------------------------------------------------------------------
