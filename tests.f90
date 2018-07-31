program tests
 implicit none
 !call test_mat_sqrt()
 call test_werner()
end

subroutine test_mat_sqrt()
  implicit none
  complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
  complex(8) :: Asr(2,2)
  call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
  call mat_sqrt(2, sigma_0, Asr)
  call array_display(2, 2, Asr)
end


subroutine test_werner()
  implicit none
  complex(8) :: rho(4,4)
  real(8) :: w, discord_he
  w = 1.d0
  call  rho_bds(-w, -w, -w, rho)
  write(*,*) discord_he(2,2,rho)
end
