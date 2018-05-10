!-----------------------------------------------------------------------------------------------------------------------------------
program dipolar
  !call Et()
  !call Emax()
  !call Et_rho()
  !call Emax_rho()
  !call EEpsi()
  call EErho()
end
!-----------------------------------------------------------------------------------------------------------------------------------
! pure-prouct initial state, changing t
subroutine Et()
implicit none
real(8) :: D, delta, tta, ttb, t, E, C, dt, dtt, concurrence_2qb
complex(8) :: psia(2), psib(2), psiI(4), psiF(4), rho(4,4), U(4,4)
complex(8) :: psip(4), psim(4), phip(4), phim(4), proj1(4,4), proj2(4,4), proj3(4,4), proj4(4,4)
real(8) :: pi
open(unit = 13, file = 'Et_tbpi2.dat', status = 'unknown')

pi = 4.d0*datan(1.d0);  D = 1.d0
!delta = -2.d0
call bell_basis(psip, psim, phip, phim)
call projector(psim, 4, proj1);  call projector(psip, 4, proj2)
call projector(phim, 4, proj3);  call projector(phip, 4, proj4)
!call psi_qubit(pi/2.d0, pi/2.d0, psib);  write(*,*) psib(1), psib(2); stop

dtt = 0.01d0;  dt = dtt
ttb = pi/2.d0;  call psi_qubit(ttb, 0.d0, psib)
tta = 0.d0 - dtt
do
  tta = tta + dtt;  if(tta > pi) exit
  call psi_qubit(tta, 0.d0, psia);  call kronecker_product_c(psia, 2, 1, psib, 2, 1, psiI)
  !write(*,*) dble(psiI(1)), dble(psiI(2)), dble(psiI(3)), dble(psiI(4)); stop
  t = 0.d0 - dt
  do
    t = t + dt;  if(t > pi) exit
    U = proj1!*(cos(D*(delta+2.d0)*t/2.d0) + (0.d0,1.d0)*sin(D*(delta+2.d0)*t/2.d0))    ! the unitary operator
    U = U + (cos(-2.d0*D*t) - (0.d0,1.d0)*sin(2.d0*D*t))*proj2
    U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj3
    U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj4
    psiF = matmul(U,psiI)!;  write(*,*) dble(psiF(1)), dble(psiF(2)), dble(psiF(3)), dble(psiF(4)); stop
    call projector(psiF, 4, rho)
    !write(*,*) dble(rho(1,1)), dble(rho(1,2)), dble(rho(1,3)), dble(rho(1,4))
    !write(*,*) dble(rho(2,1)), dble(rho(2,2)), dble(rho(2,3)), dble(rho(2,4))
    !write(*,*) dble(rho(3,1)), dble(rho(3,2)), dble(rho(3,3)), dble(rho(3,4))
    !write(*,*) dble(rho(4,1)), dble(rho(4,2)), dble(rho(4,3)), dble(rho(4,4)); stop
    E = concurrence_2qb(rho)
    if(E > 1.d-8) write(13,*) t, tta, E
  enddo
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
! pure-prouct initial state, t fixed
subroutine Emax()
implicit none
real(8) :: D, delta, tta, ttb, t, E, C, dt, dtt, concurrence_2qb, tm, Em
complex(8) :: psia(2), psib(2), psiI(4), psiF(4), rho(4,4), U(4,4)
complex(8) :: psip(4), psim(4), phip(4), phim(4), proj1(4,4)
complex(8) :: proj2(4,4), proj3(4,4), proj4(4,4)
real(8) :: pi
open(unit = 13, file = 'Emax_D2.dat', status = 'unknown')

pi = 4.d0*datan(1.d0);  delta = -2.d0
call bell_basis(psip, psim, phip, phim)
call projector(psim, 4, proj1);  call projector(psip, 4, proj2)
call projector(phim, 4, proj3);  call projector(phip, 4, proj4)
D = 2.d0
dtt = 0.05d0;  dt = 0.005d0
ttb = 0.d0 - dtt
do
  ttb = ttb + dtt;  call psi_qubit(ttb, 0.d0, psib)
  tta = 0.d0 - dtt
  do
    tta = tta + dtt
    call psi_qubit(tta, 0.d0, psia);  call kronecker_product_c(psia, 2, 1, psib, 2, 1, psiI)
    tm = 0.d0;  Em = 0.d0
    t = 0.d0 - dt
    do
      t = t + dt
      U = proj1!*(cos(D*(delta+2.d0)*t/2.d0) + (0.d0,1.d0)*sin(D*(delta+2.d0)*t/2.d0))    ! the unitary operator
      U = U + (cos(-2.d0*D*t) + (0.d0,1.d0)*sin(-2.d0*D*t))*proj2
      U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj3
      U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj4
      psiF = matmul(U,psiI);  call projector(psiF, 4, rho);  E = concurrence_2qb(rho)
      if(E > Em)then;  Em = E;  tm = t; endif
      if(t > pi) exit
    enddo
    write(13,*) tta, ttb, Em, tm
    if(tta > pi) exit
  enddo
  if(ttb > pi) exit
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
! mixed-prouct initial state, changing t
subroutine Et_rho()
implicit none
real(8) :: D, a, b, t, Ec, dt, dtt, concurrence_2qb, Dis
real(8) :: discord_hs_2qb, negativity, pi
complex(8) :: rhoa(2,2), rhob(2,2), rhoI(4,4), rhot(4,4), U(4,4), Ua(4,4)
complex(8) :: psip(4), psim(4), phip(4), phim(4), proj1(4,4), proj2(4,4)
complex (8) :: proj3(4,4), proj4(4,4), rho_pt(4,4)
open(unit = 13, &
     file = '/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/ERhor1b1.dat', &
     status = 'unknown')

pi = 4.d0*datan(1.d0);
D = 1.d0  ! Distance-related parameter
call bell_basis(psip, psim, phip, phim)
call projector(psim, 4, proj1);  call projector(psip, 4, proj2)
call projector(phim, 4, proj3);  call projector(phip, 4, proj4)
dtt = 2.d0/200.d0 - 1.d-10;  dt = 1.d0/350.d0 - 1.d-10
b = 1.d0
call rho_qubit(b,0.d0,0.d0,rhob)
!call rho_qubit(0.d0,0.d0,b,rhob)
a = -1.d0 - dtt
do
  a = a + dtt;  if (a > 1.d0) exit
  t = 0.d0 - dt
  call rho_qubit(a,0.d0,0.d0,rhoa)
  call kronecker_product_c(rhoa,2,2,rhob,2,2,rhoI)  ! system initial state
  do
    t = t + dt;  if(t > pi) exit
    ! the unitary operator
    U = proj1!*(cos(D*(delta+2.d0)*t/2.d0) + (0.d0,1.d0)*sin(D*(delta+2.d0)*t/2.d0))
    U = U + (cos(-2.d0*D*t) + (0.d0,1.d0)*sin(-2.d0*D*t))*proj2
    U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj3
    U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj4
    call adjoint(4, 4, U, Ua)
    rhot = matmul(matmul(U,rhoI),Ua)  ! system evolved state
    !call partial_transpose_a(2, 2, rhot, rho_pt);  E =  negativity(4, rho_pt)
    Ec = concurrence_2qb(rhot);
    !Dis = discord_hs_2qb('a', rhot)
    if (Ec > 1.d-12) write(13,*) a, t, Ec!, Dis
  enddo
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
! mixed-product initial state, t fixed
subroutine Emax_rho()
implicit none
real(8) :: D, a, b, t, E, C, dt, dtt, concurrence_2qb, tm, Em!, Dis, Dm, discord_hs_2qb
complex(8) :: rhoa(2,2), rhob(2,2), rhoI(4,4), rhot(4,4), U(4,4), Ua(4,4)
complex(8) :: psip(4), psim(4), phip(4), phim(4)
complex(8) :: proj1(4,4), proj2(4,4), proj3(4,4), proj4(4,4)
real(8) :: pi
open(unit = 13, &
    file = '/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/ERho3tpi4.dat', &
     status = 'unknown')

pi = 4.d0*datan(1.d0)
!delta = -2.d0
D = 1.d0
call bell_basis(psip, psim, phip, phim)
call projector(psim, 4, proj1);  call projector(psip, 4, proj2)
call projector(phim, 4, proj3);  call projector(phip, 4, proj4)
dtt = 2.d0/400.d0 - 1.d-10
b = -1.d0 - dtt
dt = 0.001d0
do
  b = b + dtt;  if(b > 1.d0) exit
  call rho_qubit(0.d0,0.d0,b,rhob)
  a = -1.d0 - dtt
  do
    a = a + dtt;  if(a > 1.d0) exit
    call rho_qubit(0.d0,0.d0,a,rhoa)
    call kronecker_product_c(rhoa, 2, 2, rhob, 2, 2, rhoI)
    !tm = 0.d0;
    !Em = 0.d0!;  Dm = 0.d0
    !t = 0.d0 - dt
    !do
      t = 1.d0*pi/4.d0
      !t = t + dt;  if(t > pi) exit
      U = proj1!*(cos(D*(delta+2.d0)*t/2.d0) + (0.d0,1.d0)*sin(D*(delta+2.d0)*t/2.d0))
      U = U + (cos(-2.d0*D*t) + (0.d0,1.d0)*sin(-2.d0*D*t))*proj2
      U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj3
      U = U + (cos(D*t) + (0.d0,1.d0)*sin(D*t))*proj4
      call adjoint(4, 4, U, Ua)
      rhot = matmul(matmul(U,rhoI),Ua)  ! system evolved state
      E = concurrence_2qb(rhot)
      !if (E > Em) then
      !  Em = E
        !tm = t
      !endif
      !Dis = discord_hs_2qb('a', rhot)
      !if (Dis > Dm) then;  Dm = Dis;  endif
    !enddo
    write(13,*) a, b, E
  enddo
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
! partially entangle initial state, pure
subroutine EEpsi()
implicit none
real(8) :: w, D, t, dt, dw, pi, E
open(unit = 13, &
    file = '/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/EEpsi.dat', &
     status = 'unknown')

pi = 4.d0*datan(1.d0)
D = 1.d0
dw = 1.d0/200.d0
dt = pi/200.d0

w = -dw
do
  w = w + dw;  if(w > 1.d0) exit
  t = -dt
  do
    t = t + dt;  if(t > pi) exit
    E = sqrt((sin(2.d0*D*t))**2.d0 + 4.d0*w*(1.d0-w)*(cos(2.d0*D*t))**2.d0)
    write(13,*) w, t, E
  enddo
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
! partially entangle initial state, mixed
subroutine EErho()
implicit none
real(8) :: p, dp, w, D, t, dt, dw, pi, E, concurrence_2qb
complex(8) :: rho(4,4), id(4,4), proj(4,4), vec(4)
open(unit = 13, &
    file = '/home/jonasmaziero/Dropbox/Research/qnesses/interplay/dipolar/dipolarCalc/EErhotpi2.dat', &
     status = 'unknown')

pi = 4.d0*atan(1.d0)
D = 1.d0
t = pi/2.d0
dw = 1.d0/101.d0
dp = 1.d0/101.d0

w = -dw
do
  w = w + dw;  if(w > 1.d0) exit
  p = -dp
  do
    p = p + dp;  if(p > 1.d0) exit
    call identity_c(4, id)
    vec(1) = 0; vec(2) = sqrt(w)*cos(D*t) - (0.d0,1.d0)*sqrt(1.d0-w)*sin(D*t)
    vec(3) = sqrt(1.d0-w)*cos(D*t) - (0.d0,1.d0)*sqrt(w)*sin(D*t); vec(4) = 0.d0
    call projector(vec, 4, proj)
    rho = ((1.d0-p)/4.d0)*id + p*proj
    E = concurrence_2qb(rho)
    write(13,*) w, p, E
  enddo
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
