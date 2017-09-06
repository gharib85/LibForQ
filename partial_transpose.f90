!-----------------------------------------------------------------------------------------------------------------------------------
!                                                     PARTIAL TRANSPOSE
! Ref: J. Maziero, Computing partial transposes and related entanglement measures, Braz. J. Phys. 46, 605 (2016),  arXiv:1609.00323
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose(d, rho, rho_pt, nss, di, ssys)  ! Returns the partial transpose (PT) for general multipartite systems
implicit none
integer :: nss  ! Number of sub-systems
integer :: di(1:nss)  ! Vector specifying the dimensions of the sub-systems
integer :: ssys(1:nss)  ! Vector (with components equal to 0 or 1) specifying the sub-systems in which the PT is to be applied. 
                        ! If ssys(j)=0, then the PT is applied to the j-th subsystem. If ssys(j)=1 it is not.
integer :: d ! Total dimension (is the product of the sub-systems dimensions)
complex(8) :: rho(1:d,1:d), rho_pt(1:d,1:d)  ! Density matrix (given as input) and its partial transpose (the output)
complex(8), allocatable :: mat1(:,:), mat2(:,:), rhopt(:,:)  ! Auxiliary matrices
integer :: j, k, l  ! Auxiliary variables for counters
integer :: da, db, dc  ! Auxiliary variables for the dimensions

allocate( rhopt(1:d,1:d) )

! For bipartite systems
if ( nss == 2 ) then
  if ( ssys(1) == 0 ) then ; call partial_transpose_a(di(1), di(2), rho, rhopt)
  else if ( ssys(2) == 0 ) then ; call partial_transpose_b(di(1), di(2), rho, rhopt) ;   endif
! For multipartite systems
else if ( nss >= 3 ) then ;   allocate( mat1(1:d,1:d), mat2(1:d,1:d) ) 
  ! Left partial transposes
    l = 0 ;   do ;   if ( ssys(l+1) == 1 ) exit ;   l = l + 1 ;   enddo  ! l defines up to which position we shall apply the PT
    if ( l == 0 ) then ;   mat1 = rho  ! This matrix shall be used in the sequence if l = 0
    else if ( l > 0 ) then
      if ( l == 1 ) then ;   da = di(1) ;   else ;   da = product(di(1:l)) ;  endif ;   db = d/da
      call partial_transpose_a(da, db, rho, mat1)    
    endif
  ! Right partial transposes
    k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo  ! k defines up to which position we shall apply the PT
    if ( k == (nss+1) ) then ;   mat2 = mat1  ! This matrix shall be used in the sequence if k = nss+1
    else if ( k < (nss+1) ) then
      if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;  endif ;   da = d/db
      call partial_transpose_b(da, db, mat1, mat2)
    endif
  ! Inner partial transposes
    if ( (k-l) > 3 ) then  ! If (k-l) <= 3 there is no need to take inner partial transposes
      do j = (l+2), (k-2)
        if ( ssys(j) == 0 ) then ;   mat1 = mat2
          db = di(j) ;   da = product(di(1:j-1)) ;   dc = d/(da*db) ;   call partial_transpose_3(da, db, dc, mat1, mat2)
        endif
      enddo
    endif
    rhopt = mat2 ;   deallocate( mat1, mat2)
endif

rho_pt = rhopt ;   deallocate( rhopt )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_a(da, db, rho, rho_ta)  ! Returns its PT with relation to system A
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_ta(1:da*db,1:da*db)  ! Bipartite, original and transposed, states
integer :: ja, jb, ka, kb  ! Auxiliary variable for counters

forall ( ja = 1:da, ka = 1:da, jb = 1:db, kb = 1:db )
  rho_ta((ka-1)*db+jb,(ja-1)*db+kb) = rho((ja-1)*db+jb,(ka-1)*db+kb)
end forall

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_b(da, db, rho, rho_tb)  ! Returns the PT with relation to system B
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_tb(1:da*db,1:da*db)  ! Bipartite, original and transposed, states
integer :: ja, jb, ka, kb  ! Auxiliary variable for counters

forall ( ja = 1:da, ka = 1:da, jb = 1:db, kb = 1:db )
  rho_tb((ja-1)*db+kb,(ka-1)*db+jb) = rho((ja-1)*db+jb,(ka-1)*db+kb)
end forall

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_3(da, db, dc, rho, rho_tb)  ! Returns the PT with relation to system B, for a 3-partite state
implicit none
integer :: da, db, dc ! Dimensions of the subsystems
complex(8) :: rho(1:da*db*dc,1:da*db*dc), rho_tb(1:da*db*dc,1:da*db*dc)  ! Bipartite original and transposed states
integer :: ja, jb, jc, ka, kb, kc  ! Auxiliary variable for counters

forall ( ja = 1:da, ka = 1:da, jb = 1:db, kb = 1:db, jc = 1:dc, kc = 1:dc )
  rho_tb((ja-1)*db*dc+(kb-1)*dc+jc,(ka-1)*db*dc+(jb-1)*dc+kc) = rho((ja-1)*db*dc+(jb-1)*dc+jc,(ka-1)*db*dc+(kb-1)*dc+kc)
end forall

end
!-----------------------------------------------------------------------------------------------------------------------------------