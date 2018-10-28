!-----------------------------------------------------------------------------------------------------------------------------------
subroutine bin2dec(vec_bin, nd, dec)
! Given a vector whose componets are the digits of a binary number, this function returns a base 10 integer.
! OBS. The components of the vector vary from the first digit of the base-two number (vec_bin(1)) to the last one (vec_bin(nd))
implicit none
integer :: nd  ! Number of digits
integer :: vec_bin(1:nd) ! Vector with the digits of the binary number
integer :: dec  ! The decimal corresponding to the binary number given as input
integer :: i  !  Auxiliary variable for counters

dec = 0
do i = 1, nd
  !dec = dec + vec_bin(i)*2**(nd - i)
  if (vec_bin(i) == 1) dec = dec + 2**(nd - i)
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine dec2bin(dec, nd, vec_bin)
! Given base 10 number, this subroutine returns a vector whose components are the digits of a binary number
implicit none
integer :: dec  ! Base 10 number to be converted to a binary number that is in its turn stored in a vector
integer :: nd  ! Number of digits needed to deal with the decimal numbers at play
integer :: vec_bin(1:nd)  ! Vector with the digits of the binary number
integer :: i  ! Auxiliary variable for counters

vec_bin = 0  ! initialization
do i = 1, nd
  if (2**(nd-i) <= dec) then ;   vec_bin(i) = 1 ;   dec = dec - 2**(nd-i) ;   end if
end do

end
!-----------------------------------------------------------------------------------------------------------------------------------
