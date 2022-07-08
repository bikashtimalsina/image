module poscar
implicit none
type lattice
real(8), dimension(3,3) :: direct
real(8), dimension(3,3) :: reciprocal
real(8) volume
real(8), dimension(3) :: alpha
end type lattice
contains
  subroutine findword(word,line,found,i)
  logical found
  character(*) line, word
  integer(4) i,l,k
  l = len_trim(line)
  k=len_trim(word)
  found = .False.
  do i=1,l-k+1
    if(line(i:i+k-1) .eq. word(1:k)) then
        found = .True.
        exit
    endif
  enddo
  if (line .eq. word) then
    found=.True.
  endif
  end subroutine findword
end module poscar

program image
use poscar
type(lattice) :: lat
integer(8) maxl, i, j, k, posread, reason, natom
real(8), allocatable:: position(:,:)
character line*199, word*59
logical found
posread=456
maxl=9999999
open(posread,file="./POSCAR")
read(posread,'(A)',iostat=reason) line
read(posread,'(A)',iostat=reason) line
read(posread,*,iostat=reason) lat%direct(1,1), lat%direct(1,2), lat%direct(1,3)
read(posread,*,iostat=reason) lat%direct(2,1), lat%direct(2,2), lat%direct(2,3)
read(posread,*,iostat=reason) lat%direct(3,1), lat%direct(3,2), lat%direct(3,3)
write(*,*) lat%direct(1,1), lat%direct(1,2), lat%direct(1,3)
write(*,*) lat%direct(2,1), lat%direct(2,2), lat%direct(2,3)
write(*,*) lat%direct(3,1), lat%direct(3,2), lat%direct(3,3)
read(posread,'(A)',iostat=reason) line
write(*,*) line
read(posread,'(A)',iostat=reason) line
write(*,*) line
read(posread,'(A)',iostat=reason) line
write(*,*) line
do i=1,maxl
read(posread,'(A)',iostat=reason) line
if (reason < 0) then
write(*,*) "END OF FILE"
exit
endif
enddo
close(posread)
end program image
