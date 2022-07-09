module poscar
implicit none
type lattice
real(8), dimension(3,3) :: direct
real(8), dimension(3,3) :: reciprocal
real(8) volume
real(8), dimension(3) :: a1
real(8), dimension(3) :: a2
real(8), dimension(3) :: a3
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
integer(4) p
integer(8) maxl, i, j, k, posread, reason, natom, atoms_len, &
 space_num, atom_type, atom_number, total_atom
 integer(8), allocatable :: atom_type_num(:)
real(8), allocatable:: position(:,:), pos_cart(:,:)
character line*199, word*59, space*1, atom_name*2
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
! evaluate the reciprocal lattice vector
lat%reciprocal(1,1)=lat%direct(2,2)*lat%direct(3,3)-lat%direct(3,2)*lat%direct(2,3)
lat%reciprocal(1,2)=lat%direct(3,1)*lat%direct(2,3)-lat%direct(2,1)*lat%direct(3,3)
lat%reciprocal(1,3)=lat%direct(2,1)*lat%direct(3,2)-lat%direct(3,1)*lat%direct(2,2)
!-------------------------------------------------------------------------------
lat%reciprocal(2,1)=lat%direct(3,2)*lat%direct(1,3)-lat%direct(3,3)*lat%direct(1,2)
lat%reciprocal(2,2)=lat%direct(1,1)*lat%direct(3,3)-lat%direct(3,1)*lat%direct(1,3)
lat%reciprocal(2,3)=lat%direct(3,1)*lat%direct(1,2)-lat%direct(1,1)*lat%direct(3,2)
!-------------------------------------------------------------------------------
lat%reciprocal(3,1)=lat%direct(1,2)*lat%direct(2,3)-lat%direct(2,2)*lat%direct(1,3)
lat%reciprocal(3,2)=lat%direct(2,1)*lat%direct(1,3)-lat%direct(1,1)*lat%direct(2,3)
lat%reciprocal(3,3)=lat%direct(1,1)*lat%direct(2,2)-lat%direct(2,1)*lat%direct(1,2)
lat%volume=lat%direct(1,1)*lat%reciprocal(1,1)+lat%direct(1,2)*lat%reciprocal(1,2)*&
lat%direct(1,3)*lat%reciprocal(1,3)
!-------------------------------------------------------------------------------
do i=1,3
do j=1,3
lat%reciprocal(i,j)=lat%reciprocal(i,j)/lat%volume
enddo
enddo
read(posread,'(A)',iostat=reason) line
space=' '
space_num=iachar(space)
atom_type=0
do i=1,len_trim(line)-1
if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
atom_name=line(i:i+1)
atom_type=atom_type+1
write(*,*) atom_name, atom_type
endif
enddo
allocate(atom_type_num(atom_type))
read(posread,*,iostat=reason) atom_type_num
write(*,*) atom_type_num
total_atom=0
do i=1,atom_type
total_atom=total_atom+atom_type_num(i)
enddo
read(posread,'(A)',iostat=reason) line
write(*,*) line
allocate(position(total_atom,3))
do i=1,maxl
read(posread,*,iostat=reason) position(i,:)
write(*,*) position(i,:)
if (reason < 0) then
write(*,*) "END OF FILE"
exit
endif
enddo
close(posread)
! Check for direct coordinates
word='D'
call findword(word,line,found,p)
if ( found ) then
write(*,*) "Found D"
allocate(pos_cart(total_atom,3))
pos_cart=0.0d0
do i=1,total_atom
do j=1,3
pos_cart(i,1)=pos_cart(i,1)+position(i,j)*lat%direct(j,1)
pos_cart(i,2)=pos_cart(i,2)+position(i,j)*lat%direct(j,2)
pos_cart(i,3)=pos_cart(i,3)+position(i,j)*lat%direct(j,3)
enddo
enddo
endif
! Should be noted it could be either D or d for Direct
word='d'
call findword(word,line,found,p)
if ( found ) then
write(*,*) "Found D"
allocate(pos_cart(total_atom,3))
pos_cart=0.0d0
do i=1,total_atom
do j=1,3
pos_cart(i,1)=pos_cart(i,1)+position(i,j)*lat%direct(j,1)
pos_cart(i,2)=pos_cart(i,2)+position(i,j)*lat%direct(j,2)
pos_cart(i,3)=pos_cart(i,3)+position(i,j)*lat%direct(j,3)
enddo
enddo
endif
do i=1,total_atom
write(*,*) pos_cart(i,1), pos_cart(i,2), pos_cart(i,3)
enddo
deallocate(atom_type_num,position,pos_cart)
end program image
