module poscar
implicit none
type lattice
real(8), dimension(3,3) :: direct
real(8), dimension(3,3) :: reciprocal
real(8) volume
real(8), dimension(3) :: a
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
 space_num, atom_type, atom_number, total_atom, outread
 integer(8), allocatable :: atom_type_num(:)
real(8), allocatable:: position(:,:), pos_cart(:,:), pos_outcar(:,:), force_outcar(:,:)
real(8), allocatable :: pos_outcar_new(:,:), rij(:,:)
character line*199, word*59, space*1, atom_name*2
logical found
posread=456
outread=789
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
write(*,*) "Found Direct/direct coordinates"
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
write(*,*) "Found Direct/direct coordinates"
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
! Read OUTCAR file beyond this point
allocate(pos_outcar(total_atom,3),force_outcar(total_atom,3))
allocate(pos_outcar_new(total_atom,3))
allocate(rij(total_atom,3))
open(outread,file="./OUTCAR")
do i=1,maxl
read(outread,'(A)',iostat=reason) line
word=' POSITION'
call findword(word,line,found,p)
if(found) then
read(outread,'(A)',iostat=reason) line
do j=1,total_atom
read(outread,*,iostat=reason) pos_outcar(j,1),pos_outcar(j,2), &
pos_outcar(j,3), force_outcar(j,1), force_outcar(j,2), force_outcar(j,3)
!write(*,*) pos_outcar(j,1),pos_outcar(j,2),pos_outcar(j,3), &
!force_outcar(j,1),force_outcar(j,2),force_outcar(j,3)
enddo
endif
if(reason < 0) then
write(*,*) ".......END OF OUTCAR......."
exit
endif
enddo
! Get the distance of each atom from POSCAR and OUTCAR file
do i=1,total_atom
do j=1,3
rij(i,j)=pos_outcar(i,j)-pos_cart(i,j)
enddo
enddo
! Try to change the cartesian coordinates from OUTCAR to reduce coordinates
write(*,*) "--------------------------------------------------------------------"
do i=1,total_atom
lat%a=0.0d0
pos_outcar_new(i,1)=0.0d0
pos_outcar_new(i,2)=0.0d0
pos_outcar_new(i,3)=0.0d0
do j=1,3
lat%a(1)=lat%a(1)+pos_outcar(i,j)*lat%reciprocal(1,j)
lat%a(2)=lat%a(2)+pos_outcar(i,j)*lat%reciprocal(2,j)
lat%a(3)=lat%a(3)+pos_outcar(i,j)*lat%reciprocal(3,j)
enddo
write(*,*) lat%a(1),lat%a(2),lat%a(3)
lat%a(1)=lat%a(1)-nint(lat%a(1))
lat%a(2)=lat%a(2)-nint(lat%a(2))
lat%a(3)=lat%a(3)-nint(lat%a(3))
write(*,*) lat%a(1),lat%a(2),lat%a(3)
!lat%a(1)=anint(lat%a(1))
!lat%a(2)=anint(lat%a(2))
!lat%a(3)=anint(lat%a(3))
write(*,*) position(i,1), position(i,2), position(i,3)
write(*,*) "-------------------------------------------------------------------"
do j=1,3
pos_outcar_new(i,j)=lat%a(1)*lat%direct(1,j)+lat%a(2)*lat%direct(2,j)+ &
lat%a(3)*lat%direct(3,j)
enddo
!write(*,*) lat%a(1),lat%a(2),lat%a(3)
enddo
close(outread)
write(*,*) "--------------------------------------------------------------------"
do i=1,total_atom
write(*,*) pos_cart(i,1),pos_cart(i,2),pos_cart(i,3)
write(*,*) pos_outcar(i,1),pos_outcar(i,2),pos_outcar(i,3)
write(*,*) pos_cart(i,1)-pos_outcar(i,1),pos_cart(i,2)-pos_outcar(i,2), &
pos_cart(i,3)-pos_outcar(i,3)
write(*,*) "--------------------------------------------------------------------"
enddo
deallocate(atom_type_num,position,pos_cart)
deallocate(pos_outcar,force_outcar,pos_outcar_new)
end program image
