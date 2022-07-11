module poscar
implicit none
type lattice
real(8), dimension(3,3) :: direct
real(8), dimension(3,3) :: reciprocal
character(8), dimension(:), allocatable :: atoms_type
integer(8), dimension(:), allocatable :: atoms_num
integer(8) total_atom
real(8), allocatable :: position(:,:)
real(8) volume
end type lattice
! user defined constructor
interface lattice
  module procedure get_structure
end interface lattice
contains
  ! Routine to find a word, it could be anything
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
  ! Routine to get header information from POSCAR and extract lattice size, the number of atoms and type_of_atoms
  function get_structure(filename) result(self)
    type(lattice) :: self
    character(8), intent(in) :: filename*6
    character(8) line*199, word*19
    character space*1, atom_name*2
    logical found_capD, found_smallD
    integer(4) posread, i, j, k, p, n, space_num, atom_type, reason
    posread=123
    open(posread,file=filename)
      read(posread,'(A)') line
      read(posread,'(A)') line
      read(posread,*) self%direct(1,1),self%direct(1,2),self%direct(1,3)
      read(posread,*) self%direct(2,1),self%direct(2,2),self%direct(2,3)
      read(posread,*) self%direct(3,1),self%direct(3,2),self%direct(3,3)
      ! get the reciprocal lattice vector
      self%reciprocal(1,1)=self%direct(2,2)*self%direct(3,3)-self%direct(3,2)*self%direct(2,3)
      self%reciprocal(1,2)=self%direct(3,1)*self%direct(2,3)-self%direct(2,1)*self%direct(3,3)
      self%reciprocal(1,3)=self%direct(2,1)*self%direct(3,2)-self%direct(3,1)*self%direct(2,2)
      !-------------------------------------------------------------------------------
      self%reciprocal(2,1)=self%direct(3,2)*self%direct(1,3)-self%direct(3,3)*self%direct(1,2)
      self%reciprocal(2,2)=self%direct(1,1)*self%direct(3,3)-self%direct(3,1)*self%direct(1,3)
      self%reciprocal(2,3)=self%direct(3,1)*self%direct(1,2)-self%direct(1,1)*self%direct(3,2)
      !-------------------------------------------------------------------------------
      self%reciprocal(3,1)=self%direct(1,2)*self%direct(2,3)-self%direct(2,2)*self%direct(1,3)
      self%reciprocal(3,2)=self%direct(2,1)*self%direct(1,3)-self%direct(1,1)*self%direct(2,3)
      self%reciprocal(3,3)=self%direct(1,1)*self%direct(2,2)-self%direct(2,1)*self%direct(1,2)
      self%volume=self%direct(1,1)*self%reciprocal(1,1)+self%direct(1,2)*self%reciprocal(1,2)*&
      self%direct(1,3)*self%reciprocal(1,3)
      !-------------------------------------------------------------------------------
      do i=1,3
      do j=1,3
      self%reciprocal(i,j)=self%reciprocal(i,j)/self%volume
      enddo
      enddo
      read(posread,'(A)',iostat=reason) line
      write(*,*) line
      space=' '
      space_num=iachar(space)
      atom_type=0
      do i=1,len_trim(line)-1
      if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
      atom_type=atom_type+1
      endif
      enddo
      allocate(self%atoms_type(atom_type))
      allocate(self%atoms_num(atom_type))
      k=0
      do i=1,len_trim(line)-1
      if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
      atom_name=line(i:i+1)
      atom_type=atom_type+1
      k=k+1
      self%atoms_type(k)=atom_name
      endif
      enddo
      read(posread,*,iostat=reason) self%atoms_num(:)
      self%total_atom=sum(self%atoms_num)
      read(posread,*,iostat=reason) line
      word='D'
      call findword(word,line,found_capD,n)
      word='d'
      call findword(word,line,found_smallD,n)
      if (.not. found_capD .and. .not. found_smallD) then
        allocate(self%position(self%total_atom,3))
        do i=1,self%total_atom
          if (i .le. self%total_atom) then
            read(posread,*,iostat=reason) self%position(i,:)
          endif
          if (i .gt. self%total_atom .and. reason < 0) then
            exit
          endif
        enddo
      endif
      if(found_capD .or. found_smallD) then
        allocate(self%position(self%total_atom,3))
        do i=1,self%total_atom
          if (i .le. self%total_atom) then
            read(posread,*,iostat=reason) self%position(i,:)
          endif
          ! this is for the first coordinate
          self%position(i,1)=self%position(i,1)*self%direct(1,1)+ &
          self%position(i,1)*self%direct(2,1)+ &
          self%position(i,1)*self%direct(3,1)
          ! this is for the second coordinate
          self%position(i,2)=self%position(i,2)*self%direct(1,2)+ &
          self%position(i,2)*self%direct(2,2)+ &
          self%position(i,2)*self%direct(3,2)
          ! this is for the third coordinate
          self%position(i,3)=self%position(i,3)*self%direct(1,3)+ &
          self%position(i,3)*self%direct(2,3)+ &
          self%position(i,3)*self%direct(3,3)
          if (i .gt. self%total_atom .and. reason < 0) then
            exit
          endif
        enddo
      endif
    close(posread)
  end function get_structure
end module poscar

program image
use poscar
type(lattice) lat
integer(4) num_args, p, i
character(8) word*6
character(8), dimension(:), allocatable :: args
logical found
num_args=command_argument_count()
if (num_args .ne. 2) then
  write(*,*) "Please provide POSCAR and OUTCAR file during execution"
endif
if (num_args .eq. 2) then
  allocate(args(num_args))
  do i=1,num_args
    call get_command_argument(i,args(i))
    !Check if POSCAR is present and read it
    if(i .eq. 1) then
      word='POSCAR'
      call findword(word,args(i),found,p)
      if (.not. found) then
        write(*,*) "Provide POSCAR first as an argument"
      endif
      if (found) then
        inquire(file=args(i),exist=found)
        if (.not. found) then
          write(*,*) "POSCAR NOT FOUND IN THE CURRENT DIRECTORY"
        endif
        if (found) then
          write(*,*) "POSCAR FOUND IN THE CURRENT DIRECTORY"
          lat=get_structure(args(i))
          deallocate(lat%atoms_type)
          deallocate(lat%atoms_num)
          deallocate(lat%position)
        endif
      endif
    endif
    !Check if OUTCAR is present and read it
    if(i .eq. 2) then
      word='OUTCAR'
      call findword(word,args(i),found,p)
      if (.not. found) then
        write(*,*) "Provide OUTCAR as second argument"
      endif
      if (found) then
        inquire(file=args(i),exist=found)
        if (.not. found) then
          write(*,*) "OUTCAR NOT FOUND IN THE CURRENT DIRECTORY"
        endif
        if (found) then
          write(*,*) "OUTCAR FOUND IN THE CURRENT DIRECTORY"
        endif
      endif
    endif
    write(*,*) "FILE PROVIDED ARE: ",args(i)
  enddo
  deallocate(args)
endif
! type(lattice) :: lat
! integer(4) p
! integer(8) maxl, i, j, k, posread, reason, natom, atoms_len, &
!  space_num, atom_type, atom_number, total_atom, outread
!  integer(8), allocatable :: atom_type_num(:)
! real(8), allocatable:: position(:,:), pos_cart(:,:), pos_outcar(:,:), force_outcar(:,:)
! real(8), allocatable :: pos_outcar_new(:,:), rij(:,:)
! character line*199, word*59, space*1, atom_name*2
! logical found
! posread=456
! outread=789
! maxl=9999999
! open(posread,file="./POSCAR")
! read(posread,'(A)',iostat=reason) line
! read(posread,'(A)',iostat=reason) line
! read(posread,*,iostat=reason) lat%direct(1,1), lat%direct(1,2), lat%direct(1,3)
! read(posread,*,iostat=reason) lat%direct(2,1), lat%direct(2,2), lat%direct(2,3)
! read(posread,*,iostat=reason) lat%direct(3,1), lat%direct(3,2), lat%direct(3,3)
! write(*,*) lat%direct(1,1), lat%direct(1,2), lat%direct(1,3)
! write(*,*) lat%direct(2,1), lat%direct(2,2), lat%direct(2,3)
! write(*,*) lat%direct(3,1), lat%direct(3,2), lat%direct(3,3)
! ! evaluate the reciprocal lattice vector
! lat%reciprocal(1,1)=lat%direct(2,2)*lat%direct(3,3)-lat%direct(3,2)*lat%direct(2,3)
! lat%reciprocal(1,2)=lat%direct(3,1)*lat%direct(2,3)-lat%direct(2,1)*lat%direct(3,3)
! lat%reciprocal(1,3)=lat%direct(2,1)*lat%direct(3,2)-lat%direct(3,1)*lat%direct(2,2)
! !-------------------------------------------------------------------------------
! lat%reciprocal(2,1)=lat%direct(3,2)*lat%direct(1,3)-lat%direct(3,3)*lat%direct(1,2)
! lat%reciprocal(2,2)=lat%direct(1,1)*lat%direct(3,3)-lat%direct(3,1)*lat%direct(1,3)
! lat%reciprocal(2,3)=lat%direct(3,1)*lat%direct(1,2)-lat%direct(1,1)*lat%direct(3,2)
! !-------------------------------------------------------------------------------
! lat%reciprocal(3,1)=lat%direct(1,2)*lat%direct(2,3)-lat%direct(2,2)*lat%direct(1,3)
! lat%reciprocal(3,2)=lat%direct(2,1)*lat%direct(1,3)-lat%direct(1,1)*lat%direct(2,3)
! lat%reciprocal(3,3)=lat%direct(1,1)*lat%direct(2,2)-lat%direct(2,1)*lat%direct(1,2)
! lat%volume=lat%direct(1,1)*lat%reciprocal(1,1)+lat%direct(1,2)*lat%reciprocal(1,2)*&
! lat%direct(1,3)*lat%reciprocal(1,3)
! !-------------------------------------------------------------------------------
! do i=1,3
! do j=1,3
! lat%reciprocal(i,j)=lat%reciprocal(i,j)/lat%volume
! enddo
! enddo
! read(posread,'(A)',iostat=reason) line
! write(*,*) line
! space=' '
! space_num=iachar(space)
! atom_type=0
! do i=1,len_trim(line)-1
! if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
! atom_name=line(i:i+1)
! atom_type=atom_type+1
! write(*,*) atom_name, atom_type
! endif
! enddo
! allocate(atom_type_num(atom_type))
! read(posread,*,iostat=reason) atom_type_num(:)
! write(*,*) atom_type_num
! write(*,*) shape(atom_type_num)
! total_atom=0
! do i=1,atom_type
! total_atom=total_atom+atom_type_num(i)
! enddo
! write(*,*) "total atom is: ",total_atom
! read(posread,'(A)',iostat=reason) line
! write(*,*) line
! allocate(position(total_atom,3))
! do i=1,total_atom
! if (reason < 0) then
! write(*,*) "END OF FILE"
! exit
! endif
! read(posread,'(A)',iostat=reason) position(i,1),position(i,2),position(i,3)
! enddo
! close(posread)
! ! Check for direct coordinates
! word='D'
! call findword(word,line,found,p)
! if ( found ) then
! write(*,*) "Found Direct/direct coordinates"
! allocate(pos_cart(total_atom,3))
! pos_cart=0.0d0
! do i=1,total_atom
! do j=1,3
! pos_cart(i,1)=pos_cart(i,1)+position(i,j)*lat%direct(j,1)
! pos_cart(i,2)=pos_cart(i,2)+position(i,j)*lat%direct(j,2)
! pos_cart(i,3)=pos_cart(i,3)+position(i,j)*lat%direct(j,3)
! enddo
! enddo
! endif
! ! Should be noted it could be either D or d for Direct
! word='d'
! call findword(word,line,found,p)
! if ( found ) then
! write(*,*) "Found Direct/direct coordinates"
! allocate(pos_cart(total_atom,3))
! pos_cart=0.0d0
! do i=1,total_atom
! do j=1,3
! pos_cart(i,1)=pos_cart(i,1)+position(i,j)*lat%direct(j,1)
! pos_cart(i,2)=pos_cart(i,2)+position(i,j)*lat%direct(j,2)
! pos_cart(i,3)=pos_cart(i,3)+position(i,j)*lat%direct(j,3)
! enddo
! enddo
! endif
! ! Read OUTCAR file beyond this point
! allocate(pos_outcar(total_atom,3),force_outcar(total_atom,3))
! allocate(pos_outcar_new(total_atom,3))
! allocate(rij(total_atom,3))
! open(outread,file="./OUTCAR")
! do i=1,maxl
! read(outread,'(A)',iostat=reason) line
! word=' POSITION'
! call findword(word,line,found,p)
! if(found) then
! read(outread,'(A)',iostat=reason) line
! do j=1,total_atom
! read(outread,*,iostat=reason) pos_outcar(j,1),pos_outcar(j,2), &
! pos_outcar(j,3), force_outcar(j,1), force_outcar(j,2), force_outcar(j,3)
! !write(*,*) pos_outcar(j,1),pos_outcar(j,2),pos_outcar(j,3), &
! !force_outcar(j,1),force_outcar(j,2),force_outcar(j,3)
! enddo
! endif
! if(reason < 0) then
! write(*,*) ".......END OF OUTCAR......."
! exit
! endif
! enddo
! ! Get the distance of each atom from POSCAR and OUTCAR file
! do i=1,total_atom
! do j=1,3
! rij(i,j)=pos_outcar(i,j)-pos_cart(i,j)
! enddo
! enddo
! ! Try to change the cartesian coordinates from OUTCAR to reduce coordinates
! write(*,*) "--------------------------------------------------------------------"
! do i=1,total_atom
! lat%a=0.0d0
! pos_outcar_new(i,1)=0.0d0
! pos_outcar_new(i,2)=0.0d0
! pos_outcar_new(i,3)=0.0d0
! do j=1,3
! lat%a(1)=lat%a(1)+pos_outcar(i,j)*lat%reciprocal(1,j)
! lat%a(2)=lat%a(2)+pos_outcar(i,j)*lat%reciprocal(2,j)
! lat%a(3)=lat%a(3)+pos_outcar(i,j)*lat%reciprocal(3,j)
! enddo
! write(*,*) lat%a(1),lat%a(2),lat%a(3)
! lat%a(1)=lat%a(1)-nint(lat%a(1))
! lat%a(2)=lat%a(2)-nint(lat%a(2))
! lat%a(3)=lat%a(3)-nint(lat%a(3))
! write(*,*) lat%a(1),lat%a(2),lat%a(3)
! !lat%a(1)=anint(lat%a(1))
! !lat%a(2)=anint(lat%a(2))
! !lat%a(3)=anint(lat%a(3))
! write(*,*) position(i,1), position(i,2), position(i,3)
! write(*,*) "-------------------------------------------------------------------"
! do j=1,3
! pos_outcar_new(i,j)=lat%a(1)*lat%direct(1,j)+lat%a(2)*lat%direct(2,j)+ &
! lat%a(3)*lat%direct(3,j)
! enddo
! !write(*,*) lat%a(1),lat%a(2),lat%a(3)
! enddo
! close(outread)
! write(*,*) "--------------------------------------------------------------------"
! write(*,*) position(200,3)
! deallocate(atom_type_num,position,pos_cart)
! deallocate(pos_outcar,force_outcar,pos_outcar_new)
end program image
