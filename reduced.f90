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
! user defined constructor for lattice derived type
interface lattice
  module procedure get_structure
end interface lattice
type vaspForce
  integer(4) natom
real(8), allocatable :: position(:,:), force(:,:)
end type vaspForce
! user defined constructor for vaspForce derived type
interface vaspForce
module procedure get_vaspstructure
end interface vaspForce
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
    integer(4) posread, i, j, k, p, n, space_num, atom_type, reason, slword, lword
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
      slword=len_trim(line)-1
      lword=len_trim(line)
      do i=1,len_trim(line)-1
      if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
      atom_type=atom_type+1
      endif
      enddo
      if(iachar(line(slword:slword)) .ne. space_num .and. iachar(line(lword:lword)) .ne. space_num) then
        atom_type=atom_type+1
      endif
      if(iachar(line(slword:slword)) .eq. space_num .and. iachar(line(lword:lword)) .ne. space_num) then
        atom_type=atom_type+1
      endif
      allocate(self%atoms_type(atom_type))
      allocate(self%atoms_num(atom_type))
      k=0
      do i=1,len_trim(line)-1
      if(iachar(line(i:i)) .ne. space_num .and. iachar(line(i+1:i+1)) .ne. space_num) then
      atom_name=line(i:i+1)
      k=k+1
      self%atoms_type(k)=atom_name
      endif
      enddo
      if(iachar(line(slword:slword)) .ne. space_num .and. iachar(line(lword:lword)) .ne. space_num) then
        atom_name=line(len_trim(line)-1:len_trim(line))
        k=k+1
        self%atoms_type(k)=atom_name
      endif
      if(iachar(line(slword:slword)) .eq. space_num .and. iachar(line(lword:lword)) .ne. space_num) then
        atom_name=line(len_trim(line):len_trim(line))
        k=k+1
        self%atoms_type(k)=atom_name
      endif
      write(*,*) "The total type of atom is: ",atom_type
      write(*,*) "The total atom is: ",self%atoms_type
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
! Routine to get position and force(in any case) from VASP OUTCAR
  function get_vaspstructure(filename) result(self)
    type(vaspForce) :: self
    character(8) filename*6, line*199, word*29
    integer(4) i,j,maxl,p,outread,reason,natom
    logical found
    outread=456
    maxl=9999999
    open(outread,file=filename)
    do i=1,maxl
        read(outread,'(A)',iostat=reason) line
        word='NIONS'
        call findword(word,line,found,p)
        if ( found ) then
          read(line(p+12:),*) natom
        endif
        word=' POSITION'
        call findword(word,line,found,p)
        if (found) then
          allocate(self%position(natom,3))
          allocate(self%force(natom,3))
          read(outread,'(A)',iostat=reason) line
          do j=1,natom
            read(outread,*,iostat=reason) self%position(j,:),self%force(j,:)
          enddo
          exit
        endif
    enddo
    self%natom=natom
    close(outread)
  end function get_vaspstructure
end module poscar

program image
use poscar
type(lattice) lat
type(vaspForce) vasp
integer(4) num_args, p, i, pf
character(8) word*6
character(8), dimension(:), allocatable :: args
real(8) da,db,dc
real(8),allocatable :: rij(:,:), a(:),position(:,:), alat(:)
logical found,foundP,foundPF,foundO,foundOF
pf=789
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
      call findword(word,args(i),foundP,p)
      if (.not. foundP) then
        write(*,*) "Provide POSCAR first as an argument"
      endif
      if (foundP) then
        inquire(file=args(i),exist=foundPF)
        if (.not. foundPF) then
          write(*,*) "POSCAR NOT FOUND IN THE CURRENT DIRECTORY"
        endif
        if (foundPF) then
          write(*,*) "POSCAR FOUND IN THE CURRENT DIRECTORY"
          lat=get_structure(args(i))
        endif
      endif
    endif
    !Check if OUTCAR is present and read it
    if(i .eq. 2) then
      word='OUTCAR'
      call findword(word,args(i),foundO,p)
      if (.not. foundO) then
        write(*,*) "Provide OUTCAR as second argument"
      endif
      if (foundO) then
        inquire(file=args(i),exist=foundOF)
        if (.not. foundOF) then
          write(*,*) "OUTCAR NOT FOUND IN THE CURRENT DIRECTORY"
        endif
        if (foundOF) then
          write(*,*) "OUTCAR FOUND IN THE CURRENT DIRECTORY"
          vasp=get_vaspstructure(args(i))
        endif
      endif
    endif
    write(*,*) "FILE PROVIDED ARE: ",args(i)
  enddo
  deallocate(args)
endif

if (num_args .eq. 2) then
  if(foundP .and. foundPF) then
    allocate(rij(lat%total_atom,3))
    allocate(a(3),alat(3))
    allocate(position(lat%total_atom,3))
    write(*,*) "****************************************************************************"
    write(*,*) "****************************************************************************"
    write(*,*) "The total number of atom is: ",lat%total_atom
    do i=1,lat%total_atom
      write(*,*) lat%position(i,1), lat%position(i,2),lat%position(i,3)
      write(*,*) vasp%position(i,1),vasp%position(i,2),vasp%position(i,3)
      rij(i,1)=lat%position(i,1)-vasp%position(i,1)
      rij(i,2)=lat%position(i,2)-vasp%position(i,2)
      rij(i,3)=lat%position(i,3)-vasp%position(i,3)
      alat(1)=dot_product(lat%position(i,:),lat%reciprocal(1,:))
      alat(2)=dot_product(lat%position(i,:),lat%reciprocal(2,:))
      alat(3)=dot_product(lat%position(i,:),lat%reciprocal(3,:))
      a(1)=dot_product(vasp%position(i,:),lat%reciprocal(1,:))
      a(2)=dot_product(vasp%position(i,:),lat%reciprocal(2,:))
      a(3)=dot_product(vasp%position(i,:),lat%reciprocal(3,:))
      da=alat(1)-a(1)
      db=alat(2)-a(2)
      dc=alat(3)-a(3)
      ! This is for checking x-position in lattice
      !------------------------------******-------------------------------------
      if(da .lt. -0.8 .and. da .gt. -1) then
        a(1)=1-a(1)
      endif
      if(da .lt. 1 .and. da .gt. 0.8) then
        a(1)=a(1)-1
      endif
      if(da .lt. -0.5 .and. da .gt. -0.8) then
        a(1)=a(1)-0.5
      endif
      if(da .lt. 0.8 .and. da .gt. 0.5) then
        a(1)=a(1)+0.5
      endif
      if(da .lt. -0.1 .and. da .gt. -0.5) then
        a(1)=a(1)-0.5
        if(a(1) .lt. 0) then
          a(1)=-a(1)
        endif
      endif
      if(da .lt. 0.5 .and. da .gt. 0.1) then
        a(1)=a(1)+0.5
      endif
      ! This is for checking y-position in lattice
      !------------------------------******-------------------------------------
      if(db .lt. -0.8 .and. db .gt. -1) then
        a(2)=1-a(2)
      endif
      if(db .lt. 1 .and. db .gt. 0.8) then
        a(2)=a(2)-1
      endif
      if(db .lt. -0.5 .and. db .gt. -0.8) then
        a(2)=a(2)-0.5
      endif
      if(db .lt. 0.8 .and. db .gt. 0.5) then
        a(2)=a(2)+0.5
      endif
      if(db .lt. -0.1 .and. db .gt. -0.5) then
        a(2)=a(2)-0.5
        if(a(2) .lt. 0) then
          a(2)=-a(2)
        endif
      endif
      if(db .lt. 0.5 .and. db .gt. 0.1) then
        a(2)=a(2)+0.5
      endif
      !------------------------------******-------------------------------------
      ! This is for checking z-position in lattice
      if(dc .lt. -0.8 .and. dc .gt. -1) then
        a(3)=1-a(3)
      endif
      if(dc .lt. 1 .and. dc .gt. 0.8) then
        a(3)=a(3)-1
      endif
      if(dc .lt. -0.5 .and. dc .gt. -0.8) then
        a(3)=a(3)-0.5
      endif
      if(dc .lt. 0.8 .and. dc .gt. 0.5) then
        a(3)=a(3)+0.5
      endif
      if(dc .lt. -0.1 .and. dc .gt. -0.5) then
        a(3)=a(3)-0.5
        if(a(3) .lt. 0) then
          a(3)=-a(3)
        endif
      endif
      if(dc .lt. 0.5 .and. dc .gt. 0.1) then
        a(3)=a(3)+0.5
      endif
      !------------------------------******-------------------------------------
      position(i,1)=a(1)*lat%direct(1,1)+a(2)*lat%direct(2,1)+a(3)*lat%direct(3,1)
      position(i,2)=a(1)*lat%direct(1,2)+a(2)*lat%direct(2,2)+a(3)*lat%direct(3,2)
      position(i,3)=a(1)*lat%direct(1,3)+a(2)*lat%direct(2,3)+a(3)*lat%direct(3,3)
      write(*,*) position(i,1), position(i,2), position(i,3)
      write(*,*) "------------------------------------------------------------------------"
    enddo
  endif
endif
write(*,*) "****************************************************************************"
write(*,*) "****************************************************************************"
if (num_args .eq. 2) then
  if(foundP .and. foundPF) then
    open(pf,file='pos-force.txt')
    do i=1,lat%total_atom
      write(pf,*) position(i,1), position(i,2), position(i,3), vasp%force(i,1), vasp%force(i,2), vasp%force(i,3)
    enddo
    close(pf)
    deallocate(lat%atoms_type)
    deallocate(lat%atoms_num)
    deallocate(lat%position)
    deallocate(vasp%position)
    deallocate(vasp%force)
    deallocate(rij)
    deallocate(a,alat,position)
  endif
endif
end program image
