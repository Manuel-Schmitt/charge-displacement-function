module cubeClass
	implicit none
	save
	
	type::cube
		character(len=200)::header1,header2
		integer::natoms,nx,ny,nz
		real(kind=8)::origin_x,origin_y,origin_z
		real(kind=8),dimension(3)::axis_x,axis_y,axis_z
		integer,allocatable,dimension(:)::elements
		real(kind=8),allocatable,dimension(:)::charges
		real(kind=8),allocatable,dimension(:,:)::coords
		real(kind=8),allocatable,dimension(:,:,:)::values
	end type cube
  
  
	contains
		type(cube) function readCube(path)
		
			character(len=200),intent(in)::path
			integer:: i1,i2,i3		
			logical :: exists
			
			INQUIRE(FILE=path, EXIST=exists) 		
			if(.not. exists) then
				write(*,*) 'The input file "', trim(path), '" does not exist!' 
				stop
			end if
			
			open(unit=11,file=path,status='old')
	
			read(11,"(A)") readCube%header1
			read(11,"(A)") readCube%header2
			
			
			read(11,*) readCube%natoms,readCube%origin_x,readCube%origin_y,readCube%origin_z

			if (readCube%natoms <0 ) then
				write(*,*) 'The input file "', trim(path), '" contains MOs, which is not supported!' 
				stop
			end if 
			
			read(11,*) readCube%nx, readCube%axis_x
			read(11,*) readCube%ny, readCube%axis_y
			read(11,*) readCube%nz, readCube%axis_z
			
			allocate(readCube%elements(readCube%natoms))
			allocate(readCube%charges(readCube%natoms))
			allocate(readCube%coords(readCube%natoms,3))
			allocate(readCube%values(readCube%nx,readCube%ny,readCube%nz))
			
			do i1 = 1, readCube%natoms
				read(11,*) readCube%elements(i1),readCube%charges(i1),readCube%coords(i1,:)
			end do
			
			do i1 = 1, readCube%nx
				do i2 = 1, readCube%ny
					read(11,*) (readCube%values(i1,i2,i3),i3=1,readCube%nz)
				end do
			end do
					
			close(11)
			
		end function readCube
		
		subroutine writeCube(output,path)
			type(cube), intent (in)  :: output 
			character(len=200),intent(in)::path
			
			integer:: i1,i2,i3		
			open(unit=12,file=path,status='REPLACE')
			
			write(12,"(A)") trim(output%header1)
			write(12,"(A)") trim(output%header2)
			write(12,"(I5,3F12.6)") output%natoms,output%origin_x,output%origin_y,output%origin_z 

			write(12,"(I5,3F12.6)") output%nx, output%axis_x(:)
			write(12,"(I5,3F12.6)") output%ny, output%axis_y(:)
			write(12,"(I5,3F12.6)") output%nz, output%axis_z(:)

			do i1 = 1, output%natoms
				write(12,"(I5,4F12.6)") output%elements(i1),output%charges(i1),output%coords(i1,:)
			end do
			
			do i1 = 1, output%nx
				do i2 = 1, output%ny
					write(12,"(6E18.10)") (output%values(i1,i2,i3),i3=1,output%nz)
				end do
			end do

			close(12)
		end subroutine
		
end module cubeClass
