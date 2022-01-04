program alignToZ
	implicit none
	integer :: i,io
	integer::atomCount
	character(len=100)::input,output,atomCountString
	character(len=200)::line
	character(len=6),allocatable,dimension(:)::elements
	real(kind=8)::n1dist,n2dist,alphaX,alphaY
	real(kind=8),allocatable,dimension(:,:)::coords
	real(kind=8),dimension(3)::metalCoords,n1Coords,n2Coords,nCoords,tmpCoords
	real(kind=8),dimension(3,3)::Rx,Ry
	
	
	call getarg(1,input)
	call getarg(2,output)
	
	if (trim(output) == "")then
		output=trim(input)
	end if 
	
	open(unit=11,file=input,status='old')
	
	read(11,"(A)") atomCountString
	
	if (trim(atomCountString) == "")then
		write(*,*) 'The input file could not be read!' 
		stop		
	end if 
	read ( atomCountString, '(I4)', iostat=io ) atomCount
	if ( io /= 0 ) THEN
		write(*,*)  'The input file could not be read! Exit!' 
		stop
	end if
	
	
	allocate(coords(atomCount,3))
	allocate(elements(atomCount))
	
	!Read comment line
	read(11,"(A)") line

	do i = 1, atomCount
			read(11,*) elements(i),coords(i,1),coords(i,2),coords(i,3)
	end do
	close(11)
	
	!The N2 ligand has to be on position 1 & 2 and the metal atom on position 3
	!There is no check for this!
	
	!Put metal in the origin.
	metalCoords(1) = coords(3,1) !x
	metalCoords(2) = coords(3,2) !y
	metalCoords(3) = coords(3,3) !z
	
	do i = 1, atomCount
			coords(i,1) = coords(i,1) - metalCoords(1)
			coords(i,2) = coords(i,2) - metalCoords(2)
			coords(i,3) = coords(i,3) - metalCoords(3)
	end do
	
	n1Coords(1) = coords(1,1) !x
	n1Coords(2) = coords(1,2) !y
	n1Coords(3) = coords(1,3) !z
	n2Coords(1) = coords(2,1) !x
	n2Coords(2) = coords(2,2) !y
	n2Coords(3) = coords(2,3) !z
	metalCoords(1) = coords(3,1) !x should be 0
	metalCoords(2) = coords(3,2) !y should be 0
	metalCoords(3) = coords(3,3) !z should be 0
	
	!Determine N atom next to metal.	
	n1dist = dsqrt ( n1Coords(1)**2 + n1Coords(2)**2 + n1Coords(3)**2 )
	n2dist = dsqrt ( n2Coords(1)**2 + n2Coords(2)**2 + n2Coords(3)**2 )
	
	if ( n1dist < n2dist ) then
		nCoords = n1Coords
	else	
		nCoords = n2Coords
	end if
	
	!Spin molecule so that the metalCoords und nCoords are on the z-axis.
	
	!Angle and matrix for x-axis
	if ( nCoords(3) == 0 .AND. nCoords(2) == 0 ) then
		alphaX = 0.0_8
	else
		alphaX = dacos( nCoords(3) / ( dsqrt( nCoords(2)**2 + nCoords(3)**2) ) )
	end if

	if ( nCoords(2) < 0 ) then
		alphaX = alphaX * (-1_8)
	end if
	Rx(1,1) = 1.0_8
	Rx(2,1) = 0.0_8
	Rx(3,1) = 0.0_8
	Rx(1,2) = 0.0_8
	Rx(2,2) = dcos(alphaX)
	Rx(3,2) = dsin(alphaX)
	Rx(1,3) = 0.0_8
	Rx(2,3) = -dsin(alphaX)
	Rx(3,3) = dcos(alphaX)
	
	do i = 1, atomCount
			tmpCoords(1) = coords(i,1)	
			tmpCoords(2) = coords(i,2)	
			tmpCoords(3) = coords(i,3)	
			
			tmpCoords = matmul(Rx,tmpCoords)
			
			coords(i,1) = tmpCoords(1)
			coords(i,2) = tmpCoords(2)
			coords(i,3) = tmpCoords(3)			
	end do
	
	n1Coords(1) = coords(1,1) !x
	n1Coords(2) = coords(1,2) !y
	n1Coords(3) = coords(1,3) !z
	n2Coords(1) = coords(2,1) !x
	n2Coords(2) = coords(2,2) !y
	n2Coords(3) = coords(2,3) !z
	if ( n1dist < n2dist ) then
		nCoords = n1Coords
	else	
		nCoords = n2Coords
	end if
	
	!Angle and matrix for y-axis
	if ( nCoords(3) == 0 .AND. nCoords(1) == 0 ) then
		alphaY = 0.0_8
	else
		alphaY = - dacos( nCoords(3) / ( dsqrt( nCoords(1)**2 + nCoords(3)**2) ) )
	end if
	if ( nCoords(1) < 0 ) then
		alphaY = alphaY * (-1_8)
	end if
	
	Ry(1,1) = dcos(alphaY)
	Ry(2,1) = 0.0_8
	Ry(3,1) = -dsin(alphaY)
	Ry(1,2) = 0.0_8
	Ry(2,2) = 1.0_8
	Ry(3,2) = 0.0_8
	Ry(1,3) = dsin(alphaY)
	Ry(2,3) = 0.0_8
	Ry(3,3) = dcos(alphaY)
	
	do i = 1, atomCount
			tmpCoords(1) = coords(i,1)	
			tmpCoords(2) = coords(i,2)	
			tmpCoords(3) = coords(i,3)	
			
			tmpCoords = matmul(Ry,tmpCoords)
			
			coords(i,1) = tmpCoords(1)
			coords(i,2) = tmpCoords(2)
			coords(i,3) = tmpCoords(3)			
	end do
	
	open(unit=13,file=output,status='replace')
	
	write(13,"(I4)") atomCount
	write(13,*) 'The coordinates in "', trim(input) ,'" were aligned to the z-axis and saved in "',trim(output),'"'
	do i = 1, atomCount
		write(13,'(A,A,F16.10,F16.10,F16.10)') elements(i),'  ',coords(i,1),coords(i,2),coords(i,3)
	end do
	write(13,"(A)") ' '
	
	close(13)
	
	write(*,*) 'The coordinates in "', trim(input) ,'" were aligned to the z-axis"'
	
	deallocate(coords)
	deallocate(elements)


end program alignToZ



