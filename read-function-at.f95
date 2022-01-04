program readCdfAt
	implicit none
	integer::i,io,nz,z0,numberOfDataLines
	character(len=100)::input,output,zValueString
	character(len=200)::dataString
	real(kind=8)::zValue,functionValue,z,zStart,slope,deltaZ
	real(kind=8),allocatable,dimension(:)::zData,functionData
	
	
	
	
	call getarg(1,input)
	call getarg(2,zValueString)
	call getarg(3,output)
	
	if (trim(zValueString) == "" )then
		write(*,*) 'No value was given as second parameter' 
		stop		
	end if 
	read ( zValueString, '(F100.0)', iostat=io ) zValue
	if ( io /= 0 ) then
		write(*,*) 'The second parameter "', trim(zValueString), '" is not a number!' 
		stop
	end if
	
	if (trim(output) == "" )then
		output="CDF.dat"
	end if
	
	nz = numberOfDataLines(input)
	allocate(zData(nz))
	allocate(functionData(nz))

	open(unit=11,file=input,status='old')
	
	do i = 1, nz
		read(11,'(A)') dataString
		dataString=trim(dataString) 
		if ( trim(dataString) == "") then
			cycle
		end if		
		if ( dataString(1:1) == "#") then
			cycle
		end if		
		
		read(dataString,*) zData(i),functionData(i)
	end do
	close(11)
	
	!Determine the data point z0 after zValue, resulting in:  zdate(z0-1)<zValue<zData(z0)
	
	z = zData(1) - zValue
	zStart = z 	
	do i = 1, nz
		if ((zStart<0) .NEQV. (z<0))then
			z0 = i
			exit
		end if
		z = zData(i+1) - zValue
	end do

	!Interpolate functionValue at zValue
	slope = (functionData(z0) - functionData(z0-1)) / (zData(z0) - zData(z0-1))
	deltaZ = zValue - zData(z0)
	
	functionValue = functionData(z0) + slope * deltaZ

	!write(*,*) z0
	!write(*,*) slope
	!write(*,*) deltaZ
	!write(*,*) zData(z0)
	!write(*,*) functionData(z0)

	open(unit=13,file=output,action='write',position='append')
	
	write(13,*) trim(input)//',', zValue,',',functionValue
	
	close(13)
	
	deallocate(zData)
	deallocate(functionData)

	
end program readCdfAt



integer function numberOfDataLines(path)
	implicit none
	character(len=100),intent(in)::path
	character(len=200)::line
	numberOfDataLines = 0
	open(11,file=path,status='old')
	do 
		read (11,*, end=10) line
		line=trim(line)
		if ( trim(line) == "") then
			cycle
		end if		
		if ( line(1:1) == "#") then
			cycle
		end if		
		
		numberOfDataLines = numberOfDataLines + 1 
	end do 
	10 close (11) 			
end function numberOfDataLines