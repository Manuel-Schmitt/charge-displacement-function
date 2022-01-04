program integrateCube
	use cubeClass
	implicit none
	character(len=200)::input,output,integralOutput,cdfOutput
	type(cube)::cube1

	real(kind=8)::z,cdf,dx,dy,dz,bohrToAngstrom,integralToAngstrom
	integer::i,lastSlash,lastDot
	real(kind=8),allocatable,dimension(:)::integral
	
	
	call getarg(1,input)
	call getarg(2,output)
	
	cube1 = readCube(input)
	
	if (trim(output) == "")then
		lastSlash = scan(trim(input),"/", BACK= .true.)
		lastDot = scan(trim(input),".", BACK= .true.)
		if (lastDot == 0)then
			lastDot=len(input)
		else
			lastDot=lastDot-1
		end if
		integralOutput = trim(input(1:lastSlash))//'integral-'//trim(input(lastSlash+1:lastDot))//".dat"
		cdfOutput = trim(input(1:lastSlash))//"cdf-"//trim(input(lastSlash+1:lastDot))//".dat"
	else	
		integralOutput="integral-"//trim(output)
		cdfOutput="cdf-"//trim(output)
	end if

	dx=cube1%axis_x(1)
	dy=cube1%axis_y(2)
	dz=cube1%axis_z(3)
	z=cube1%origin_z
	allocate(integral(cube1%nz))
	
	integral=sum(sum(cube1%values,1),1) * dx * dy


	open(unit=11,file=integralOutput,status='replace')
	open(unit=12,file=cdfOutput,status='replace')
	
	cdf = 0_8
	do i = 1, cube1%nz
		
		write(11,*) bohrToAngstrom(z),integralToAngstrom(integral(i))
		
		cdf = cdf + (integral(i) * dz)
		write(12,*) bohrToAngstrom(z),cdf
		
		z = z + dz
	end do
	
	close(11)
	close(12)
	
	deallocate(integral)
	
end program integrateCube

real(kind=8) function bohrToAngstrom(a)
	implicit none
	real(kind=8),intent(in)::a
	bohrToAngstrom = a * 0.529177210903_8		
end function bohrToAngstrom

real(kind=8) function integralToAngstrom(a)
	implicit none
	real(kind=8),intent(in)::a
	integralToAngstrom = a / 0.529177210903_8		
end function integralToAngstrom

