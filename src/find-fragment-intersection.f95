program findFragmentIntersection
	use cubeClass
	implicit none
	type(cube)::cube1,cube2
	character(len=200)::input1,input2,output
	
	real(kind=8)::x,y,z,dx,dy,dz,zFragmentIntersection,densToAngstrom,bohrToAngstrom
	real(kind=8)::x0Value,y0Value,zTmp,tmpval1,tmpval2,deltaZ,zDensity1AtIntersection,zDensity2AtIntersection
	integer::x0,y0,z0Length,z0SP
	real(kind=8),allocatable,dimension(:)::zDensity1,zDensity2,difDensity,z0Value,z0Density
	integer(kind=8),allocatable,dimension(:)::z0
	
	
	integer :: i1,i2,i3
	
	
	call getarg(1,input1)
	call getarg(2,input2)
	call getarg(3,output)
	
	if (trim(output) == "")then
		output="intersection.dat"
	end if
	
	cube1 = readCube(input1)
	cube2 = readCube(input2)
	
	if (cube1%nx.NE.cube2%nx.OR.cube1%ny.NE.cube2%ny.OR.cube1%nz.NE.cube2%nz.OR. & 
		cube1%origin_x.NE.cube2%origin_x.OR.cube1%origin_y.NE.cube2%origin_y.OR.cube1%origin_z.NE.cube2%origin_z.OR. &
		all(cube1%axis_x.NE.cube2%axis_x).OR.all(cube1%axis_y.NE.cube2%axis_y).OR.all(cube1%axis_z.NE.cube2%axis_z)) then
			write(*,*) 'The grids of the input files "',trim(input1),'" and "',trim(input2),'" is not equal!' 
			stop
	end if
	
	x=cube1%origin_x
	y=cube1%origin_y
	z=cube1%origin_z
	dx=cube1%axis_x(1)
	dy=cube1%axis_y(2)
	dz=cube1%axis_z(3)
	
	allocate(zDensity1(cube1%nz))
	allocate(zDensity2(cube1%nz))
	allocate(difDensity(cube1%nz))
	allocate(z0(cube1%nz))
	allocate(z0Value(cube1%nz))
	allocate(z0Density(cube1%nz))
	
	!!!!!!!!!!!
	!It cannot be assumed that a data point lies exactly on the z-axis.
	!Therefore, the data points (x0,y0), (x0-1,y0), (x0,y0-1) and (x0-1,y0-1) around the z-axis are determined.
	!x0 and y0 are the array indexes.
	!x0Value and y0Value the x and y coordinates of point (x0,y0).
	!!!!!!!!!!!
	
	x0Value = x
	do i1 = 1, cube1%nx
		if ((x0Value<0) .NEQV. (x<0))then
			x0 = I1
			exit
		end if
		x0Value = x0Value + dx
	end do
	
	y0Value = y
	do i2 = 1, cube1%ny
		if ((y0Value<0) .NEQV. (y<0))then
			y0 = i2
			exit
		end if
		y0Value = y0Value + dy
	end do
	
	!Interpolation of the electron density on the z-axis for each fragment
	!First along the x-axis  -> tmpval1 and tmpval2
	!Then along the y-axis -> zDensity

	
	do i3 = 1, cube1%nz
		tmpval1 = cube1%values(x0,y0,i3) - ( x0Value * (cube1%values(x0,y0,i3) - cube1%values(x0-1,y0,i3)) / dx)	
		tmpval2 = cube1%values(x0,y0-1,i3) - ( x0Value * (cube1%values(x0,y0-1,i3) - cube1%values(x0-1,y0-1,i3)) / dx)
		zDensity1(i3) = tmpval1 - ( y0Value * (tmpval1 - tmpval2) / dy)
		
		tmpval1 = cube2%values(x0,y0,i3) - ( x0Value * (cube2%values(x0,y0,i3) - cube2%values(x0-1,y0,i3)) / dx)	
		tmpval2 = cube2%values(x0,y0-1,i3) - ( x0Value * (cube2%values(x0,y0-1,i3) - cube2%values(x0-1,y0-1,i3)) / dx)
		zDensity2(i3) = tmpval1 - ( y0Value * (tmpval1 - tmpval2) / dy)
	end do
	
	!Calculate difference electron density between both fragments on the z-axis 
	difDensity = zDensity1 - zDensity2
	
	
	!Determine all sign changes of difDensity and save them in z0 (array indexes), z0Values (z coordinates)
    !and z0Density (Sum of electron density of both fragment at this point) 	
	zTmp = z + dz
	z0Length = 0
	do i3 = 2, cube1%nz	
		if ((difDensity(i3-1) < 0) .NEQV. (difDensity(i3)<0))then
			z0Length = z0Length + 1
			z0(z0Length) = i3
			z0Value(z0Length) = zTmp
			z0Density(z0Length) = zDensity1(i3) + zDensity2(i3)
		end if 
		zTmp = zTmp + dz
	end do
	
	!The point on the z-axis, where the electron density of both fragments is the same, 
	!has the highes value of z0Density of all sign chane points
	if (z0Length == 0) then
		write(*,*) 'No fragment intersection found! Exit!'
		stop
	end if	
	
	z0SP = maxloc(z0Density(1:z0Length),1)
	
	!write(*,*) 'Found sign changes:'
	!write(*,*) '#zValue = ',z0Value(1:z0Length)
	!write(*,*) '#zDensity = ',z0Density(1:z0Length)
	!write(*,*) '#Index with highest density = ',z0SP
	
	!The real point is probably between the data points z0(z0SP) and z0(zoSP)-1 -> Interpolation 
	
	deltaZ = - ( difDensity(z0(z0SP))/ ((difDensity(z0(z0SP)) - difDensity(z0(z0SP)-1)) / dz))
	zFragmentIntersection = z0Value(z0SP) + deltaZ
	
	zDensity1AtIntersection= zDensity1(z0(z0SP)) + ( deltaZ * (zDensity1(z0(z0SP)) - zDensity1(z0(z0SP)-1)) / dz)
	zDensity2AtIntersection= zDensity2(z0(z0SP)) + ( deltaZ * (zDensity2(z0(z0SP)) - zDensity2(z0(z0SP)-1)) / dz)
	
	!write results
	
	open(unit=13,file=output,status='replace')
	write(13,*) '#Intersection of both fragments/ point of same e-density on the z-axis):'
	write(13,*) '#z[Angstrom] = ',bohrToAngstrom(zFragmentIntersection)
	write(13,*) '#e-Density(Frag1)[e Angstrom^-3] = ',densToAngstrom(zDensity1AtIntersection)
	write(13,*) '#e-Density(Frag2)[e Angstrom^-3] = ',densToAngstrom(zDensity2AtIntersection)
	write(13,*)	'#Values on the z-axis: z / e-Density(Frag1) / e-Density(Frag2) / e-Density difference'
	
	zTmp = z
	do i3 = 1, cube1%nz	
	 write(13,*) bohrToAngstrom(zTmp),densToAngstrom(zDensity1(i3)),densToAngstrom(zDensity2(i3)),densToAngstrom(difDensity(i3))
	 zTmp=zTmp+dz
	end do
	
	close(13)
	
	
end program findFragmentIntersection

real(kind=8) function densToAngstrom(a)
	implicit none
	real(kind=8),intent(in)::a
	densToAngstrom = a / (0.529177210903_8 ** 3) 		
end function densToAngstrom

real(kind=8) function bohrToAngstrom(a)
	implicit none
	real(kind=8),intent(in)::a
	bohrToAngstrom = a * 0.529177210903_8		
end function bohrToAngstrom

