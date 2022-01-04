program motodens
	use cubeClass
	implicit none
	character(len=200)::input,occupationString,output
	real(kind=8)::occupation,dx,dy,dz,sumOverCompleteSpace
	type(cube)::cube1
	integer::io,lastSlash
	
	
	call getarg(1,input)
	call getarg(2,occupationString)
	call getarg(3,output)
	
	cube1 = readCube(input)
	
	dx=cube1%axis_x(1)
	dy=cube1%axis_y(2)
	dz=cube1%axis_z(3)

	if (trim(occupationString) == "")then
		write(*,*) 'The occupation of the MO hast to be given as second parameter!' 
		stop		
	end if 
	read ( occupationString, '(F50.0)', iostat=io ) occupation
	if ( io /= 0 ) THEN
		write(*,*) 'The second parameter "',trim(occupationString) , '"could not be converted to a number!' 
		stop
	end if


	if (trim(output) == "")then
		lastSlash = scan(trim(input),"/", BACK= .true.)
		output = trim(input(1:lastSlash))//"dens"//trim(occupationString)//"-"//trim(input(lastSlash+1:))
	end if
	
	sumOverCompleteSpace= sum(cube1%values) 
	write(*,*) sumOverCompleteSpace
	
	sumOverCompleteSpace= sum(cube1%values) * dx * dy *dz
	write(*,*) sumOverCompleteSpace
	
	cube1%values = cube1%values * cube1%values 
	
	sumOverCompleteSpace= sum(cube1%values) * dx * dy *dz
	write(*,*) sumOverCompleteSpace
	
	cube1%values = cube1%values* occupation
	
	!cube1%values = cube1%values * cube1%values
	!sumOverCompleteSpace= sum(cube1%values) * dx * dy *dz
	
	!cube1%values = cube1%values / sumOverCompleteSpace * occupation

	call writeCube(cube1,output)
	
end program motodens
