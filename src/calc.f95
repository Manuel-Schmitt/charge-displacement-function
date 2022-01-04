program calc		
	use cubeClass
	implicit none
	character(len=200)::input1,input2,output
	character(len=1)::op
	type(cube)::cube1, cube2, result_cube
	
	logical::exists,input2IsNumber=.false.
	integer::io
	real(kind=8)::value
	
	
	call getarg(1,input1)
	call getarg(2,op)
	call getarg(3,input2)
	call getarg(4,output)
	
	cube1 = readCube(input1)
		
	INQUIRE(FILE=input2, EXIST=exists) 		
	if(.not. exists) then
		if (trim(input2) == "")then
			write(*,*) 'No secend input file was specified!' 
			stop		
		end if 
		read ( input2, '(F50.0)', iostat=io ) value
		if ( io /= 0 ) THEN
			write(*,*) 'The input file "', trim(input2), '" does not exist and is not a number!' 
			stop
		end if
		input2IsNumber=.true.
	else
		cube2 = readCube(input2)
	
		if (cube1%nx.NE.cube2%nx.OR.cube1%ny.NE.cube2%ny.OR.cube1%nz.NE.cube2%nz.OR. & 
		cube1%origin_x.NE.cube2%origin_x.OR.cube1%origin_y.NE.cube2%origin_y.OR.cube1%origin_z.NE.cube2%origin_z.OR. &
		all(cube1%axis_x.NE.cube2%axis_x).OR.all(cube1%axis_y.NE.cube2%axis_y).OR.all(cube1%axis_z.NE.cube2%axis_z)) then
			write(*,*) 'The grids of the input files "',trim(input1),'" and "',trim(input2),'" is not equal!' 
			stop
		end if
	end if	
	
	!result_cube gets header, grid and atoms of cube1
	result_cube=cube1
	
	  select case( op )
		case( "+" )
			if (input2IsNumber) then
				result_cube%values = cube1%values + value
			else
				result_cube%values = cube1%values + cube2%values
			end if
		case( "-" ) 
			if (input2IsNumber) then
				result_cube%values = cube1%values - value
			else
				result_cube%values = cube1%values - cube2%values
			end if    
		case( "*" )
			if (input2IsNumber) then
				result_cube%values = cube1%values * value
			else
				result_cube%values = cube1%values * cube2%values
			end if
		case( "/" )
			if (input2IsNumber) then
				result_cube%values = cube1%values / value
			else
				result_cube%values = cube1%values / cube2%values
			end if
		case default
			write(*,*) 'The specified operator "', op ,'" was not recognized!' 
			stop
		end select 

	result_cube%header1 = "Result of "//trim(input1)//" "//op//" "//trim(input2)
	call writeCube(result_cube,output)
	
end program calc