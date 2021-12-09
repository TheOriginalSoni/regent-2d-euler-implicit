import "regent"
local c = terralib.includec("stdio.h")

fspace Point
{
    x : double,
    y : double,
    nx : double,
    ny : double,
    left : int,
    right : int,
    status : int,
    geo_index : int,
    min_dist : double,
    vor_area : double,
    nbhs : int,
    conn : int[25],
    xpos_nbhs : int,
    xneg_nbhs : int,
    ypos_nbhs : int,
    yneg_nbhs : int,
    prim : double[4],
    flux_res : double[4],
    q : double[4],
    qm : double[4][2],
    delt : double,
    entropy : double,
    alias : int,
    point_with_alias : int,
    U : double[4],
    delUp : double[4],
    delUn : double[4],
    D_inv : double
}

fspace pointvars
{
    max_points : int,
    wall_points : int,
    outer_points : int,
    interior_points : int,
    shape_points : int,
    no_of_shapes : int,

    res_old : double,
    res_new : double,
    residue : double,
    max_res : double,
    max_res_point : int,
    cfv : double,
    Cl : double,
    Cd : double,
    Cm : double,
    total_entropy : double,
    entropy : double
}

fspace pointvarsarrays
{
    wall_points_index : int,
    outer_points_index : int,
    interior_points_index : int,
    shape_points_index : int
}

task init()
    var max_points = 100
    var points = region ( ispace (int1d , max_points ) , Point )
    return points
end

task input(max_points : int, points : region(ispace(int1d), Point))
where reads writes(points) do
    var input_file = "testcases/naca0012-38400/2d_input_data"
    var file = c.fopen(input_file, "r")
    c.fscanf(file, "%i", max_points)

    for i = 1, max_points do
        c.fscanf(file,"%i", points[i].x)
        c.fscanf(file,"%i", points[i].y)
        c.fscanf(file,"%i", points[i].left)
        c.fscanf(file,"%i", points[i].right)
        c.fscanf(file,"%i", points[i].status)
        c.fscanf(file,"%i", points[i].geo_index)
        c.fscanf(file,"%i", points[i].min_dist)
        c.fscanf(file,"%i", points[i].vor_area)
        c.fscanf(file,"%i", points[i].nbhs)
        for j = 1,points[i].nbhs do
            c.fscanf(file,"%i", points[i].conn[j])
        end
    end

    c.printf("%i\n",max_points)
    c.fclose(file)

    var no_of_shapes = 0
    var shape_points = 0
    var interior_points = 0
    var wall_points = 0
    var outer_points = 0

    for i = 1,max_points do
        if(points[i].geo_index > no_of_shapes) then
            no_of_shapes = points[i].geo_index
        end    
        
        if(points[i].geo_index > 0) then
            shape_points = shape_points + 1
        end

        if(points[i].status==1) then
            interior_points = interior_points + 1
        end
        if(points[i].status==0) then
            wall_points = wall_points + 1
        end
        if(points[i].status==2) then
            outer_points = outer_points + 1
        end
    end

    var shape_points_index = region(ispace(int1d , shape_points),int)
    var interior_points_index = region(ispace(int1d , interior_points),int)
    var wall_points_index = region(ispace(int1d , wall_points),int)
    var outer_points_index = region(ispace(int1d , outer_points),int)
    
    var shape_count = 0
    var interior_count = 0
    var wall_count = 0
    var outer_count = 0

    for i= 1,max_points do
        if(points[i].geo_index > 0) then
            shape_count = shape_count + 1
            shape_points_index[shape_count] = i
        end

        if(points[i].status==1) then 
            interior_count = interior_count + 1
            interior_points_index[interior_count] = i
        end
        
        if(points[i].status==0) then
            wall_count = wall_count + 1
            wall_points_index[wall_count] = i
        end
        
        if(points[i].status==2) then
            outer_count = outer_count + 1
            outer_points_index[outer_count] = i
        end
    end
end

task parameter_mod()
    var max_iters=30000
    var Mach=1.2
    var aoa = 0.0

    val rho_inf = 1.0
    val pr_inf = 1.0/1.4
    val gamma = 1.4
    val pi=4.0d0*datan(1.0d0)
    val theta = aoa*pi/180.0

    val u1_inf = Mach*dcos(theta)
    val u2_inf = Mach*dsin(theta)

    --CFL is the CFL number for stability ..
    val CFL = 1.0

    --The parameter power is used to specify the weights in the LS formula for the derivatives power = 0.0d0, -2.0d0, -4.0d0, -6.0d0 ..
    --For example, power = -2.0 implies that power = -2.0 => weights = 1/d^2
    --power = -4.0 => weights = 1/d^4
    val power = -0.0d0 

    --The following parameter decides the spatial order of accuracy. If second_order_flag = 0 then we get first-order solution. If second_order_flag = 1, then it leads to second-order scheme.
    val second_order_flag = 1

    --Venkatakrishnan limiter constant ..
    val VL_CONST = 10.0d0

    --Used to get accurate q-derivatives ..
    val inner_iterations = 3

    --flag for initial conditions. 
    --if flag = 0 => free stream conditions ..
    --if flag = 1 => read initial conditions from a restart file ..
    val initial_conditions_flag = 0

    --flag for the implcit scheme based on LUSGS ..
    --If flag = 1 => implicit solver based on LUSUS.
    val implicit_scheme_flag = 1

    --Parameters for the adjoint code ..

    --If Cl_flag = 1.0, then it is the objective function. 
    --Otherwise if Cl_flag = 0.0, then it is not the objective function 
    val Cl_flag = 0.0
    val Cd_flag = 0.0
    val Cm_flag = 0.0
    val Cl_Cd_flag = 0.0
    val entropy_flag = 0.0
    val enstrophy_flag = 0.0
end

task hello_world()
  c.printf("Hello World!\n")
end

task main()
  hello_world()
  var points = init()
  var max_points = 100
  c.printf("%i\n",max_points)
  input(max_points,points)
  c.printf("%i\n",max_points)
end

regentlib.start(main)


