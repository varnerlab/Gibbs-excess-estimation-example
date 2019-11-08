function compute_Pxy_raoults_law(temperature_in_K::Float64)

    # compute single component saturation pressures -
    T = temperature_in_K
    sat_pressure_array = compute_saturation_pressure_array(T)

    # setup compostion array -
    x1_array = collect(0:0.01:1.0)

    # Get the saturation pressures -
    P1_sat = sat_pressure_array[1]   # kPa
    P2_sat = sat_pressure_array[2]   # KPa

    # main loop -
    P_array = Float64[]
    y1_array = Float64[]
    for x1_value in x1_array

        # x2 -
        x2_value = (1 - x1_value)

        # compute the pressure -
        P_value = x1_value*P1_sat + x2_value*P2_sat

        # compute y1 -
        y1_value = (1/P_value)*x1_value*P1_sat

        # grab -
        push!(P_array, P_value)
        push!(y1_array,y1_value)
    end


    # return -
    return (x1_array, y1_array, P_array)
end

function compute_Pxy_mod_raoults_law(temperature_in_K::Float64, parameter_array::Array{Float64,1})

    # compute single component saturation pressures -
    T = temperature_in_K
    sat_pressure_array = compute_saturation_pressure_array(T)

    # setup compostion array -
    x1_array = collect(0:0.01:1.0)

    # Get the saturation pressures -
    P1_sat = sat_pressure_array[1]   # kPa
    P2_sat = sat_pressure_array[2]   # KPa

    # compute single component saturation pressures -
    T = temperature_in_K
    sat_pressure_array = compute_saturation_pressure_array(T)

    # setup compostion array -
    x_array = collect(0:0.01:1.0)

    # Get the saturation pressures -
    P1_sat = sat_pressure_array[1]   # kPa
    P2_sat = sat_pressure_array[2]   # KPa

    # compute A-term -
    A = parameter_array[1]
    A_term = (A/(R*temperature_in_K))

    # initialize some storage -
    P_array = Float64[]
    y1_array = Float64[]
    for x1_value in x_array

        # get x2 -
        x2_value = 1.0 - x1_value

        # compute gamma -
        gamma_1 = exp(A_term*(x2_value)^2)
        gamma_2 = exp(A_term*(x1_value)^2)

        # compute the pressure -
        P_value = gamma_1*x1_value*P1_sat + gamma_2*x2_value*P2_sat

        # compute y1 -
        y1_value = (1/P_value)*gamma_1*x1_value*P1_sat

        # grab -
        push!(P_array, P_value)
        push!(y1_array, y1_value)
    end

    return (x_array, y1_array, P_array)
end
