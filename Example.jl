include("Include.jl")

# set the system T -
temperature_in_K = 283.15
R = 8.314

# components -
# 1 => benzene
# 2 => cyclohexane

function compute_saturation_pressure_array(temperature::Float64 = 273.15)

    # initialize -
    T = temperature - 273.15    # T has to be in units of C
    saturation_pressure_array = zeros(2)

    # Antoine parameters (SVN Table B.1)
    # benzene
    A1 = 13.7819
    B1 = 2726.81
    C1 = 217.572

    # cyclohexane -
    A2 = 13.6568
    B2 = 2723.44
    C2 = 220.618

    # Antoine -
    P1_sat = A1 - (B1/(T+C1))
    P2_sat = A2 - (B2/(T+C2))

    # grab -
    saturation_pressure_array[1] = exp((P1_sat))
    saturation_pressure_array[2] = exp((P2_sat))

    # return -
    return saturation_pressure_array
end

function compute_model_pressure(parameter_guess_array::Array{Float64,1}, data_table::DataFrame, saturation_pressure_array::Array{Float64,1})

    # setup compostion array -
    x_array = data_table[!,:x]

    # Get the saturation pressures -
    P1_sat = saturation_pressure_array[1]   # kPa
    P2_sat = saturation_pressure_array[2]   # KPa

    # compute A-term -
    A = parameter_guess_array[1]
    A_term = (A/(R*temperature_in_K))

    # initialize some storage -
    P_array = Float64[]
    for x1_value in x_array

        # get x2 -
        x2_value = 1.0 - x1_value

        # compute gamma -
        gamma_1 = exp(A_term*(x2_value)^2)
        gamma_2 = exp(A_term*(x1_value)^2)

        # compute the pressure -
        P_value = gamma_1*x1_value*P1_sat + gamma_2*x2_value*P2_sat

        # grab -
        push!(P_array, P_value)
    end

    # return computed pressure -
    return P_array
end

function full_objective_function(parameter_guess_array::Array{Float64,1}, data_table::DataFrame, saturation_pressure_array::Array{Float64,1})::Float64

    # get the actual pressure (Pa) -
    actual_pressure = data_table[!,:P]

    # compute the model pressure -
    model_pressure = compute_model_pressure(parameter_guess_array, data_table, saturation_pressure_array)

    # compute the MSE -
    MSE = sum((model_pressure .- actual_pressure).^2)

    # return -
    return MSE
end

function solve()

    # load the load_experimental_data_file -
    data_table = load_experimental_data_file("./data/Data.csv")

    # compute single component saturation pressures -
    T = temperature_in_K
    sat_pressure_array = compute_saturation_pressure_array(T)

    # Setup some bounds and an initial guess -
    lower_bound = 0.0
    upper_bound = 4000.0
    initial_guess = [1000.0]

    # alias the objective_function -
    objective_function(x) = full_objective_function(x,data_table,sat_pressure_array)

    # call the solver -
    # make a call to the optim package -
    result = optimize(objective_function,lower_bound,upper_bound,initial_guess,Fminbox(LBFGS()))

    # get the min energy composition array -
    parameter_value = Optim.minimizer(result)

    # return -
    return parameter_value
end



function evaluate_model(parameter_array::Array{Float64,1})

    # load the load_experimental_data_file -
    data_table = load_experimental_data_file("./data/Data.csv")

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
    for x1_value in x_array

        # get x2 -
        x2_value = 1.0 - x1_value

        # compute gamma -
        gamma_1 = exp(A_term*(x2_value)^2)
        gamma_2 = exp(A_term*(x1_value)^2)

        # compute the pressure -
        P_value = gamma_1*x1_value*P1_sat + gamma_2*x2_value*P2_sat

        # grab -
        push!(P_array, P_value)
    end

    return (x_array, P_array)
end
