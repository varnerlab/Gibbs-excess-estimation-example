function load_experimental_data_file(path_to_data_file::String)::DataFrame

    # use CSV to load data -
    dt = CSV.read(path_to_data_file)

    # return -
    return dt
end

# function compute_saturation_pressure(parameter_object::AntoineParameters, temperature_in_K::Float64)::Float64
#
#     # Get the parameters -
#     A = parameter_object.A
#     B = parameter_object.B
#     C = parameter_object.C
#
#     # compute -
#     T = temperature_in_K - 273.15 # need to convert to C, varies by source of parameters -
#     ln_PSat = A - B/(T+C)
#
#     # return -
#     return exp(ln_PSat)
# end
