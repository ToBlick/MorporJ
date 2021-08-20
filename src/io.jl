using DelimitedFiles

"""
load_from_txt(path)

Returns the content of the delimited file at path (generated e.g. using writedlm in Julia or np.savetxt in Python)
Checks if all entires are non-negative
"""
function load_from_txt(path)

    S = readdlm(path)

    return readdlm(path)
end
