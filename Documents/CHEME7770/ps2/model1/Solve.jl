include("Include.jl")

# load the data dictionary -
#data_dictionary = maximize_acetate_data_dictionary(0,0,0)
#data_dictionary = maximize_atp_data_dictionary(0,0,0)
#data_dictionary = maximize_cellmass_data_dictionary(0,0,0)

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)

epsilon = 0.1

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)

#show_flux_profile(flux_array, epsilon, data_dictionary)

path_to_atom_file = "/Users/karenna/Documents/CHEME7770/ps2/model1/Atom.txt"

#X = generate_atom_matrix(path_to_atom_file, data_dictionary)

#Y = transpose(X)

#unbalanced_array = *(Y,uptake_array)

checkAllBalances()
