# ----------------------------------------------------------------------------------- #
# Copyright (c) 2016 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2016-09-18T10:39:47
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar) 
# time_stop::Float64 => Simulation stop time value (scalar) 
# time_step::Float64 => Simulation time step (scalar) 
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs 
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)

	# stoichiometric_matrix and dilution_matrix - 
	stoichiometric_matrix = readdlm("./Network.dat")
	dilution_matrix = readdlm("./Dilution.dat")
	degradation_matrix = readdlm("./Degradation.dat")

	# array of gene lengths - 
	gene_coding_length_array = [
		15000	;	# 1	gene_1
		15000	;	# 2	gene_2
		15000	;	# 3	gene_3
	]

	# array of mRNA coding lengths - 
	mRNA_coding_length_array = [
		gene_coding_length_array[1]	;	# 4	1	mRNA_gene_1
		gene_coding_length_array[2]	;	# 5	2	mRNA_gene_2
		gene_coding_length_array[3]	;	# 6	3	mRNA_gene_3
	]

	# array of mRNA coding lengths - 
	protein_coding_length_array = [
		round((0.33)*mRNA_coding_length_array[1])	;	# 7	1	protein_gene_1
		round((0.33)*mRNA_coding_length_array[2])	;	# 8	2	protein_gene_2
		round((0.33)*mRNA_coding_length_array[3])	;	# 9	3	protein_gene_3
	]

	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 12                  # mum
	number_of_rnapII = 75000            # copies/cells
	number_of_ribosome = 1e6            # copies/cells
	mRNA_half_life_TF = 2               # hrs
	protein_half_life = 10              # hrs
	doubling_time_cell = 19.5           # hrs
	max_translation_rate = 5            # aa/sec
	max_transcription_rate = 6.0        # nt/sec
	average_transcript_length = 15000   # nt
	average_protein_length = 5000       # aa
	fraction_nucleus = 0.49             # dimensionless
	av_number = 6.02e23                 # number/mol
	avg_gene_number = 2                 # number of copies of a gene
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)
	
	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9               # nM
	
	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                    # hr^-1
	
	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1
	kcat_translation = max_translation_rate*(3600/average_protein_length)             # hr^-1
	
	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1
	
	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM
	
	# How fast do my cells die?
	death_rate_constant = 0.2*maximum_specific_growth_rate                            # hr^-1
	
	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                           # nM
	saturation_translation = 100000*(1/av_number)*(1/V)*1e9                           # nM
	# -------------------------------------------------------------------------------------------#

	# initial condition array - 
	initial_condition_array = [
		avg_gene_concentration	;	# 1	gene_1
		avg_gene_concentration	;	# 2	gene_2
		avg_gene_concentration	;	# 3	gene_3
		0.0	;	# 4	mRNA_gene_1
		0.0	;	# 5	mRNA_gene_2
		0.0	;	# 6	mRNA_gene_3
		0.0	;	# 7	protein_gene_1
		0.0	;	# 8	protein_gene_2
		0.0	;	# 9	protein_gene_3
	]

	binding_parameter_dictionary = Dict{AbstractString,Float64}()
	binding_parameter_dictionary["n_gene_1_gene_3_gene_2"] = 1.0
	binding_parameter_dictionary["K_gene_1_gene_3_gene_2"] = 10.0
	binding_parameter_dictionary["n_gene_1_gene_2"] = 1.0
	binding_parameter_dictionary["K_gene_1_gene_2"] = 10.0

	# Alias the control function parameters - 
	control_parameter_dictionary = Dict{AbstractString,Float64}()
	control_parameter_dictionary["W_gene_1_RNAP"] = 0.01
	control_parameter_dictionary["W_gene_1_gene_3_gene_2"] = 1.0
	control_parameter_dictionary["W_gene_1_gene_2"] = 1.0
	control_parameter_dictionary["W_gene_2_RNAP"] = 0.01
	control_parameter_dictionary["W_gene_3_RNAP"] = 0.01

	# Parameter name index array - 
	parameter_name_mapping_array = [
		"n_gene_1_gene_3_gene_2"	;	# 1
		"K_gene_1_gene_3_gene_2"	;	# 2
		"n_gene_1_gene_2"	;	# 3
		"K_gene_1_gene_2"	;	# 4
		"W_gene_1_RNAP"	;	# 5
		"W_gene_1_gene_3_gene_2"	;	# 6
		"W_gene_1_gene_2"	;	# 7
		"W_gene_2_RNAP"	;	# 8
		"W_gene_3_RNAP"	;	# 9
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["initial_condition_array"] = initial_condition_array
	data_dictionary["average_transcript_length"] = average_transcript_length
	data_dictionary["average_protein_length"] = average_protein_length
	data_dictionary["gene_coding_length_array"] = gene_coding_length_array
	data_dictionary["mRNA_coding_length_array"] = mRNA_coding_length_array
	data_dictionary["protein_coding_length_array"] = protein_coding_length_array
	data_dictionary["rnapII_concentration"] = rnapII_concentration  # muM 
	data_dictionary["ribosome_concentration"] = ribosome_concentration # muM 
	data_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1 
	data_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1 
	data_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1 
	data_dictionary["kcat_translation"] = kcat_translation  # hr^-1 
	data_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1 
	data_dictionary["death_rate_constant"] = death_rate_constant 
	data_dictionary["avg_gene_concentration"] = avg_gene_concentration 
	data_dictionary["saturation_constant_transcription"] = saturation_transcription 
	data_dictionary["saturation_constant_translation"] = saturation_translation 

	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["dilution_matrix"] = dilution_matrix
	data_dictionary["degradation_matrix"] = degradation_matrix

	data_dictionary["binding_parameter_dictionary"] = binding_parameter_dictionary
	data_dictionary["control_parameter_dictionary"] = control_parameter_dictionary
	data_dictionary["parameter_name_mapping_array"] = parameter_name_mapping_array
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
