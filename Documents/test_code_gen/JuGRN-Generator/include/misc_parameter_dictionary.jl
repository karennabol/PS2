# Alias the misc parameters -
misc_parameter_dictionary = Dict{AbstractString,Float64}()
misc_parameter_dictionary["rnapII_concentration"] = rnapII_concentration  # muM
misc_parameter_dictionary["ribosome_concentration"] = ribosome_concentration # muM
misc_parameter_dictionary["degradation_constant_mRNA"] = degradation_constant_mRNA  # hr^-1
misc_parameter_dictionary["degradation_constant_protein"] = degradation_constant_protein  # hr^-1
misc_parameter_dictionary["kcat_transcription"] = kcat_transcription  # hr^-1
misc_parameter_dictionary["kcat_translation"] = kcat_translation  # hr^-1
misc_parameter_dictionary["maximum_specific_growth_rate"] = maximum_specific_growth_rate  # hr^-1
misc_parameter_dictionary["avg_gene_concentration"] = avg_gene_concentration
misc_parameter_dictionary["saturation_constant_transcription"] = saturation_transcription
misc_parameter_dictionary["saturation_constant_translation"] = saturation_translation
