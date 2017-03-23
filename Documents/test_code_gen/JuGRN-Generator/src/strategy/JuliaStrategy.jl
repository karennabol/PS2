function build_function_header_buffer(comment_dictionary)

  # initialize -
  buffer = ""

  # get some data from the comment_dictionary -
  function_name = comment_dictionary["function_name"]
  function_description = comment_dictionary["function_description"]
  input_arg_array = comment_dictionary["input_args"]
  output_arg_array = comment_dictionary["output_args"]

  buffer *= "# ----------------------------------------------------------------------------------- #\n"
  buffer *= "# Function: $(function_name)\n"
  buffer *= "# Description: $(function_description)\n"
  buffer *= "# Generated on: $(now())\n"
  buffer *= "#\n"
  buffer *= "# Input arguments:\n"

  for argument_dictionary in input_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end

  buffer *= "#\n"
  buffer *= "# Output arguments:\n"
  for argument_dictionary in output_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end
  buffer *= "# ----------------------------------------------------------------------------------- #\n"

  # return the buffer -
  return buffer
end

function build_copyright_header_buffer(problem_object::ProblemObject)

  # What is the current year?
  current_year = string(Dates.year(now()))

  # Get comment data from
  buffer = ""
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  buffer*= "# Copyright (c) $(current_year) Varnerlab\n"
  buffer*= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
  buffer*= "# Cornell University, Ithaca NY 14850\n"
  buffer*= "#\n"
  buffer*= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
  buffer*= "# of this software and associated documentation files (the \"Software\"), to deal\n"
  buffer*= "# in the Software without restriction, including without limitation the rights\n"
  buffer*= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
  buffer*= "# copies of the Software, and to permit persons to whom the Software is\n"
  buffer*= "# furnished to do so, subject to the following conditions:\n"
  buffer*= "#\n"
  buffer*= "# The above copyright notice and this permission notice shall be included in\n"
  buffer*= "# all copies or substantial portions of the Software.\n"
  buffer*= "#\n"
  buffer*= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
  buffer*= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  buffer*= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
  buffer*= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  buffer*= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
  buffer*= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
  buffer*= "# THE SOFTWARE.\n"
  buffer*= "# ----------------------------------------------------------------------------------- #\n"

  # return -
  return buffer

end

function iterate_binding_control_connection(gene_object::SpeciesObject,list_of_connections::Array{ConnectionObject})

  # What is my gene_symbol -
  gene_symbol = gene_object.species_symbol


  buffer = ""
  for connection_object in list_of_connections

    # actor -
    connection_symbol = connection_object.connection_symbol

    buffer *= "\tn_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"]\n"
    buffer *= "\tK_$(gene_symbol)_$(connection_symbol) = binding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"]\n"

  end

  return buffer

end

function iterate_control_control_connection(gene_object::SpeciesObject,list_of_connections::Array{ConnectionObject})

  # What is my gene_symbol -
  gene_symbol = gene_object.species_symbol


  buffer = ""
  for connection_object in list_of_connections

    # actor -
    connection_symbol = connection_object.connection_symbol
    buffer *= "\tW_$(gene_symbol)_$(connection_symbol) = control_parameter_dictionary[\"W_$(gene_symbol)_$(connection_symbol)\"]\n"
  end

  return buffer

end

function build_control_buffer(problem_object::ProblemObject)

  filename = "Control.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["control_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # extract list of genes -
  list_of_species = problem_object.list_of_species
  list_of_genes = extract_species_of_type(list_of_species,:gene)

  # get list of connections -
  list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function Control(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "\t# initialize the control - \n"
  buffer *= "\tcontrol_array = zeros($(length(list_of_genes)))\n"
  buffer *= "\n"

  buffer *= "\t# Alias the species - \n"
  for (index,species_object) in enumerate(list_of_species)

    # Grab the symbol -
    species_symbol = species_object.species_symbol

    # write the record -
    buffer *= "\t$(species_symbol) = x[$(index)]\n"
  end

  buffer *= "\n"
  buffer *= "\t# Alias the binding parameters - \n"
  buffer *= "\tbinding_parameter_dictionary = data_dictionary[\"binding_parameter_dictionary\"]\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # activating connections -
    buffer *= iterate_binding_control_connection(gene_object,activating_connections)

    # inhibiting_connections -
    buffer *= iterate_binding_control_connection(gene_object,inhibiting_connections)
  end

  buffer *= "\n"
  buffer *= "\t# Alias the control function parameters - \n"
  buffer *= "\tcontrol_parameter_dictionary = data_dictionary[\"control_parameter_dictionary\"]\n"
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # get the RNAp binding symbol out -
    buffer *="\tW_$(gene_symbol)_RNAP = control_parameter_dictionary[\"W_$(gene_symbol)_RNAP\"]\n"

    # activating -
    buffer *= iterate_control_control_connection(gene_object,activating_connections)

    # inhibting -
    buffer *= iterate_control_control_connection(gene_object,inhibiting_connections)
  end

  buffer *= "\n"
  # get list of genes -
  for (gene_index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # generate the binding functions -
    list_of_all_connections = ConnectionObject[]
    append!(list_of_all_connections,activating_connections)
    append!(list_of_all_connections,inhibiting_connections)
    for (index,connection_object) in enumerate(list_of_all_connections)

      # actor -
      actor_list = connection_object.connection_actor_set
      connection_symbol = connection_object.connection_symbol

      buffer *= "\t# Transfer function target:$(gene_symbol) actor:$(connection_symbol)\n"
      buffer *= "\tactor_set_$(gene_symbol)_$(connection_symbol) = [\n"
      for actor_object in actor_list

        actor_symbol = actor_object.species_symbol
        buffer *= "\t\tprotein_$(actor_symbol)\n"
      end

      buffer *="\t]\n"
      buffer *="\tactor = prod(actor_set_$(gene_symbol)_$(connection_symbol))\n"
      buffer *= "\tb_$(gene_symbol)_$(connection_symbol) = (actor^(n_$(gene_symbol)_$(connection_symbol)))/"
      buffer *= "(K_$(gene_symbol)_$(connection_symbol)^(n_$(gene_symbol)_$(connection_symbol))+actor^(n_$(gene_symbol)_$(connection_symbol)))\n"
      buffer *= "\n"
    end

    buffer *= "\t# Control function for $(gene_symbol) - \n"
    buffer *= "\tcontrol_array[$(gene_index)] = ("
    numerator = ""

    if (isempty(activating_connections) == true)
      buffer *= "W_$(gene_symbol)_RNAP"
    else

      numerator *= "W_$(gene_symbol)_RNAP+"
      for connection_object in activating_connections
        # actor -
        connection_symbol = connection_object.connection_symbol
        numerator *= "W_$(gene_symbol)_$(connection_symbol)*b_$(gene_symbol)_$(connection_symbol)+"
      end
      buffer *= numerator[1:end-1]
    end

    buffer *= ")/(1+W_$(gene_symbol)_RNAP"

    if (isempty(activating_connections) == false)
      demoninator = ""
      for connection_object in activating_connections
        # actor -
        connection_symbol = connection_object.connection_symbol
        demoninator *= "+W_$(gene_symbol)_$(connection_symbol)*b_$(gene_symbol)_$(connection_symbol)"
      end

      buffer *= demoninator[1:end]
    end

    # ok - do we have inhibitory statements?
    if (isempty(inhibiting_connections) == true)
      buffer *= ")\n"
      buffer *= "\n"
    else

      demoninator = ""
      for connection_object in inhibiting_connections
        # actor -
        connection_symbol = connection_object.connection_symbol
        demoninator *= "+W_$(gene_symbol)_$(connection_symbol)*b_$(gene_symbol)_$(connection_symbol)"
      end

      buffer *= demoninator[1:end]
      buffer *= ")\n"
      buffer *= "\n"
    end
  end

  buffer *= "\t# return - \n"
  buffer *= "\treturn control_array\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_data_dictionary_buffer(problem_object::ProblemObject,host_flag::Symbol)

  filename = "DataDictionary.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # default parameter default dictionary -
  parameter_value_default_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]

  # get list of species from the po -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function DataDictionary(time_start::Float64,time_stop::Float64,time_step_size::Float64)\n"
  buffer *= "\n"
  buffer *= "\t# stoichiometric_matrix and dilution_matrix - \n"
  buffer *= "\tstoichiometric_matrix = readdlm(\"./Network.dat\")\n"
  buffer *= "\tdilution_matrix = readdlm(\"./Dilution.dat\")\n"
  buffer *= "\tdegradation_matrix = readdlm(\"./Degradation.dat\")\n"
  buffer *= "\n"
  buffer *= "\t# number of states, and rates - \n"
  buffer *= "\t(number_of_states,number_of_rates) = size(stoichiometric_matrix)\n"
  buffer *= "\n"
  buffer *= "\t# array of gene lengths - \n"
  buffer *= "\tgene_coding_length_array = [\n"

  # write out the length of genes -
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :gene)
      buffer *= "\t\t15000\t;\t# $(index)\t$(species_symbol)\n"
    end
  end
  buffer *= "\t]\n"
  buffer *= "\n"
  buffer *= "\t# array of mRNA coding lengths - \n"
  buffer *= "\tmRNA_coding_length_array = [\n"

  # write out the length of genes -
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :mrna)
      buffer *= "\t\tgene_coding_length_array[$(counter)]\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
      counter = counter + 1
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"

  buffer *= "\t# array of mRNA coding lengths - \n"
  buffer *= "\tprotein_coding_length_array = [\n"

  # write out the length of genes -
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :protein)
      buffer *= "\t\tround((0.33)*mRNA_coding_length_array[$(counter)])\t;\t# $(index)\t$(counter)\t$(species_symbol)\n"
      counter = counter + 1
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"

  if host_flag == :bacteria
    buffer *= @include_function("txtl_constants_ecoli","\t")
  else
    buffer *= @include_function("txtl_constants_hl60","\t")
  end
  buffer *= "\n"

  buffer *= "\t# initial condition array - \n"
  buffer *= "\tinitial_condition_array = [\n"

  # write out species -
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    if (species_type == :gene)
      buffer *= "\t\tavg_gene_concentration\t;\t# $(index)\t$(species_symbol)\n"
    elseif (species_type == :mrna || species_type == :protein)
      buffer *= "\t\t0.0\t;\t# $(index)\t$(species_symbol)\n"
    end
  end

  buffer *= "\t]\n"
  buffer *= "\n"

  # Setup parameters dictionaries -
  list_of_genes = extract_species_of_type(list_of_species,:gene)

  # get list of connections -
  list_of_connections::Array{ConnectionObject} = problem_object.list_of_connections
  buffer *= "\tbinding_parameter_dictionary = Dict{AbstractString,Float64}()\n"
  default_TF_KD_parameter = parameter_value_default_dictionary["default_TF_KD_parameter"]
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # Iterate through all connections -
    list_of_all_connections = ConnectionObject[]
    append!(list_of_all_connections,activating_connections)
    append!(list_of_all_connections,inhibiting_connections)
    for connection_object in list_of_all_connections

      # connection -
      connection_symbol = connection_object.connection_symbol

      # write the line -
      buffer *= "\tbinding_parameter_dictionary[\"n_$(gene_symbol)_$(connection_symbol)\"] = 1.0\n"
      buffer *= "\tbinding_parameter_dictionary[\"K_$(gene_symbol)_$(connection_symbol)\"] = $(default_TF_KD_parameter)\n"

    end
  end

  buffer *= "\n"
  buffer *= "\t# Alias the control function parameters - \n"
  buffer *= "\tcontrol_parameter_dictionary = Dict{AbstractString,Float64}()\n"
  default_W_gene_symbol_RNAP = parameter_value_default_dictionary["default_background_expression_parameter"]
  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # generate an RNAP term -
    buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_RNAP\"] = $(default_W_gene_symbol_RNAP)\n"

    # process all connections -
    list_of_all_connections = ConnectionObject[]
    append!(list_of_all_connections,activating_connections)
    append!(list_of_all_connections,inhibiting_connections)
    for connection_object in list_of_all_connections

      # actor -
      connection_symbol = connection_object.connection_symbol
      buffer *= "\tcontrol_parameter_dictionary[\"W_$(gene_symbol)_$(connection_symbol)\"] = 1.0\n"
    end
  end

  # put the misc dictionary -
  buffer *= "\n"
  buffer *= @include_function("misc_parameter_dictionary","\t")
  buffer *= "\n"

  # parameter name mapping -
  buffer *= "\n"
  buffer *= "\t# Parameter name index array - \n"
  name_parameter_mapping_buffer = generate_parameter_name_mapping(list_of_genes,list_of_connections::Array{ConnectionObject})
  buffer *= name_parameter_mapping_buffer

  # return block -
  buffer *= "\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = Dict{AbstractString,Any}()\n"
  buffer *= "\tdata_dictionary[\"number_of_states\"] = number_of_states\n"
  buffer *= "\tdata_dictionary[\"initial_condition_array\"] = initial_condition_array\n"
  buffer *= "\tdata_dictionary[\"average_transcript_length\"] = average_transcript_length\n"
  buffer *= "\tdata_dictionary[\"average_protein_length\"] = average_protein_length\n"
  buffer *= "\tdata_dictionary[\"gene_coding_length_array\"] = gene_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"mRNA_coding_length_array\"] = mRNA_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"protein_coding_length_array\"] = protein_coding_length_array\n"
  buffer *= "\tdata_dictionary[\"rnapII_concentration\"] = rnapII_concentration  # \muM \n"
  buffer *= "\tdata_dictionary[\"ribosome_concentration\"] = ribosome_concentration # \muM \n"
  buffer *= "\tdata_dictionary[\"degradation_constant_mRNA\"] = degradation_constant_mRNA  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"degradation_constant_protein\"] = degradation_constant_protein  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"kcat_transcription\"] = kcat_transcription  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"kcat_translation\"] = kcat_translation  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"maximum_specific_growth_rate\"] = maximum_specific_growth_rate  # hr^-1 \n"
  buffer *= "\tdata_dictionary[\"death_rate_constant\"] = death_rate_constant \n"
  buffer *= "\tdata_dictionary[\"avg_gene_concentration\"] = avg_gene_concentration \n"
  buffer *= "\tdata_dictionary[\"saturation_constant_transcription\"] = saturation_transcription \n"
  buffer *= "\tdata_dictionary[\"saturation_constant_translation\"] = saturation_translation \n"
  buffer *= "\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"dilution_matrix\"] = dilution_matrix\n"
  buffer *= "\tdata_dictionary[\"degradation_matrix\"] = degradation_matrix\n"
  buffer *= "\n"
  buffer *= "\tdata_dictionary[\"binding_parameter_dictionary\"] = binding_parameter_dictionary\n"
  buffer *= "\tdata_dictionary[\"control_parameter_dictionary\"] = control_parameter_dictionary\n"
  buffer *= "\tdata_dictionary[\"parameter_name_mapping_array\"] = parameter_name_mapping_array\n"
  buffer *= "\tdata_dictionary[\"misc_parameter_dictionary\"] = misc_parameter_dictionary\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_kinetics_buffer(problem_object::ProblemObject)

  filename = "Kinetics.jl"

  # get list of species from the po -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer

  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_transcription_rates"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function calculate_transcription_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_genes = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :gene)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_genes = number_of_genes + 1
    end
  end
  buffer *="\n"
  buffer *="\t# Initialize the transcription rate - \n"
  buffer *="\ttranscription_rate_array = zeros($(number_of_genes))\n"
  buffer *="\tKSAT = data_dictionary[\"saturation_constant_transcription\"]\n"
  buffer *="\tkcat_transcription = data_dictionary[\"kcat_transcription\"]\n"
  buffer *="\trnapII_concentration = data_dictionary[\"rnapII_concentration\"]\n"
  buffer *="\taverage_transcript_length = data_dictionary[\"average_transcript_length\"]\n"
  buffer *="\tgene_coding_length_array = data_dictionary[\"gene_coding_length_array\"]\n"

  buffer *="\n"
  buffer *="\t# Populate the transcription rate array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :gene)
      buffer *= "\t# Gene: $(species_symbol)\n"
      buffer *= "\tgene_length = gene_coding_length_array[$(index)]\n"
      buffer *= "\tscale_factor = (average_transcript_length/gene_length)\n"
      buffer *= "\ttranscription_rate_array[$(counter)] = scale_factor*kcat_transcription*(rnapII_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
      buffer *= "\n"
      counter = counter + 1
    end
  end

  buffer *="\n"
  buffer *= "\t# return transcription_rate_array - \n"
  buffer *= "\treturn transcription_rate_array\n"

  buffer *= "end\n"
  buffer *= "\n"


  # calculate_background_transcription_rates -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_background_transcription_rates"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
  buffer *= function_comment_buffer
  buffer *= "function calculate_background_transcription_rates(t::Float64,x::Array{Float64,1},transcription_rate_array::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\treturn zeros(length(x))\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"


  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_translation_rates"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
  buffer *= function_comment_buffer
  buffer *= "function calculate_translation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_mRNA = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_mRNA = number_of_mRNA + 1
    end
  end
  buffer *="\n"
  buffer *="\t# Initialize the translation rate - \n"
  buffer *="\ttranslation_rate_array = zeros($(number_of_mRNA))\n"
  buffer *="\tKSAT = data_dictionary[\"saturation_constant_translation\"]\n"
  buffer *="\tkcat_translation = data_dictionary[\"kcat_translation\"]\n"
  buffer *="\tribosome_concentration = data_dictionary[\"ribosome_concentration\"]\n"
  buffer *="\taverage_protein_length = data_dictionary[\"average_protein_length\"]\n"
  buffer *="\tprotein_coding_length_array = data_dictionary[\"protein_coding_length_array\"]\n"
  buffer *="\n"
  buffer *="\t# Populate the translation rate array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\t# Transcript: $(species_symbol)\n"
      buffer *= "\tprotein_length = protein_coding_length_array[$(counter)]\n"
      buffer *= "\tscale_factor = (average_protein_length/protein_length)\n"
      buffer *= "\ttranslation_rate_array[$(counter)] = scale_factor*kcat_translation*(ribosome_concentration)*(($(species_symbol))/(KSAT+$(species_symbol)))\n"
      buffer *= "\n"
      counter = counter + 1
    end
  end

  buffer *="\n"
  buffer *= "\t# return translation array - \n"
  buffer *= "\treturn translation_rate_array\n"
  buffer *= "end\n"
  buffer *= "\n"


  # calculate_mRNA_degradation_rates -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_mRNA_degradation_rates"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
  buffer *= function_comment_buffer
  buffer *= "function calculate_mRNA_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"

  buffer *="\t# Alias the species - \n"
  number_of_mRNA = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_mRNA = number_of_mRNA + 1
    end
  end
  buffer *= "\n"
  buffer *="\t# Initialize the degrdation array - \n"
  buffer *="\tdegradation_rate_array = zeros($(number_of_mRNA))\n"
  buffer *="\tmRNA_degrdation_constant = data_dictionary[\"degradation_constant_mRNA\"]\n"
  buffer *= "\n"
  buffer *="\t# Calculate the degradation_rate_array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :mrna)
      buffer *= "\tdegradation_rate_array[$(counter)] = (mRNA_degrdation_constant)*$(species_symbol)\n"
      counter = counter + 1
    end
  end
  buffer *= "\n"
  buffer *= "\t# return the degrdation rate array - \n"
  buffer *= "\treturn degradation_rate_array\n"
  buffer *= "end\n"
  buffer *= "\n"

  # calculate_protein_degradation_rates
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["calculate_protein_degradation_rates"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)
  buffer *= function_comment_buffer
  buffer *= "function calculate_protein_degradation_rates(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *="\t# Alias the species - \n"
  number_of_proteins = 0
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :protein)
      buffer *= "\t$(species_symbol) = x[$(index)]\n"
      number_of_proteins = number_of_proteins + 1
    end
  end
  buffer *= "\n"
  buffer *="\t# Initialize the degrdation array - \n"
  buffer *="\tdegradation_rate_array = zeros($(number_of_proteins))\n"
  buffer *="\tprotein_degrdation_constant = data_dictionary[\"degradation_constant_protein\"]\n"
  buffer *= "\n"
  buffer *="\t# Calculate the degradation_rate_array - \n"
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    # grab the species -
    species_symbol = species_object.species_symbol
    species_type = species_object.species_type

    # grab -
    if (species_type == :protein)
      buffer *= "\tdegradation_rate_array[$(counter)] = (protein_degrdation_constant)*$(species_symbol)\n"
      counter = counter + 1
    end
  end
  buffer *= "\n"
  buffer *= "\t# return the degrdation rate array - \n"
  buffer *= "\treturn degradation_rate_array\n"
  buffer *= "end\n"
  buffer *= "\n"
  buffer *= "\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function build_inputs_buffer(problem_object::ProblemObject)

  filename = "Inputs.jl"

  # build the header -
  license_header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["input_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # initialize the buffer -
  buffer = ""
  buffer *= license_header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function calculate_input_array(t::Float64,x::Array{Float64,1},data_dictionary::Dict{AbstractString,Any})\n"
  buffer *= "\n"
  buffer *= "\t# Initialize default - \n"
  buffer *= "\tu_array = zeros(length(x))\n"
  buffer *= "\n"
  buffer *= "\t# return - \n"
  buffer *= "\treturn u_array\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer

  # return -
  return (program_component)
end

function generate_parameter_name_mapping(list_of_genes::Array{SpeciesObject},list_of_connections::Array{ConnectionObject})

  # iterate through the list of genes - add names to a string -
  list_of_names = AbstractString[]

  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # grab the list of connections -
    list_of_all_connections = ConnectionObject[]
    append!(list_of_all_connections,activating_connections)
    append!(list_of_all_connections,inhibiting_connections)
    for connection_object in list_of_all_connections

      # actor -
      connection_symbol = connection_object.connection_symbol

      tmp_name = "n_$(gene_symbol)_$(connection_symbol)"
      push!(list_of_names,tmp_name)

      tmp_name = "K_$(gene_symbol)_$(connection_symbol)"
      push!(list_of_names,tmp_name)
    end
  end

  for (index,gene_object) in enumerate(list_of_genes)

    # get gene symbol -
    gene_symbol = gene_object.species_symbol

    # connections -
    activating_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:activate)
    inhibiting_connections = is_species_a_target_in_connection_list(list_of_connections,gene_object,:inhibit)

    # generate an RNAP term -
    tmp_name = "W_$(gene_symbol)_RNAP"
    push!(list_of_names,tmp_name)

    # grab the list of connections -
    list_of_all_connections = ConnectionObject[]
    append!(list_of_all_connections,activating_connections)
    append!(list_of_all_connections,inhibiting_connections)
    for connection_object in list_of_all_connections

      # actor -
      connection_symbol = connection_object.connection_symbol
      tmp_name = "W_$(gene_symbol)_$(connection_symbol)"
      push!(list_of_names,tmp_name)
    end
  end

  # push some "hard coded" constants onto the end -
  # we should encode this list in a config file ...
  # "rnapII_concentration"	; 				# 1
  # "ribosome_concentration"	;				# 2
  # "degradation_constant_mRNA"	;			# 3
  # "degradation_constant_protein"	;	# 4
  # "kcat_transcription"	;						# 5
  # "kcat_translation"	;							# 6
  # "maximum_specific_growth_rate"	;	# 7
  # "saturation_constant_transcription"	;	# 8
  # "saturation_constant_translation"	;	# 9
  push!(list_of_names,"rnapII_concentration")
  push!(list_of_names,"ribosome_concentration")
  push!(list_of_names,"degradation_constant_mRNA")
  push!(list_of_names,"degradation_constant_protein")
  push!(list_of_names,"kcat_transcription")
  push!(list_of_names,"kcat_translation")
  push!(list_of_names,"maximum_specific_growth_rate")
  push!(list_of_names,"saturation_constant_transcription")
  push!(list_of_names,"saturation_constant_translation")

  buffer = ""
  buffer *= "\tparameter_name_mapping_array = [\n"
  for (index,symbol) in enumerate(list_of_names)

    buffer *= "\t\t\"$(symbol)\"\t;\t# $(index)\n"

  end
  buffer *= "\t]\n"

  # return -
  return buffer
end
