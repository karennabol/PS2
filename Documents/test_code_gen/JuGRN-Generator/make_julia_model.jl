using ArgParse
using JSON
using Dates

include("./src/Types.jl")
include("./src/Macros.jl")
include("./src/Parser.jl")
include("./src/Problem.jl")
include("./src/Common.jl")
include("./src/strategy/JuliaStrategy.jl")

# Grab the required functions for code generation -
# const parser_function::Function = parse_vff_file

function parse_commandline()
    settings_object = ArgParseSettings()
    @add_arg_table settings_object begin
      "-o"
        help = "Directory where the Julia model files will be written."
        arg_type = AbstractString
        default = "."

      "-m"
        help = "Path to the grn file written in the vgrn format."
        arg_type = AbstractString
        required = true

        "-s"
        help = "Host type: bacteria or mammalian?"
        arg_type = Symbol
        default = :bacteria
    end

    # return a dictionary w/args -
    return parse_args(settings_object)
end


function main()

  # Build the arguement dictionary -
  parsed_args = parse_commandline()

  # Load the statement_vector -
  path_to_model_file = parsed_args["m"]
  statement_vector::Array{VGRNSentence} = parse_grn_file(path_to_model_file)

  # Generate the problem object -
  problem_object = generate_problem_object(statement_vector)

  # initialize the program set -
  component_set = Set{ProgramComponent}()

  # Write the Inputs -
  program_component_inputs = build_inputs_buffer(problem_object)
  push!(component_set,program_component_inputs)

  # Write the data_dictionary -
  host_type = parsed_args["s"]
  program_component_data_dictionary = build_data_dictionary_buffer(problem_object,host_type)
  push!(component_set,program_component_data_dictionary)

  # Write the Kinetics -
  program_component_kinetics = build_kinetics_buffer(problem_object)
  push!(component_set,program_component_kinetics)

  # Write the Control -
  program_component_control = build_control_buffer(problem_object)
  push!(component_set,program_component_control)

  # Write the dilution_matrix --
  program_component_dilution_matrix = generate_dilution_matrix_buffer(problem_object)
  push!(component_set,program_component_dilution_matrix)

  # Write the degradation_matrix --
  program_component_degradation_matrix = generate_degradation_matrix_buffer(problem_object)
  push!(component_set,program_component_degradation_matrix)

  # Write the stoichiometric_matrix --
  program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
  push!(component_set,program_component_stoichiometric_matrix)

  # Dump the component_set to disk -
  path_to_output_file = parsed_args["o"]
  write_program_components_to_disk(path_to_output_file,component_set)

  # Transfer distrubtion jl files to the output -
  transfer_distribution_files("./distribution",path_to_output_file,".jl")

  # transfer the README files -
  path_to_readme_file = splitdir(path_to_model_file)[1]
  transfer_distribution_files("./distribution",path_to_readme_file,".md")

end

main()
