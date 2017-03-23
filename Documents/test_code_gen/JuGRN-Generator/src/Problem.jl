function generate_problem_object(statement_vector::Array{VGRNSentence})

  # Initilize an empty problem object -
  problem_object::ProblemObject = ProblemObject()

  # Load the JSON configuration file -
  config_dict = JSON.parsefile("./config/Configuration.json")
  problem_object.configuration_dictionary = config_dict

  # construct the array of species -
  species_array::Array{SpeciesObject} = build_species_list(statement_vector)

  # construct the array of reactions -
  connection_array::Array{ConnectionObject} = build_connection_list(statement_vector,config_dict)

  # set data on problem_object -
  problem_object.list_of_species = species_array
  problem_object.list_of_connections = connection_array

  # return#the problem_object -
  return problem_object
end

function build_connection_list(statement_vector::Array{VGRNSentence},configuration_dictionary::Dict{String,Any})

  # initialize -
  list_of_connections = ConnectionObject[]

  # iterate through the statement vector -
  for vgrn_sentence in statement_vector

    # Create conenction object -
    connection_object = ConnectionObject()

    # Who are my actors?
    list_of_actor_symbols = AbstractString[]
    actor_string = vgrn_sentence.sentence_actor_clause
    recursive_species_parser!(reverse(collect(actor_string)),list_of_actor_symbols)
    actor_set::Array{SpeciesObject} = species_object_factory(list_of_actor_symbols)

    # Who are my targets?
    list_of_target_symbols = AbstractString[]
    target_string = vgrn_sentence.sentence_target_clause
    recursive_species_parser!(reverse(collect(target_string)),list_of_target_symbols)
    target_set::Array{SpeciesObject} = species_object_factory(list_of_target_symbols)

    # set them on the connection object -
    # ok - do we have a compound set of actors -
    if (length(actor_set) == 1)

      actor_object = actor_set[1]
      connection_object.connection_symbol = actor_object.species_symbol
    else


      local_buffer = ""
      number_of_actors = length(actor_set)
      for (index,actor_object::SpeciesObject) in enumerate(actor_set)

        local_buffer *= "$(actor_object.species_symbol)"
        if (index <= number_of_actors - 1)
          local_buffer *= "_"
        end
      end

      connection_object.connection_symbol = local_buffer
    end


    connection_object.connection_actor_set = actor_set
    connection_object.connection_target_set = target_set

    # create a set for list_of_induction_synonyms -
    induction_synonyms = Set{String}()
    list_of_induction_synonyms = configuration_dictionary["list_of_induction_synonyms"]
    for (index,local_dictionary) in enumerate(list_of_induction_synonyms)

      # grab the symbol -
      symbol = local_dictionary["symbol"]
      push!(induction_synonyms,symbol)
    end

    # create a set of list_of_repression_synonyms -
    repression_synonyms = Set{String}()
    list_of_repression_synonyms = configuration_dictionary["list_of_repression_synonyms"]
    for (index,local_dictionary) in enumerate(list_of_repression_synonyms)

      # grab the symbol -
      symbol = local_dictionary["symbol"]
      push!(repression_synonyms,symbol)
    end

    # What is my type?
    connection_action_string = vgrn_sentence.sentence_action_clause
    if (in(connection_action_string,induction_synonyms))
      connection_object.connection_type = :activate
    elseif (in(connection_action_string,repression_synonyms))
      connection_object.connection_type = :inhibit
    end

    # cache this connection -
    push!(list_of_connections,connection_object)
  end

  # return -
  return list_of_connections
end

function species_object_factory(list_of_symbols::Array{AbstractString})

  species_set = SpeciesObject[]
  for symbol_text in list_of_symbols
    push!(species_set,SpeciesObject(:gene,symbol_text))
  end

  return species_set
end

function build_species_list(statement_vector::Array{VGRNSentence})

  list_of_symbols = AbstractString[]
  for vgrn_sentence in statement_vector

    # build species string -
    species_string = vgrn_sentence.sentence_actor_clause*" "*vgrn_sentence.sentence_target_clause
    recursive_species_parser!(reverse(collect(species_string)),list_of_symbols)
  end

  # ok, so we need to convert this to a set, and then convert back to an array (set is unique)
  species_set = Set{AbstractString}()
  for symbol in list_of_symbols
    push!(species_set,symbol)
  end

  # unique list -
  unique_list_of_species = AbstractString[]
  for symbol in species_set
    push!(unique_list_of_species,symbol)
  end

  # sort -
  sort!(unique_list_of_species)

  # ok, so we have a sorted list of *genes* (that is what is in the file)
  # create a list of genes, mRNA and protein -
  list_of_species_objects = SpeciesObject[]

  # build the species objects -
  for species_symbol in unique_list_of_species

    # objects -
    gene_object = SpeciesObject(:gene,species_symbol)

    # mRNA -
    mRNA_symbol = "mRNA_"*species_symbol
    mRNA_object = SpeciesObject(:mrna,mRNA_symbol)

    # protein -
    protein_symbol = "protein_"*species_symbol
    protein_object = SpeciesObject(:protein,protein_symbol)

    # push -
    push!(list_of_species_objects,gene_object)
    push!(list_of_species_objects,mRNA_object)
    push!(list_of_species_objects,protein_object)
  end

  # partition by type -
  return partition!(list_of_species_objects)
end

function recursive_species_parser!(sentence_char_array::Array{Char},list_of_symbols::Array{AbstractString})

  if (isempty(sentence_char_array) == true)
    return
  end

  # process each char -
  test_char = pop!(sentence_char_array)
  if (test_char == '(' || test_char == ' ')
    recursive_species_parser!(sentence_char_array,list_of_symbols)

  elseif (isalpha(test_char) == true)

    # cache this letter -
    biological_symbol_cache = Char[]
    push!(biological_symbol_cache,test_char)

    # When should we stop?
    stop_set = Set{Char}()
    push!(stop_set,' ')
    push!(stop_set,',')
    push!(stop_set,')')
    push!(stop_set,'\n')
    push!(stop_set,'\r')

    # ok, so lets read until we hit a space or a comma -
    if (isempty(sentence_char_array))
      stop_flag = true
    else
      stop_flag = false
    end

    while (stop_flag == false)

      if (isempty(sentence_char_array) == false)
        next_test_char = pop!(sentence_char_array)

        if (in(next_test_char,stop_set) == false)
          push!(biological_symbol_cache,next_test_char)
        else
          stop_flag = true
        end
      else
        stop_flag = true
      end
    end

    # Store the symbol -
    if (isempty(biological_symbol_cache) == false)
      push!(list_of_symbols,String(biological_symbol_cache))
    end

    @show sentence_char_array

    # ok, so we should be ready for another dive -
    recursive_species_parser!(sentence_char_array,list_of_symbols)
  end
end
