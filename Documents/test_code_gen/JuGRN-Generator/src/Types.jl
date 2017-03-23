type VGRNSentence

  original_sentence::AbstractString
  sentence_actor_clause::AbstractString
  sentence_action_clause::AbstractString
  sentence_target_clause::AbstractString
  sentence_delimiter::Char

  function VGRNSentence()
    this = new()
  end
end

type ProgramComponent

  filename::AbstractString
  buffer::AbstractString

  function ProgramComponent()
    this = new()
  end

end

type SpeciesObject

  species_type::Symbol
  species_symbol::AbstractString

end

type ConnectionObject

  connection_symbol::AbstractString
  connection_actor_set::Array{SpeciesObject}
  connection_target_set::Array{SpeciesObject}
  connection_type::Symbol

  function ConnectionObject()
    this = new()
  end
end

type ProblemObject

  configuration_dictionary::Dict{AbstractString,Any}
  list_of_species::Array{SpeciesObject}
  list_of_connections::Array{ConnectionObject}

  function ProblemObject()
    this = new()
  end
end
