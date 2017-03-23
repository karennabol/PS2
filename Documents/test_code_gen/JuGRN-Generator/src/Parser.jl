function parse_grn_file(path_to_model_file::AbstractString)

  # We are going to load the sentences in the file into a vector
  # if not a valid model file, then throw an error -
  sentence_vector = VGRNSentence[]
  tmp_array::Array{AbstractString} = AbstractString[]

  try

    # Open the model file, and read each line into a vector -
    open(path_to_model_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    for sentence in tmp_array

      # Ok, so now we have the array for sentences -
      grn_sentence = VGRNSentence()
      grn_sentence.original_sentence = sentence

      # split the sentence -
      split_array = split(sentence," ")

      # sentence_actor_clause::AbstractString
      # sentence_action_clause::AbstractString
      # sentence_target_clause::AbstractString
      # sentence_delimiter::Char
      grn_sentence.sentence_actor_clause = split_array[1]
      grn_sentence.sentence_action_clause = split_array[2]
      grn_sentence.sentence_target_clause = split_array[3]
      grn_sentence.sentence_delimiter = ' '


      # add sentence to sentence_vector =
      push!(sentence_vector,grn_sentence)
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  # return - (I know we don't need the return, but I *** hate *** the normal Julia convention)
  return sentence_vector
end
