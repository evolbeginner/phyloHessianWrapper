#! /bin/env julia


################################################
function fasta_to_phylip(input_file::String)
    headers = String[]
    sequences = String[]
    current_seq = ""

    # Read and parse FASTA
    open(input_file, "r") do file
        for line in eachline(file)
            line = strip(line)
            if startswith(line, '>')
                push!(headers, strip(line[2:end]))
                if !isempty(current_seq)
                    push!(sequences, current_seq)
                    current_seq = ""
                end
            else
                current_seq *= line
            end
        end
    end
    push!(sequences, current_seq)  # Add last sequence

    # Validate sequence lengths
    seq_lengths = map(length, sequences)
    if length(unique(seq_lengths)) != 1
        error("All sequences must be the same length for PHYLIP format.")
    end

    # Write PHYLIP file
    output_file = input_file * ".phy"
    open(output_file, "w") do file
        println(file, " ", length(sequences), " ", seq_lengths[1])
        for (header, seq) in zip(headers, sequences)
            padded_header = rpad(header[1:min(end, 10)], 10)
			println(padded_header)
            println(file, padded_header, " ", seq)
        end
    end

    println("✅ PHYLIP file created: $output_file")
end

# Run from command line
if length(ARGS) == 1
    fasta_to_phylip(ARGS[1])
else
    println("Usage: julia MFAtoPHY.jl <input.fasta>")
end


