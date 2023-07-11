module Parsers
export NcbiGenomeAnnotationParser

    module NcbiGenomeAnnotationParser

        """
            parse(path::String)

        Return dictionary containing annotated genome data parsed from a NCBI exported .txt file.
        """
        function parse(path::String)
            retval = Dict{String, Dict}()

            # Check to make sure the file has a .txt extension
            extension = split(path, '.')[end]
            if extension != "txt"
                error("NcbiGenomeAnnotationParser: Incompatible file format found ($(extension)). Only txt file extensions are supported.")
            end

            # Open and read file
            f = open(path, "r")
            raw_txt = read(f, String)
            close(f)

            # Split the raw text on the `>` character to get the information pertaining to each gene
            raw_genes = split(raw_txt, '>')

            # Iterate through each gene and store annotation data into a dictionary
            for gene_str in raw_genes[2:end]

                if gene_str != ""
                    # Split up the raw gene strings by newlines and strip extra spaces
                    lines = [strip(line, ' ') for line in split(gene_str, "\n") if line != ""]

                    # Use the first line for metadata
                    metadata_pairs = split(lines[1], " [")
                    gene_id = metadata_pairs[1]
                    retval[gene_id] = Dict{String, Any}()
                    for metadata_pair in metadata_pairs[2:end]
                        kv_split = split(metadata_pair, '=')
                        retval[gene_id][kv_split[1]] = join(kv_split[2:end], '=')[1:end-1]
                    end

                    # Construct the gene translation by appending the rest of the lines together
                    gene_translation = ""
                    for line in lines[2:end]
                        gene_translation = gene_translation * line
                    end
                    retval[gene_id]["gene_translation"] = gene_translation
                    retval[gene_id]["genome"] = [String(split(path, Base.Filesystem.path_separator)[end])]
                    retval[gene_id]["group"] = split(path, Base.Filesystem.path_separator)[end-1]
                end
                
            end

            return retval

        end

    end


end