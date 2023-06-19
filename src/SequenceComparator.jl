module SequenceComparator
export Parsers, SimilarityMetrics, main
include("similarity_metrics.jl")
include("utils.jl")
include("parsers.jl")
import JSON

    """
        load_properties(config_path::String)

    Load application properties from json file.
    """
    function load_properties(config_path::String)
        return JSON.parsefile(config_path)
    end


    """
        load_group(input_path::String)

    Returns dictionary containing all the parsed data in a given dir using the filenames as keys.
    """
    function load_group(input_path::String)
        retval = Dict{Any, Any}()
        for file_path in readdir(input_path, join=true)
            filename = split(file_path, '/')[end]
            if filename[1] != '.'
                retval[filename] = Parsers.NcbiGenomeAnnotationParser.parse(file_path)
            end
        end
        return retval
    end


    """
        get_common_seq_elements(seq1, seq2, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

    Returns dictionary of elements common to both seq1 and seq2 using the similarity_metric provided 
    as the basis for comparison against the acceptance_threshold.

    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - acceptance_threshold: Threshold for the similarity metric on what should be counted as similar items.
    - data_field: Name of the field in the element that contains the data that will be used in the element-wise comparison.
    """
    function get_common_seq_elements(seq1, seq2, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")
        retval = Dict{String, Any}()

        # Get all combination of keys in each sequence using dot product and make each comparison
        # TODO: optimize this loop by adding some parallelization. Maybe multi-threading?
        for comparison_keys in Iterators.product([key for key in keys(seq1)], [key2 for key2 in keys(seq2)])
            similarity = SimilarityMetrics.get_similarity(seq1[comparison_keys[1]][data_field], seq2[comparison_keys[2]][data_field], similarity_metric)
            if similarity >= acceptance_threshold
                retval[comparison_keys[1]] = seq1[comparison_keys[1]]
                retval[comparison_keys[2]] = seq2[comparison_keys[2]]
            end

        end
        return retval
    end


    """
        get_common_group_elements(group, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

    Returns dictionary of elements common between all sequences in the group using the similarity_metric provided 
    as the basis for comparison against the acceptance_threshold.

    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - acceptance_threshold: Threshold for the similarity metric on what should be counted as similar items.
    - data_field: Name of the field in the element that contains the data that will be used in the element-wise comparison.
    """
    function get_common_group_elements(group, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")
        if length(keys(group)) == 1
            return group[collect(keys(group))[1]]
        end

        # Set the common elements to the entire first sequence in the group to begin with
        seq_keys = collect(keys(group))
        retval = group[seq_keys[1]]

        # Iteratively reduce the common elements by going through the other sequences and only keeping elements 
        # that pass the element-wise comparison
        for seq_key in seq_keys[2:end]
            retval = get_common_seq_elements(retval, group[seq_key], similarity_metric, acceptance_threshold, data_field)
        end
        return retval

    end


    """
        main(config_path::String)

    Compare two groups of seqeunces and return the individual elements that are common between both groups.
    """
    function main(config_path::String)
        app_properties = load_properties(config_path)

        # Parse the sensitive group
        sensitive_group = load_group("src/input/sensitive")

        # Parse the nonsensitive group
        nonsensitive_group = load_group("src/input/nonsensitive")

        common_elements_nonsensitive = get_common_group_elements(nonsensitive_group)
        common_elements_sensitive = get_common_group_elements(sensitive_group)
        retval = get_common_seq_elements(common_elements_sensitive, common_elements_nonsensitive)
        return retval
    end

end
