module SimilarityMetrics
import FastLevenshtein

    """
        get_similarity(input1, input2, similarity_func)

    Returns similarity between two inputs using the given simliarity function if the similarity function is supported
    """
    function get_similarity(input1, input2, similarity_func::String)
        supported_metrics = Dict{String, Any}(
            "levenshtein" => levenshtien_similarity
        )

        # Throw an error if the given similarity_func is not supported
        if !(similarity_func in keys(supported_metrics))
            error("SimilarityMetrics: $(similarity_func) is not supported. Please choose from $(join(collect(keys(supported_metrics)), ',')).")
        end

        return supported_metrics[similarity_func](input1, input2)
    end

    """
        levenshtien_similarity(input1::String, input2::String)
    
    Returns inverse normalized levenshtein distanve between input1 and input2
    """
    function levenshtien_similarity(input1::String, input2::String)
        dist = FastLevenshtein.fastlevenshtein(input1, input2)
        max_len = maximum(length, [input1, input2])
        if max_len == 0
            return 0.0
        end
        return 1 - (dist / max_len)
    end

end