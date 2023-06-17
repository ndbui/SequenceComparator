module SequenceComparator
export Parsers, SimilarityMetrics, load_properties, main
include("parsers.jl")
include("similarity_metrics.jl")
import JSON


    """
        load_properties(config_path::String)

    Load application properties from json file 
    """
    function load_properties(config_path::String)
        return JSON.parsefile(config_path)
    end


    """
        main(config_path::String)

    Compare two groups of seqeunces and return the individual elements that are common between both groups
    """
    function main(config_path::String)
        app_properties = JSON.parsefile(config_path)
    end


end
