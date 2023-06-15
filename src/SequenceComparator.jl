module SequenceComparator
export NcbiGenomeAnnotationParser, app_properties
include("parsers/ncbi_genome_annotation_parser.jl")
import JSON

# Load application properties from json file 
app_properties = JSON.parsefile("config/local-config.json")


end
