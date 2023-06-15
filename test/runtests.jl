using SequenceComparator
using Test

@testset "SequenceComparator.jl" begin

    # Test that application properties load correctly
    @test SequenceComparator.app_properties["similarity-metric"]["threshold"] == 0.5
    @test SequenceComparator.app_properties["similarity-metric"]["metric"] == "levenshtein"

end

@testset "parsers/ncbi_genome_annotation_parser.jl" begin
    # Test that the ncbi genome annotation parser won't parse non txt files
    @test_throws ErrorException NcbiGenomeAnnotationParser.parse("input/test_genome.pdf")

    # Test that the ncbi genome annotation parser is parsing all the genes correctly
    genome = NcbiGenomeAnnotationParser.parse("input/test_genome.txt")
    gene1 = Dict{String, String}("gene" => "test", "locus_tag" => "test_0001", "protein" => "test protein", "gene_translation" => "AAAAAAAAAAGGGGGGGGGGAAAAAAAAAAAAAAAAAAAA")
    gene2 = Dict{String, String}("locus_tag" => "test_0002", "protein" => "test protein #2", "gene_translation" => "TTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCC")
    @test length(genome) == 2
    @test genome["gene1|identifier"] == gene1
    @test genome["gene2|identifier"] == gene2

end