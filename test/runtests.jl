using SequenceComparator
using Test

@testset "SequenceComparator" begin

    # Test that application properties load correctly
    @test SequenceComparator.load_properties("config/test-config.json")["similarity-metric"]["threshold"] == 0.5
    @test SequenceComparator.load_properties("config/test-config.json")["similarity-metric"]["metric"] == "levenshtein"

end

@testset "NcbiGenomeAnnotationParser" begin
    # Test that the ncbi genome annotation parser won't parse non txt files
    @test_throws ErrorException Parsers.NcbiGenomeAnnotationParser.parse("input/test_genome.pdf")

    # Test that the ncbi genome annotation parser is parsing all the genes correctly
    genome = Parsers.NcbiGenomeAnnotationParser.parse("input/test_genome.txt")
    gene1 = Dict{String, String}("gene" => "test", "locus_tag" => "test_0001", "protein" => "test protein", "gene_translation" => "AAAAAAAAAAGGGGGGGGGGAAAAAAAAAAAAAAAAAAAA")
    gene2 = Dict{String, String}("locus_tag" => "test_0002", "protein" => "test protein #2", "gene_translation" => "TTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCC")
    @test length(genome) == 2
    @test genome["gene1|identifier"] == gene1
    @test genome["gene2|identifier"] == gene2

end

@testset "SimilarityMetrics" begin
    # Test that the SimilarityMetric module throws an error on unsupported SimilarityMetrics
    @test_throws ErrorException SimilarityMetrics.get_similarity("test", "testing", "unsupported_metric")
end

@testset "LevenshteinMetric" begin

    # Test that levenshtein metric is symmetric and outputting correct value
    @test SimilarityMetrics.get_similarity("test", "besttest", "levenshtein") == 0.5
    @test SimilarityMetrics.get_similarity("besttest", "test", "levenshtein") == 0.5
    @test round(SimilarityMetrics.get_similarity("ATCGTAG", "AGTACCT", "levenshtein"),digits=3) == round(2/7, digits=3)
    @test round(SimilarityMetrics.get_similarity("AGTACCT", "ATCGTAG", "levenshtein"),digits=3) == round(2/7, digits=3)
    @test SimilarityMetrics.get_similarity("ATCGTAG", "ATCGTAG", "levenshtein") == 1
    @test round(SimilarityMetrics.get_similarity("ATCGTAG", "ATCTTAG", "levenshtein"),digits=3) == round(6/7, digits=3)
    @test round(SimilarityMetrics.get_similarity("ATCTTAG", "ATCGTAG", "levenshtein"),digits=3) == round(6/7, digits=3)
    @test SimilarityMetrics.get_similarity("a", "sdfg", "levenshtein") == 0
    @test SimilarityMetrics.get_similarity("sdfg", "a", "levenshtein") == 0
end