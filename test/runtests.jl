using SequenceComparator
using Test

@testset "SequenceComparator" begin

    # Test that application properties load correctly
    @test SequenceComparator.load_properties("config/test-config.json")["similarity-metric"]["threshold"] == 0.5
    @test SequenceComparator.load_properties("config/test-config.json")["similarity-metric"]["metric"] == "levenshtein"

    # Test that loading groups of input files is working correctly
    group, genomes = SequenceComparator.load_group("input/sensitive")
    @test length(group) == 4
    @test group["gene1|identifier"]["genome"] == ["test_genome1.txt","test_genome2.txt"]
    @test group["gene1|identifier"]["locus_tag"] == "test_0001"
    @test group["gene2|identifier"]["locus_tag"] == "test_0002"
    @test group["gene3|identifier"]["locus_tag"] == "test_0003"
    @test group["gene4|identifier"]["locus_tag"] == "test_0004"
    @test genomes == ["test_genome1.txt","test_genome2.txt"]

    # Test that getting all genes from a target genome works
    genes = SequenceComparator.get_genome_genes(group, "test_genome1.txt")
    @test length(genes) == 2

    #Test that the gene similarity matrix within a group is being calculated correctly
    @test SequenceComparator.get_gene_group_similarity_matrix(group, group, similarity_metric="levenshtein") == [1.0 0.275 0.0 0.0; 0.275 1.0 0.09999999999999998 0.09999999999999998; 0.0 0.09999999999999998 1.0 0.975; 0.0 0.09999999999999998 0.975 1.0]
    elements_of_interest = Dict{String,Any}("gene1|identifier" => group["gene1|identifier"], "gene2|identifier" => group["gene2|identifier"])
    @test SequenceComparator.get_gene_group_similarity_matrix(elements_of_interest, group, similarity_metric="levenshtein") == [1.0 0.0; 0.275 0.09999999999999998; 0.0 1.0; 0.0 0.975]

    # Test getting common seqeunce elements
    common_elements_seq = SequenceComparator.get_common_seq_elements(SequenceComparator.get_genome_genes(group, "test_genome1.txt"), SequenceComparator.get_genome_genes(group, "test_genome2.txt"))
    @test length(common_elements_seq) == 3
    @test sort(collect(keys(common_elements_seq))) == ["gene1|identifier", "gene2|identifier", "gene4|identifier"]

    # Test getting common group elements
    nonsensitive_group, nonsensitive_genomes = SequenceComparator.load_group("input/nonsensitive")
    common_elements_group = SequenceComparator.get_common_group_elements(nonsensitive_group, nonsensitive_genomes)
    @test length(common_elements_group) == 6
    @test sort(collect(keys(common_elements_group))) == [
        "gene10|identifier",
        "gene12|identifier",
        "gene13|identifier",
        "gene6|identifier",
        "gene7|identifier",
        "gene9|identifier"
    ]

end

@testset "NcbiGenomeAnnotationParser" begin
    # Test that the ncbi genome annotation parser won't parse non txt files
    @test_throws ErrorException Parsers.NcbiGenomeAnnotationParser.parse("input/test_genome.pdf")

    # Test that the ncbi genome annotation parser is parsing all the genes correctly
    genome = Parsers.NcbiGenomeAnnotationParser.parse("input/sensitive/test_genome1.txt")
    gene1 = Dict{String, Any}("gene" => "test", "locus_tag" => "test_0001", "protein" => "test protein", "gene_translation" => "AAAAAAAAAAGGGGGGGGGGAAAAAAAAAAAAAAAAAAAA", "genome" => ["test_genome1.txt"], "group" => "sensitive")
    gene2 = Dict{String, Any}("locus_tag" => "test_0002", "protein" => "test protein #2", "gene_translation" => "TTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCC", "genome" => ["test_genome1.txt"], "group" => "sensitive")
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

@testset "Utils" begin
    dict1 = Dict{String, String}("a" => "b", "c" => "d")
    dict2 = Dict{String, String}("e" => "f", "g" => "h")
    dict3 = Dict{String, String}("a" => "z", "i" => "j")

    # Test that appending dictionaries is working correctly and that it ignores existing keys if seen again
    SequenceComparator.append_dict!(dict1, dict2)
    @test dict1 == Dict{String, String}("a" => "b", "c" => "d", "e" => "f", "g" => "h")
    SequenceComparator.append_dict!(dict1, dict3)
    @test dict1 == Dict{String, String}("a" => "b", "c" => "d", "e" => "f", "g" => "h", "i" => "j")

    genome = Parsers.NcbiGenomeAnnotationParser.parse("input/sensitive/test_genome1.txt")
    test_dict = Dict{String, Any}()
    SequenceComparator.append_dict!(test_dict, genome)
    @test length(test_dict) == 2
end