using SequenceComparator
using Test

@testset "SequenceComparator.jl" begin
    @test SequenceComparator.test_function() == "Hello world!"
end
