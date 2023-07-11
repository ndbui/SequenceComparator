SequenceComparator
================
[![Build Status](https://github.com/ndbui/SequenceComparator.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ndbui/SequenceComparator.jl/actions/workflows/CI.yml?query=branch%3Amain)


This project will primarily be used in the aid of studying how certain genetic elements influence the sensitivity of bacterial strains (Streptococci in this instance) to certain small molecules by comparing groups of genomes, where one group is known to have its growth inhibited by the small molecule, and identifying genetic elements that are shared/distinct between them.



## Using this project

The input of this project are NCBI annotated genomes exported in .txt format. Place each exported annotated genome in the `src/input/` directory of this project where each genome is in its own .txt file
