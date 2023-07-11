module SequenceComparator
export Parsers, SimilarityMetrics, main
include("similarity_metrics.jl")
include("utils.jl")
include("parsers.jl")
import JSON, XLSX, Dates, Printf

    """
        load_properties(config_path::String)

    Load application properties from json file.
    """
    function load_properties(config_path::String)
        return JSON.parsefile(config_path)
    end


    """
        load_group(input_path::String)

    Returns dictionary containing all the genes in the group along with metadata taken from each genome file 
    as well as a list of all the genome file names
    """
    function load_group(input_path::String)
        retval = Dict{Any, Any}()
        genomes = []
        for file_path in readdir(input_path, join=true)
            filename = basename(file_path)
            if filename[1] != '.'
                push!(genomes, filename)
                genome = Parsers.NcbiGenomeAnnotationParser.parse(file_path)
                for gene_id in keys(genome)
                    if !(gene_id in keys(retval))
                        retval[gene_id] = genome[gene_id]
                    else
                        push!(retval[gene_id]["genome"], filename)
                    end
                end
            end
                
        end
        return retval, genomes
    end

    """
        get_common_seq_elements(seq1, seq2, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

    Returns dictionary of elements common to both seq1 and seq2 using the similarity_metric provided as the basis for comparison 
    against the acceptance_threshold

    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - acceptance_threshold: Threshold for the similarity metric on what should be counted as similar items.
    - data_field: Name of the field in the element that contains the data that will be used in the element-wise comparison.
    """
    function get_common_seq_elements(seq1, seq2; similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")
        retval = Dict{String, Any}()

        # Get all combination of keys in each sequence using dot product and make each comparison
        comparison_combinations = collect(Iterators.product([key for key in keys(seq1)], [key2 for key2 in keys(seq2)]))
        Threads.@threads for comparison_keys in comparison_combinations
            similarity = SimilarityMetrics.get_similarity(seq1[comparison_keys[1]][data_field], seq2[comparison_keys[2]][data_field], similarity_metric)
            if similarity >= acceptance_threshold
                retval[comparison_keys[1]] = seq1[comparison_keys[1]]
                retval[comparison_keys[2]] = seq2[comparison_keys[2]]
            end

        end
        return retval
    end


    """
        get_common_group_elements(group, genomes, similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

    Returns dictionary of elements common between all sequences in the group using the similarity_metric provided 
    as the basis for comparison against the acceptance_threshold.

    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - acceptance_threshold: Threshold for the similarity metric on what should be counted as similar items.
    - data_field: Name of the field in the element that contains the data that will be used in the element-wise comparison.
    """
    function get_common_group_elements(group, genomes; similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

        if length(genomes) == 1
            return group
        end

        # Set the common elements to the entire first sequence in the group to begin with
        retval = get_genome_genes(group, genomes[1])
    
        # Iteratively reduce the common elements by going through the other sequences and only keeping elements 
        # that pass the element-wise comparison
        for genome_id in genomes[2:end]
            retval = get_common_seq_elements(retval, get_genome_genes(group, genome_id), similarity_metric=similarity_metric, acceptance_threshold=acceptance_threshold, data_field=data_field)
        end
        return retval
    end


    """
        get_gene_group_similarity_matrix(genes, group; similarity_metric="levenshtein", acceptance_threshold=0.8, data_field="gene_translation")

    Returns similarity matrix for each gene compared against every other gene in the specified group
    
    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - data_field: Name of the field in the element that contains the data that will be used in the element-wise comparison.
    """
    function get_gene_group_similarity_matrix(genes, group; similarity_metric="levenshtein", data_field="gene_translation")
        gene_ids = collect(keys(genes))
        group_gene_ids = collect(keys(group))
        symmetric_similarity_metrics = ["levenshtein"]

        similarity_matrix = zeros(length(group_gene_ids), length(gene_ids))
        Threads.@threads for gene_enum in collect(enumerate(gene_ids))
            gene_index = gene_enum[1]
            gene_id = gene_enum[2]
            Threads.@threads for group_gene_enum in collect(enumerate(group_gene_ids))
                group_gene_index = group_gene_enum[1]
                group_gene_id = group_gene_enum[2]
                # Take advantage of symmetric similarity metrics and memoization
                if similarity_metric in symmetric_similarity_metrics
                    sorted_args = sort([group_gene_id, gene_id])
                    arg1 = sorted_args[1]
                    arg2 = sorted_args[2]
                else
                    arg1 = gene_id
                    arg2 = group_gene_id
                end
                similarity_matrix[group_gene_index, gene_index] = SimilarityMetrics.get_similarity(group[arg1][data_field], group[arg2][data_field], similarity_metric)
            end
        end
        return similarity_matrix
    end

    
    """
        calculate_gene_group_identities!(genes, group, group_sequences; similarity_metric="levenshtein", data_field="gene_translation")

    Adds all other genomes each gene is present in, all other similar genes in the group, and the group identity 
    for each gene in the genes argument. Group identity is calculated by [TODO: FINALIZE IDENTITY CALCULATION METHOD]

    
    
    The optional keyword arguments are:
    - similarity_metric: Name of the similarity metric that will be used in the element-wise comparison.
    - acceptance_threshold: Threshold for the similarity metric on what should be counted as similar items.
    - group_identity_threshold: Threshold for which genes are kept based on the group identity calculation.
    """
    function calculate_gene_group_identities!(genes, group, group_sequences; similarity_metric="levenshtein", acceptance_threshold=0.8, group_identity_threshold=0.8)
        # Get the similarity matrix comparing each of the elements of interest against the entire control group
        similarity_matrix =  get_gene_group_similarity_matrix(genes, group,similarity_metric=similarity_metric)
        
        # Add additional metadata to genes including the genomes each gene is present in and all of the other similar genes
        gene_ids = collect(keys(genes))
        gene_ids_index_map = Dict(i => gene_ids[i] for i in eachindex(gene_ids))
        group_ids = collect(keys(group))
        group_ids_index_map = Dict(i => group_ids[i] for i in eachindex(group_ids))
        for y in 1:length(similarity_matrix[1,:])
            genes[gene_ids_index_map[y]]["in_genomes"] = Set()
            genes[gene_ids_index_map[y]]["similar_genes"] = []
            for x in 1:length(similarity_matrix[:,1])
                if similarity_matrix[x,y] >= acceptance_threshold && genes[gene_ids_index_map[y]] != group[group_ids_index_map[x]]
                    # Keep track of all the genomes that contains an element that is similar to the element of interest
                    for genome in group[group_ids_index_map[x]]["genome"]
                        push!(genes[gene_ids_index_map[y]]["in_genomes"], genome)
                    end

                    # Keep track of which elements the element of interest is similar to
                    similar_gene = Dict{String,Any}()
                    similar_gene["gene_id"] = group_ids_index_map[x]
                    similar_gene["protein"] = group[group_ids_index_map[x]]["protein"]
                    similar_gene["similarity"] = similarity_matrix[x,y]
                    similar_gene["genome"] = group[group_ids_index_map[x]]["genome"]
                    push!(genes[gene_ids_index_map[y]]["similar_genes"], similar_gene)
                end
            end
        end

        # Calculate percent identity for each of the elements of interest against the group
        retval = Dict{String, Any}()
        for (gene_id, gene) in pairs(genes)
            if length(gene["similar_genes"]) > 0
                max_similarity_per_genome = Dict{String, Float64}()
                for genome in gene["genome"]
                    max_similarity_per_genome[genome] = 1
                end

                for similar_gene in gene["similar_genes"]
                    for genome in similar_gene["genome"]
                        if !haskey(max_similarity_per_genome, genome) || (haskey(max_similarity_per_genome, genome) && max_similarity_per_genome[genome] < similar_gene["similarity"])
                            max_similarity_per_genome[genome] = similar_gene["similarity"]
                        end
                    end
                end

                genes[gene_id]["group_identity"] = sum(values(max_similarity_per_genome)) / length(group_sequences)
            else
                genes[gene_id]["group_identity"] = 0
            end

            if genes[gene_id]["group_identity"] >= group_identity_threshold
                retval[gene_id] = genes[gene_id]
            end
        end

        return retval
    end


    """
        get_genome_genes(group, target_genome)

    Returns all genes belonging to a target genome from a group
    """
    function get_genome_genes(group, target_genome)
        retval = Dict{String, Any}()

        for gene_id in keys(group)
            if target_genome in group[gene_id]["genome"]
                retval[gene_id] = group[gene_id]
            end
        end
        return retval
    end


    """
        save_excel(genes)

    Saves genes of interest, control group, and the comparison group into an excel file with each having their own sheet
    """
    function save_excel(genes)
        output_file = Printf.@sprintf("%s_%s.xlsx", load_group(joinpath("src","output","seqco")), Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS"))
        sorted_gene_ids = sort([[gene_id, genes[gene_id]["group_identity"]] for gene_id in collect(keys(genes))], by = x -> x[2], rev=true)
        XLSX.openxlsx(output_file, mode="w") do xf
            # Create the first sheet with the genes of interest
            sheet = xf[1]
            sheet["A1"] = ["Protein ID", "Protein", "Source Genomes", "Group Identity", "Similar Genes"]
            XLSX.rename!(sheet, "new_sheet")
            for (gene_index, gene_id) in enumerate(sorted_gene_ids)
                gene = genes[gene_id[1]]
                similar_gene_str = join([Printf.@sprintf("%s (%f)",similar_gene["protein_id"], similar_gene["similarity"]) for similar_gene in gene["similar_genes"]], " | ")
                row = [gene["protein_id"],gene["protein"], join(gene["genome"], ", "), gene["group_identity"], similar_gene_str]
                sheet[Printf.@sprintf("A%s", gene_index+1)] = row
            end
        end
    end


    """
        main(config_path::String)

    Compare two groups of seqeunces and return the individual elements that are common between both groups.
    """
    function main(config_path::String)
        app_properties = load_properties(config_path)

        # Load the control and comparison groups according to the configs
        if app_properties["control_group"] == "sensitive"
            control_group, control_genomes  = load_group(joinpath("src","input","sensitive"))
            comparison_group, comparison_genomes = load_group(joinpath("src","input","nonsensitive"))
        else
            comparison_group, comparison_genomes = load_group(joinpath("src","input","sensitive"))
            control_group, control_genomes = load_group(joinpath("src","input","nonsensitive"))
        end

        # Determine the elements that are common between all genomes in a group and between the two groups
        common_elements_control = get_common_group_elements(control_group, control_genomes, similarity_metric=app_properties["similarity-metric"]["metric"], acceptance_threshold=app_properties["similarity-metric"]["threshold"])
        common_elements_comparison = get_common_group_elements(comparison_group, comparison_genomes, similarity_metric=app_properties["similarity-metric"]["metric"], acceptance_threshold=app_properties["similarity-metric"]["threshold"])
        common_elements = get_common_seq_elements(common_elements_control, common_elements_comparison)

        # Get all the elements in the control group that is not present in the common elements between the two groups
        elements_of_interest = Dict{String, Any}()
        common_gene_ids = collect(keys(common_elements))
        for control_gene_id in keys(control_group)
            if !(control_gene_id in common_gene_ids)
                elements_of_interest[control_gene_id] = control_group[control_gene_id]
            end
        end
            
        # Calculate group identities for each element in the elements_of_interest
        output = calculate_gene_group_identities!(elements_of_interest, control_group, control_genomes, similarity_metric=app_properties["similarity-metric"]["metric"], acceptance_threshold=app_properties["similarity-metric"]["threshold"], group_identity_threshold=app_properties["group_identity_threshold"])

        # Write output to excel file
        save_excel(output)

    end

    # main(joinpath("src","config","local-config.json"))

end

