import os
import pandas as pd
import networkx as nx
import csv
from argparse import ArgumentParser


def jaccard_index(bpm1, bpm2):
    """  
    Parameters: bpm1 is a set containing all genes in both modules of the first bpm
                bpm2 is a set containing all genes in both modules of the second bpm
    Purpose:    to compute the jaccard similarity between two bpms
    Returns:    a float representing the jaccard index of the two bpms
    """
    return len(bpm1.intersection(bpm2)) / float(len(bpm1.union(bpm2)))


def calc_gene_interaction_weight(gene, bpm):
    """  
    Parameters: gene is the name of the gene used to compute its interaction weight
                bpm is a tuple of two sets that contain the genes in each bpm module
    Purpose:    to compute the gene interaction weight for a given gene, which means
                only using the edges adjacent to that gene, as well as the number
                of negative interactions for that gene across modules
    Returns:    a float representing the gene interaction weight and an integer
                representing the number of edges across
    """
    module_a, module_b = bpm
    negative_edges_across = 0

    all_bpm_genes = set(module_a).union(set(module_b))
    bpm_graph = nx.induced_subgraph(genetic_interaction_graph, all_bpm_genes).copy()

    # Calculate gene's interaction weight
    gene_weight_across = 0
    gene_weight_within = 0

    for edge in bpm_graph.edges():
        first_node, second_node = edge

        if first_node == gene or second_node == gene:
            edge_weight = genetic_interaction_graph[first_node][second_node]["weight"]
            if within_module(edge, bpm):
                gene_weight_within += edge_weight
            else:
                gene_weight_across += edge_weight
                negative_edges_across += 1 if edge_weight < 0 else 0

    gene_iteraction_weight = (gene_weight_within - gene_weight_across) / len(
        all_bpm_genes
    )

    return gene_iteraction_weight, negative_edges_across


def within_module(edge, bpm):
    """  
    Parameters: edge is a tuple of two genes
                bpm is a tuple of two sets that contain the genes in each bpm module
    Purpose:    to determine whether the edge is within bpm modules or between modules
    Returns:    a boolean representing if the edge is within bpm modules
    """
    first_node, second_node = edge
    module_a, module_b = bpm

    if (first_node in module_a and second_node in module_a) or (
        first_node in module_b and second_node in module_b
    ):
        return True
    else:
        return False


# Prune the trimmed according to LocalCut
def sort_bpms_by_interaction_weight(bpms, bpm_nums):
    """  
    Parameters: bpms is a list of tuples of two sets that contain the genes in each bpm module
                bpm_nums is the number corresponding to the gene used to generate the bpm in bpms
    Purpose:    to sort the bpms by their interaction weight, while making sure that bpm_nums
                still corresponds to the bpm number
    Returns:    a list of tuples corresponding to the sorted bpms and a list of integers
                corresponding to the sorted bpm numbers
    """

    bpms_and_interaction_weights = []

    for i in range(len(bpms)):

        bpm = bpms[i]
        module_a, module_b = bpm
        all_bpm_genes = set(module_a).union(set(module_b))
        interaction_weight = 0

        if len(all_bpm_genes) > 0:
            bpm_graph = nx.induced_subgraph(
                genetic_interaction_graph, all_bpm_genes
            ).copy()

            # Calculate interweight
            weight_within = 0
            weight_across = 0

            for edge in bpm_graph.edges():
                first_node, second_node = edge
                edge_weight = genetic_interaction_graph[first_node][second_node][
                    "weight"
                ]

                if within_module(edge, bpm):
                    weight_within += edge_weight
                else:
                    weight_across += edge_weight

            interaction_weight = (weight_within - weight_across) / len(all_bpm_genes)

        bpms_and_interaction_weights.append((bpm, interaction_weight, bpm_nums[i]))

    sorted_bpms_by_iw = sorted(
        bpms_and_interaction_weights, key=lambda x: x[1], reverse=True
    )
    sorted_bpms, sorted_iw, sorted_bpm_numbers = zip(*sorted_bpms_by_iw)

    return list(sorted_bpms), list(sorted_bpm_numbers)


def prune(sorted_bpms, sorted_bpm_nums):
    """  
    Parameters: sorted_bpms is a list of tuples with two sets corresponding to the genes in
                    each bpm modules
                sorted_bpm_nums is a list of gene numbers corresponding to the gene used to
                    construct the bpm in sorted_bpms
    Purpose:    to prune the bpm set by ensuring that all bpms have a maximum jaccard similarity
                of 0.66 with any other bpm in the outputted set
    Returns:    a list of tuples corresponding to the pruned bpms and a list of integers
                corresponding to the pruned bpm numbers
    """
    pruned_bpms = []
    pruned_bpm_nums = []

    for i in range(len(sorted_bpms)):
        bpm = sorted_bpms[i]
        mod_a, mod_b = bpm
        if (len(mod_a) >= 3 and len(mod_b) >= 3) and all(
            jaccard_index(mod_a.union(mod_b), S1.union(S2)) < 0.66
            for S1, S2 in pruned_bpms
        ):
            pruned_bpms.append(bpm)
            pruned_bpm_nums.append(sorted_bpm_nums[i])

    return pruned_bpms, pruned_bpm_nums


def read_bpm_files(directory):
    """  
    Parameters: directory is a directory of bpm files that are named according to the gene number
                    used to generate the bpm
    Purpose:    to read in the list of bpms outputted by the ilp and the gene numbers used to generate
                    those bpms
    Returns:    a list of tuples corresponding to the bpms and a list of integers corresponding to the
                    bpm number, which is the gene used to generate the bpm
    """
    bpm_file_list = os.listdir(directory)
    bpms = []
    bpm_nums = []

    for bpm_file in bpm_file_list:
        with open(os.path.join(directory, bpm_file)) as file:
            if bpm_file[0:3] == "bpm":
                for line in file:
                    if "Module1" in line:
                        module_a = set(
                            line.split("Module1\t")[1].rstrip("\n").split("\t")
                        )
                    else:
                        module_b = set(
                            line.split("Module2\t")[1].rstrip("\n").split("\t")
                        )

                bpms.append((module_a, module_b))

                bpm_num = int(bpm_file[3:-4])
                bpm_nums.append(bpm_num)

    return bpms, bpm_nums


def read_gin_file_and_genes(gin_file_path, genes_file_path):
    """  
    Parameters: gin_file_path is the location of the genetic interaction network file
                genes_file_path is the location of the list of genes in the network
    Purpose:    to construct the genetic interaction network and a list of all genes used
                    to build the bpms
    Returns:    a genetic interaction networkx graph and a dataframe listing all genes used
                    in the network to build the bpms
    """
    gin = pd.read_csv(gin_file_path, sep="\t", header=None)
    gin.columns = ["Gene1", "Gene2", "weight"]
    gin["type"] = gin["weight"] > 0

    gin_graph = nx.from_pandas_edgelist(
        gin, source="Gene1", target="Gene2", edge_attr=["weight", "type"]
    )

    genes_df = pd.read_csv(genes_file_path, header=None, names=["Gene"])
    genes_df.index = range(1, len(genes_df) + 1)

    return gin_graph, genes_df


def trim(bpms, cutoff):
    """  
    Parameters: bpms is a list of tuples with two sets corresponding to the genes in
                    each bpm modules
                cutoff is the minimum gene interaction weight needed to be kept in a bpm
    Purpose:    to iteratively remove all genes with a genetic interaction weight below the
                cutoff or with less than 2 adjacent negative interactions across modules
    Returns:    a list of tuples corresponding to the trimmed bpms
    """

    trimmed_bpms = []

    for i in range(len(bpms)):
        bpm = bpms[i]
        module_a, module_b = bpm

        bpm_is_shrinking = True

        while bpm_is_shrinking:
            remove_from_module_a = set()
            remove_from_module_b = set()

            for gene in module_a:
                gene_iw, negative_edges_across = calc_gene_interaction_weight(gene, bpm)
                if gene_iw < cutoff or negative_edges_across < 2:
                    remove_from_module_a.add(gene)
            for gene in module_b:
                gene_iw, negative_edges_across = calc_gene_interaction_weight(gene, bpm)
                if gene_iw < cutoff or negative_edges_across < 2:
                    remove_from_module_b.add(gene)

            if len(remove_from_module_a) == 0 and len(remove_from_module_b) == 0:
                bpm_is_shrinking = False
            else:
                module_a = module_a.difference(remove_from_module_a)
                module_b = module_b.difference(remove_from_module_b)
                bpm = (module_a, module_b)

        trimmed_bpms.append((module_a, module_b))

        if i % 100 == 0:
            print("\t", i, "modules trimmed")

    return trimmed_bpms
    



def split_components(unsplit_bpms, unsplit_bpm_nums):
    """  
    Parameters: unsplit_bpms is a list of tuples with two sets corresponding to the genes in
                    each bpm modules
                unsplit_bpm_nums is a list of gene numbers corresponding to the gene used to
                    construct the bpm in unsplit_bpms
    Purpose:    to separate bpms into their connected components after removing all edges with 
                    magnitude less than 0.2, while only adding a split bpm to the final list if
                    both modules have at least three genes
    Returns:    a list of tuples corresponding to the split bpms and a list of integers
                corresponding to the split bpm numbers
    """
    split_bpms = []
    split_bpm_nums = []
    
    for k in range(len(unsplit_bpms)):
        unsplit_bpm = unsplit_bpms[k]
        unsplit_bpm_num = unsplit_bpm_nums[k]
        
        both_mod_genes = unsplit_bpm[0].union(unsplit_bpm[1])
        bpm_graph = genetic_interaction_graph.subgraph(both_mod_genes).copy()
        
        for u,v,data in list(bpm_graph.edges(data = True)):
            if abs(data['weight']) < 0.2:
                bpm_graph.remove_edge(u,v)
                
        for component in nx.connected_components(bpm_graph):
            comp_a = unsplit_bpm[0].intersection(set(component))
            comp_b = unsplit_bpm[1].intersection(set(component))
            
            if len(comp_a) >= 3 and len(comp_b) >= 3:
                split_bpms.append((comp_a,comp_b))
                split_bpm_nums.append(unsplit_bpm_num)

    return split_bpms, split_bpm_nums



def final_trim(first_trimmed_bpms, first_trimmed_bpm_nums):
    """  
    Parameters: first_trimmed_bpms is a list of tuples with two sets corresponding to the genes in
                    each bpm modules, and these bpms should have undergone the first trim already
                first_trimmed_bpm_nums is a list of gene numbers corresponding to the gene used to
                    construct the bpm in first_trimmed_bpms
    Purpose:    to split the bpms into their connected components and then prune that split set,
                ensuring the same measure of bpm diversity
    Returns:    a list of tuples corresponding to the pruned split bpms and a list of integers
                corresponding to the pruned split bpm numbers
    """
    
    split_unpruned_bpms, split_unpruned_bpm_nums = split_components(first_trimmed_bpms, first_trimmed_bpm_nums)
    sorted_split_unpruned_bpms, sorted_split_unpruned_bpm_nums = sort_bpms_by_interaction_weight(split_unpruned_bpms, split_unpruned_bpm_nums)
    pruned_split_bpms, pruned_split_bpm_nums = prune(sorted_split_unpruned_bpms, sorted_split_unpruned_bpm_nums)
    
    return pruned_split_bpms, pruned_split_bpm_nums
                
    
    





def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument(
        "-gf",
        "--all-genes-file",
        required=True,
        help="the line-separated list of all genes in the interaction network",
        type=str,
    )
    parser.add_argument(
        "-gifp",
        "--genetic-interaction-filepath",
        required=True,
        help="the genetic interaction network file",
        type=str,
    )
    parser.add_argument(
        "-iwc",
        "--interaction-weight-cutoff",
        required=False,
        help="the interaction weight requirement for each gene",
        type=float,
        default=0.015,
    )
    parser.add_argument(
        "-bpmd",
        "--raw-bpm-directory",
        required=True,
        help="the folder containing the raw BPMs (no slash!)",
        type=str,
    )
    parser.add_argument(
        "-outf",
        "--output-file",
        required=True,
        help="the location to save the BPM fle",
        type=str,
    )
    return parser



if __name__ == "__main__":

    args = parse_arguments().parse_args()

    print("About to read in raw bpm files...", flush=True)
    raw_bpms, raw_bpm_nums = read_bpm_files(args.raw_bpm_directory + "/")

    print("About to read in the gi file...", flush=True)
    genetic_interaction_graph, all_genes = read_gin_file_and_genes(
        args.genetic_interaction_filepath, args.all_genes_file
    )

    print("About to trim...", flush=True)
    bpms_trimmed = trim(raw_bpms, float(args.interaction_weight_cutoff))
    print("BPMs are trimmed...", flush=True)
    bpms_sorted, bpm_nums_sorted = sort_bpms_by_interaction_weight(
        bpms_trimmed, raw_bpm_nums
    )
    print("BPMs are sorted...", flush=True)
    bpms_pruned, bpms_num_pruned = prune(bpms_sorted, bpm_nums_sorted)
    print("BPMs are pruned...", flush=True)
    
    bpms_final_trim, bpms_num_final_trim = final_trim(bpms_pruned, bpms_num_pruned)
    
    # Output the trimmed & pruned bpms to a single file ------------------------
    trimmed_pruned_bpm_file = args.output_file
    
    with open(trimmed_pruned_bpm_file, "w") as bpmFile:
        out = csv.writer(bpmFile, delimiter="\t")
        for k in range(len(bpms_final_trim)):
            mod = bpms_final_trim[k]
            mod1 = list(mod[0])
            mod2 = list(mod[1])
    
            bpm_number = bpms_num_final_trim[k]
    
            out.writerow(["BPM%d/Module1" % bpm_number] + mod1)
            out.writerow(["BPM%d/Module2" % bpm_number] + mod2)


