from argparse import ArgumentParser
from gurobipy import Model, GRB, quicksum
import csv


def run_ilp(nodes, edges, gene, num_threads):
    """  
    Parameters: nodes is list of all nodes in the graphs, which are numbers that can be mapped
                    back to the node name
                edges is the list of tuples for all interactions that include the two component
                    genes and their interaction weight
                gene is the number of the gene being centered on to build the bpm
                num_threads is the maximum number of threads that Gurobi should use
    Purpose:    to construct a bpm centered around a given gene using an ILP
    Returns:    a tuple containing two lists of genes that correspond to the two bpm modules
    """
    model = Model("mip")

    # Store the model variables, as well as the weights
    weights_dict = dict()
    right_node_vars = dict()
    left_node_vars = dict()
    outside_node_vars = dict()
    same_edge_vars = dict()
    opposite_edge_vars = dict()

    # Form adjacency list for efficient iteration of adjacent edges to a node
    adjacency_list = dict()

    # Add node variables and constraints for individual nodes
    for n in nodes:

        # Make left, right, and outside binary variables for each node
        left_node_vars[n] = model.addVar(vtype=GRB.BINARY, name="l_%d" % n)
        right_node_vars[n] = model.addVar(vtype=GRB.BINARY, name="r_%d" % n)
        outside_node_vars[n] = model.addVar(vtype=GRB.BINARY, name="o_%d" % n)

        # Fix the node of interest to be in the left module
        if n == gene:
            left_node_vars[n].lb = 1
            left_node_vars[n].ub = 1
            right_node_vars[n].lb = 0
            right_node_vars[n].ub = 0
            outside_node_vars[n].lb = 0
            outside_node_vars[n].ub = 0

        # Ensure that the node is either in the left module, in the right module, or "outside" of the BPM
        model.addConstr(
            left_node_vars[n] + right_node_vars[n] + outside_node_vars[n] == 1,
            name="node_%d_constr" % n,
        )

        # Add all genes to the adjacency list
        adjacency_list[n] = []

    #########################
    ###### Constraints ######
    #########################

    # Add constraint for the sizes of each gBPM
    model.addConstr(quicksum(left_node_vars[i] for i in left_node_vars) == 25, name = "size_left_constr")
    model.addConstr(quicksum(right_node_vars[i] for i in right_node_vars) == 25, name = "size_right_constr")

    # Add edges and constraints for individual edges
    for edge in edges:
        first_edge_node, second_edge_node, edge_weight = edge

        # Make edge variables
        same_edge_vars[first_edge_node, second_edge_node] = model.addVar(
            vtype=GRB.BINARY, name="same_%d_%d" % (first_edge_node, second_edge_node)
        )
        opposite_edge_vars[first_edge_node, second_edge_node] = model.addVar(
            vtype=GRB.BINARY, name="opp_%d_%d" % (first_edge_node, second_edge_node)
        )

        # Add the square of the edge weight (with sign retained) to the dictionary
        if edge_weight < 0:
            weights_dict[first_edge_node, second_edge_node] = (
                float(edge_weight) ** 2
            ) * -1
            weights_dict[second_edge_node, first_edge_node] = (
                float(edge_weight) ** 2
            ) * -1
        else:
            weights_dict[first_edge_node, second_edge_node] = float(edge_weight) ** 2
            weights_dict[second_edge_node, first_edge_node] = float(edge_weight) ** 2

        # Add edge constraints
        model.addConstr(
            opposite_edge_vars[first_edge_node, second_edge_node]
            >= right_node_vars[first_edge_node] + left_node_vars[second_edge_node] - 1,
            name="opp_constr_rf_ls_%d_%d" % (first_edge_node, second_edge_node)
        )
        model.addConstr(
            opposite_edge_vars[first_edge_node, second_edge_node]
            >= left_node_vars[first_edge_node] + right_node_vars[second_edge_node] - 1,
            name="opp_constr_lf_rs_%d_%d" % (first_edge_node, second_edge_node)
        )

        model.addConstr(
            same_edge_vars[first_edge_node, second_edge_node]
            >= right_node_vars[first_edge_node] + right_node_vars[second_edge_node] - 1,
            name="same_constr_rf_rs_%d_%d" % (first_edge_node, second_edge_node)
        )
        model.addConstr(
            same_edge_vars[first_edge_node, second_edge_node]
            >= left_node_vars[first_edge_node] + left_node_vars[second_edge_node] - 1,
            name="same_constr_lf_ls_%d_%d" % (first_edge_node, second_edge_node)
        )

        model.addConstr(
            same_edge_vars[first_edge_node, second_edge_node]
            + opposite_edge_vars[first_edge_node, second_edge_node]
            <= 1 - outside_node_vars[first_edge_node],
            name="out_constr_of_%d_%d" % (first_edge_node, second_edge_node)
        )
        model.addConstr(
            same_edge_vars[first_edge_node, second_edge_node]
            + opposite_edge_vars[first_edge_node, second_edge_node]
            <= 1 - outside_node_vars[second_edge_node],
            name="out_constr_os_%d_%d" % (first_edge_node, second_edge_node)
        )

        # Add edges to adjacency for each node
        adjacency_list[first_edge_node].append((first_edge_node, second_edge_node))
        adjacency_list[second_edge_node].append((first_edge_node, second_edge_node))

    # Gene-centered constraints
    # When n == gene, these constraints are vacuously true
    for n in nodes:
        model.addConstr(
            quicksum(
                weights_dict[i, j] * (same_edge_vars[i, j] - opposite_edge_vars[i, j])
                for i, j in adjacency_list[n]
            )
            <= quicksum(
                weights_dict[i, j] * (same_edge_vars[i, j] - opposite_edge_vars[i, j])
                for i, j in adjacency_list[gene]
            ),
            name="most_iw_%d" % n
        )  # Gene used to build the BPM must have the largest interaction weight
        model.addConstr(
            quicksum(
                weights_dict[i, j] * opposite_edge_vars[i, j]
                for i, j in adjacency_list[gene]
            )
            <= quicksum(
                weights_dict[i, j] * opposite_edge_vars[i, j]
                for i, j in adjacency_list[n]
            ),
            name="mostOpp_%d" % n,
        )  # Gene used to build the BPM must specifically have the most negative opposite weight

    ########################################
    ###### Solve & Interpret Solution ######
    ########################################

    # Instruct ILP to maximize the interaction weight, can iterate over same_edge_vars because every edge has one
    objective = quicksum(
        float(weights_dict[i, j]) * (same_edge_vars[i, j] - opposite_edge_vars[i, j])
        for i, j in same_edge_vars
    )
    model.ModelSense = GRB.MAXIMIZE
    model.setObjective(objective)

    # Set Gurobi to use only the No Relaxation Heuristic, and for 300 work units
    model.setParam("Threads", num_threads)
    model.Params.WorkLimit = 300 # change this blah
    model.Params.NoRelHeurWork = 300 
    # model.Params.MIPFocus = 1 # remove this blah

    # Solve the ILP
    model.optimize()
    
    print()
    print("All variable values:")
    for v in model.getVars():
        print(f"{v.VarName} = {v.X}")
    print()
    
    # feasible = is_solution_feasible(model)
    print(f"Optimization status: {model.Status}")
    print(f"Solution count: {model.SolCount}")

    # Determine which nodes are in the left and right modules
    left_module_nodes = []
    right_module_nodes = []
    
    binary_tolerance = 1e-4

    for n in nodes:
        retrieved_right_var = model.getVarByName("r_%d" % n)
        if retrieved_right_var.X > 1 - binary_tolerance:
            right_module_nodes.append(n)
        retrieved_left_var = model.getVarByName("l_%d" % n)
        if retrieved_left_var.X > 1 - binary_tolerance:
            left_module_nodes.append(n)

    return (left_module_nodes, right_module_nodes)


def make_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-gl",
        "--all-genes-list",
        required=True,
        help="the list of all genes in the interaction network",
        type=str,
    )
    parser.add_argument(
        "-gif",
        "--genetic-interaction-file",
        required=True,
        help="the genetic interaction network file",
        type=str,
    )
    parser.add_argument(
        "-gnum",
        "--gene-number",
        required=True,
        help="the number of the gene to build the BPM with respect to the allGenes file",
        type=int,
    )
    parser.add_argument(
        "-outf",
        "--output-folder",
        required=True,
        help="the folder for the bpm output file (no forward slash!)",
        type=str,
    )
    parser.add_argument(
        "-ths",
        "--thread-count",
        required=True,
        help="number of threads to use",
        type=int,
    )
    return parser


if __name__ == "__main__":
    args = make_parser().parse_args()

    print("This is the gene for use: " + str(args.gene_number))

    # Dictionaries for mapping from gene name to index
    node_to_index_dict = {}
    index_to_node_dict = {}
    curr_node_index = 0

    # Parse gene list
    all_nodes = []
    with open(args.all_genes_list, "r") as f:
        for line in f:
            line = line.split("\t")
            node1 = line[0].strip()
            node_to_index_dict[node1] = curr_node_index
            index_to_node_dict[curr_node_index] = node1
            all_nodes.append(curr_node_index)
            curr_node_index += 1

    # Parse genetic interactions list
    all_edges = []
    with open(args.genetic_interaction_file, "r") as f:
        for line in f:
            line = line.split("\t")
            gene1 = line[0].strip()
            gene2 = line[1].strip()
            weight = float(line[2].strip())

            edge1 = node_to_index_dict[gene1]
            edge2 = node_to_index_dict[gene2]
            all_edges.append((edge1, edge2, weight))

    # Run the ILP
    gene_num = int(args.gene_number)
    left_module_indices, right_module_indices = run_ilp(
        all_nodes, all_edges, gene_num, int(args.thread_count)
    )

    # Translate the indices to the node names
    left_module_translated = []
    right_module_translated = []
    for s in left_module_indices:
        left_module_translated.append(str(index_to_node_dict[s]))
    for t in right_module_indices:
        right_module_translated.append(str(index_to_node_dict[t]))

    print("Bottom Genes:", left_module_translated)
    print("Top Genes:", right_module_translated)

    # Write the BPM to the file
    bpm_file_name = args.output_folder + "/bpm" + str(gene_num) + ".bpm"
    bpm_file = open(bpm_file_name, "w+")
    out = csv.writer(bpm_file, delimiter="\t")
    out.writerow(["BPM%d/Module1" % gene_num] + sorted(list(left_module_translated)))
    out.writerow(["BPM%d/Module2" % gene_num] + sorted(list(right_module_translated)))

