# GIDEON
**An Integer Linear Programming framework for discovering compensatory pathways (Between-Pathway Models, BPMs) in genetic interaction networks**

This repository contains code for the preprint  
_“A Novel ILP Framework to Identify Compensatory Pathways in Genetic Interaction Networks with GIDEON”_  
by **Jocelyn J. Garcia, Kevin Yu, Catherine Freudenreich, and Lenore J. Cowen**. The preprint is available here: https://doi.org/10.64898/2026.03.29.715009.

If you use this work, please cite the authors. We welcome questions about the method and its applications.

---

## Included Data and Outputs

We provide:

- The full set of BPMs identified by GIDEON on data from  
  _A global genetic interaction network maps a wiring diagram of cellular function_ by Costanzo et al. One file lists the BPMs, while the other translates the gene names, lists enriched terms for each module, and notes the stronger interactions included in the BPM.
- The interaction network with DI weights and filtering from the Costanzo 2016 data (interaction_network_DI_rsd2.gi)
- The list of nonessential genes in the network (nonessential_genes.txt)

---

## Quick Start

```bash
# Step 1: Reweight network
python DiWeighting.py \
  --genetic-interaction-filepath data.csv \
  --output-file weighted_network.txt

# Step 2: Run ILP for one gene (parallelizable)
python ILP.py \
  --all-genes-list genes.txt \
  --genetic-interaction-file weighted_network.txt \
  --gene-number 42 \
  --output-folder bpm_outputs \
  --thread-count 8

# Step 3: Trim and prune BPMs
python TrimmingPruning.py \
  --all-genes-list genes.txt \
  --genetic-interaction-filepath weighted_network.txt \
  --raw-bpm-directory bpm_outputs \
  --output-file final_bpms.txt
```

---

## Method Overview

The GIDEON pipeline consists of three stages:

1. **Network reweighting and filtering**  
   Converts raw genetic interaction data into a weighted GI network using DI scores.

2. **ILP-based BPM discovery**  
   Runs an ILP centered at each gene to identify candidate compensatory pathways.

3. **Trimming and pruning**  
   Refines BPMs by removing weakly connected genes and enforcing diversity.

**Input requirement:**  
Raw genetic interaction data containing:
- single mutant fitnesses (SMF)
- double mutant fitnesses (DMF)

---

## Package Requirements

Recommended package versions:

```text
gurobipy==12.0.3
networkx==3.4.2
numpy==1.23.5
pandas==2.2.3
scikit-learn==1.7.2
scipy==1.15.3
```

---

## Step-by-Step Usage

### 1. Reweighting the Network (`DiWeighting.py`)

**Required arguments**
- `--genetic-interaction-filepath`  
  CSV with columns: `GeneA`, `GeneB`, `SMFA`, `SMFB`, `DMF`
- `--output-file`  
  Output path for weighted network

**Optional arguments**
- `--rsd-cutoff-value` (default: 2)  
- `--no-filtering`

**Output format**
```text
geneA    geneB    weight
```

---

### 2. Searching for BPMs with the ILP (`ILP.py`)

Runs one ILP per gene (fully parallelizable).

**Required arguments**
- `--all-genes-list`: line-separated gene list  
- `--genetic-interaction-file`: weighted network file  
- `--gene-number`: index of target gene  
- `--output-folder`: directory for BPM outputs  
- `--thread-count`: Gurobi threads (recommended: 8)

**Notes**
- Requires a **Gurobi license** (free for academic use)  
- Outputs one BPM per gene  

---

### 3. Trimming and Pruning BPMs (`TrimmingPruning.py`)

Refines raw BPMs to improve biological relevance and diversity.

**Required arguments**
- `--all-genes-list`  
- `--genetic-interaction-filepath`  
- `--raw-bpm-directory`  
- `--output-file`

**Optional arguments**
- `--interaction-weight-cutoff` (default: 0.015)

**Output**
- Final list of BPMs (may include multiple components per gene)

---

## Notes

- BPMs are indexed by their source gene; indices may repeat after splitting into connected components.  
- The pipeline is designed for **large-scale parallel execution**.  
- The only required data are the raw single and double mutant fitness measurements from a genetic interaction network.
<!--
# GIDEON
An ILP-based method for mining Between-Pathway Models (BPMs)

This repository contains the associated code for the preprint _A Novel ILP Framework to Identify Compensatory Pathways in Genetic Interaction Networks with GIDEON_ by Jocelyn J. Garcia, Kevin Yu, Catherine Freudenreich, and Lenore J. Cowen. We ask that all uses of this method properly cite the authors, and we welcome questions concerning citation or the method as a whole.

## Associated Outputs
We provide the full set of BPMs identified by GIDEON using the data from _A global genetic interaction network maps a wiring diagram of cellular function_ by Costanzo et. al. We also include the underlying DI weighted and filtered genetic interaction network.

## Usage
The full GIDEON pipeline involves the following three steps: reweighting and filtering the genetic interaction network, running an ILP for every gene in the genetic interaction network, and last trimming and pruning the outputted set of BPMs. The only data needed for the GIDEON framework is the raw data of a genetic interaction network, specifically the single and double mutant fitnesses for genes tested in double knockout experiments. 

### Package Requirements
When running the code, we recommend the following versions:

gurobipy: 12.0.3 <br>
networkx: 3.4.2 <br>
numpy: 1.23.5 <br>
pandas: 2.2.3 <br>
scikit-learn: 1.7.2 <br>
scipy: 1.15.3 <br>

### Reweighting the Network (DiWeighting.py)
To reweight the network, the DiWeighting script requires the following arguments: <br>
--genetic-interaction-filepath: the location of a csv file containing the raw genetic interaction data, with the required columns 'GeneA','GeneB','SMFA','SMFB', and 'DMF' representing the single and double mutant fitnesses of genes A and B respectively. <br>
--output-file: the path for the outputted genetic interaction network with the DI weights.

The following optional arguments can also be used:
--rsd-cutoff-value: the distance with respect to residual standard deviation (rsd) of both component genes required for an interaction to be kept. Default is 2.
--no-filtering: a flag that keeps the genetic interaction network from being filtered.

The genetic interaction network returned by the DiWeighting script has a line for every interaction, with each interaction listed as the tab-delimited gene \t gene \t weight

### Searching for BPMs with the ILP (ILP.py)
To search for BPMs, the GIDEON framework runs the ILP once for every gene in the network, centering the resulting BPM on said gene. We note that solving the ILP script requires a license for Gurobi given the size of the ILP, and there are readily available free academic licenses. To ensure easy parallelization, the script has the gene used to build the BPM as an input and outputs that BPM into a directory, allowing many genes to be run simultaneously as computational resources allow. The required arguments are:
--all-genes-list: a simple .txt file with a line-separated list of all genes included in the network
--genetic-interaction-file: the genetic interaction file outputted by the DiWeights script
--gene-number: the number of the gene (based on the order of --all-genes-list) that should be used to search for a BPM
--output-folder: the directory to output a .bpm text file containing the genes of the BPM, with each gene's BPM file named using the number of the gene
--thread-count: the maximum number of threads for Gurobi. We found that 8 threads was sufficient on our network of 4,459 genes and ~156,000 edges

The output of the ILP.py file, once run on all genes in the network, is a directory of the BPMs identified by each gene, numbered by the gene with respect to its order in the all-genes-list argument.

### Trimming and Pruning BPMs (TrimmingPruning.py)
Once the set of raw BPMs have been generated by the ILP, each BPM is trimmed first to remove weakly related genes, pruned to ensure BPM diversity, and then a Final Trim ensures that BPMs form a single connected component (followed by another prune for BPM diversity). The required arguments for BPM trimming/pruning are:
--all-genes-list: a simple .txt file with a line-separated list of all genes included in the network
--genetic-interaction-filepath: the genetic interaction file outputted by the DiWeights script
--raw-bpm-directory: the directory of BPMs identified by the ILP, numbered by the gene with respect to all-genes-list
--output-file: the name of the resulting list of trimmed and pruned BPMs, numbered by their gene

The following optional argument can also be used:
--interaction-weight-cutoff: the interaction weight requirement for each gene to be kept in the BPM. The default is 0.015.

The output of TrimmingPruning.py is a text file containing the final list of all BPMs. Note that the numbering of these BPMs is not necessarily unique because they are numbered by the gene used to build the BPM and some BPMs are split into their connected components.
-->
