# GIDEON
An ILP-based method for mining Between-Pathway Models (BPMs)

This repository contains the associated code for the preprint _A Novel ILP Framework to Identify Compensatory Pathways in Genetic Interaction Networks with GIDEON_ by Jocelyn J. Garcia, Kevin Yu, Catherine Freudenreich, and Lenore J. Cowen. We ask that all uses of this method properly cite the authors.

## Usage
The full GIDEON framework involves the following three steps: reweighting and filtering the genetic interaction network, running an ILP for every gene in the genetic interaction network, and last trimming and pruning the outputted set of BPMs. For packages, we recommend the following versions:

gurobipy: 12.0.3
networkx: 3.4.2
numpy: 1.23.5
pandas: 2.2.3
scikit-learn: 1.7.2
scipy: 1.15.3
