import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import gmean
from argparse import ArgumentParser


def read_data(file_path):
    """  
    Parameters: file_path is the location of the csv file that contains the single and
                double mutant fitnesses, with column names "GeneA", "GeneB", "SMFA", "SMFB",
                and "DMF" that contain gene identifiers, their respective single mutant fitnesses,
                and their double mutant fitness
    Purpose:    to read in the pandas dataframe with only the relevant columns for weighting
    Returns:    pandas data frame
    """
    df_all_columns = pd.read_csv(file_path)
    df = df_all_columns[["GeneA", "GeneB", "SMFA", "SMFB", "DMF"]]
    return df


def find_smf_differences(df):
    """  
    Parameters: df is the data frame containing the single and double mutant fitnesses
    Purpose:    to identify all genes with discrepances in their single mutant fitness based
                on whether that gene was used for the array or query strain
    Returns:    an array containing all genes with different single mutant fitnesses
    """
    genes_with_smf_differences = []

    unique_genes_a_rows = df.drop_duplicates(subset=["GeneA"])
    unique_genes_b_rows = df.drop_duplicates(subset=["GeneB"])

    genes_a_smfs = dict()

    for index, row in unique_genes_a_rows.iterrows():
        genes_a_smfs[row["GeneA"]] = row["SMFA"]

    for index, row in unique_genes_b_rows.iterrows():
        if (
            row["GeneB"] in genes_a_smfs.keys()
            and genes_a_smfs[row["GeneB"]] != row["SMFB"]
        ):
            genes_with_smf_differences.append(row["GeneB"])

    return genes_with_smf_differences


def get_gene_list(df):
    """  
    Parameters: df is the data frame containing the single and double mutant fitnesses
    Purpose:    to obtain all unique genes included in the double knockout
    Returns:    a numpy array containing unique genes in the double knockout data
    """
    unique_genes_a, unique_genes_b = df["GeneA"].unique(), df["GeneB"].unique()

    return np.union1d(unique_genes_a, unique_genes_b)


def predict_using_marginals(df):
    """  
    Parameters: df is the data frame containing the single and double mutant fitnesses
    Purpose:    to compute the predict the double mutant fitness for each double knockout
                using a linear model constructed on each component gene's marginal distribution.
                The marginal is constructed separately based on the gene being in the array or query
                strain when there is a difference in single mutant fitness
    Returns:    a pandas dataframe with the predicted double mutant fitness and the squared residual standard
                deviation for each marginal
    """
    group_gene_a, group_gene_b = df.groupby("GeneA"), df.groupby("GeneB")

    all_gene_list = get_gene_list(df)

    genes_with_smf_differences = find_smf_differences(df)

    for gene in all_gene_list:
        gene_a_rows = (
            group_gene_a.get_group(gene)
            if gene in group_gene_a.groups
            else pd.DataFrame(columns=["SMFB", "DMF"])
        )
        gene_b_rows = (
            group_gene_b.get_group(gene)
            if gene in group_gene_b.groups
            else pd.DataFrame(columns=["SMFA", "DMF"])
        )

        if gene in genes_with_smf_differences:

            if len(gene_a_rows) > 0:
                model = LinearRegression()
                model.fit(gene_a_rows[["SMFB"]], gene_a_rows["DMF"])
                dmf_prediction_a = model.predict(gene_a_rows[["SMFB"]])

                residuals_a = gene_a_rows["DMF"] - dmf_prediction_a

                if len(dmf_prediction_a) > 2:
                    standard_error_a = np.square(residuals_a).sum() / (
                        len(dmf_prediction_a) - 2
                    )
                else:
                    standard_error_a = (
                        np.inf
                    )  # ensures that the edge is kept no matter what when there are few results (but this shouldn't happen on high-throughput data)

                df.loc[gene_a_rows.index, "Predicted (GeneA Marginal)"] = (
                    dmf_prediction_a
                )
                df.loc[gene_a_rows.index, "Variance (GeneA Marginal)"] = (
                    standard_error_a.item()
                )

            if len(gene_b_rows) > 0:
                model = LinearRegression()
                model.fit(gene_b_rows[["SMFA"]], gene_b_rows["DMF"])
                dmf_prediction_b = model.predict(gene_b_rows[["SMFA"]])

                residuals_b = gene_b_rows["DMF"] - dmf_prediction_b
                if len(dmf_prediction_b) > 2:
                    standard_error_b = np.square(residuals_b).sum() / (
                        len(dmf_prediction_b) - 2
                    )
                else:
                    standard_error_b = np.inf

                df.loc[gene_b_rows.index, "Predicted (GeneB Marginal)"] = (
                    dmf_prediction_b
                )
                df.loc[gene_b_rows.index, "Variance (GeneB Marginal)"] = (
                    standard_error_b.item()
                )

        else:
            opposite_gene_fitnesses = (
                pd.concat(
                    [
                        (
                            gene_a_rows["SMFB"].rename("SMF")
                            if not gene_a_rows.empty
                            else None
                        ),
                        (
                            gene_b_rows["SMFA"].rename("SMF")
                            if not gene_b_rows.empty
                            else None
                        ),
                    ]
                )
                .to_frame()
                .reset_index(drop=True)
            )

            double_knockout_fitnesses = (
                pd.concat(
                    [
                        gene_a_rows[["DMF"]] if not gene_a_rows.empty else None,
                        gene_b_rows[["DMF"]] if not gene_b_rows.empty else None,
                    ]
                )
                .reset_index(drop=True)
                .squeeze()
            )

            model = LinearRegression()
            model.fit(opposite_gene_fitnesses, double_knockout_fitnesses)

            predictions = model.predict(opposite_gene_fitnesses)

            residuals_all = predictions - double_knockout_fitnesses

            if len(residuals_all) > 2:
                standard_error_all = np.square(residuals_all).sum() / (
                    len(double_knockout_fitnesses) - 2
                )
            else:
                standard_error_all = np.inf

            gene_a_predictions = predictions[: len(gene_a_rows)]
            gene_b_predictions = predictions[len(gene_a_rows) :]

            df.loc[gene_a_rows.index, "Predicted (GeneA Marginal)"] = gene_a_predictions
            df.loc[gene_b_rows.index, "Predicted (GeneB Marginal)"] = gene_b_predictions

            df.loc[gene_a_rows.index, "Variance (GeneA Marginal)"] = standard_error_all
            df.loc[gene_b_rows.index, "Variance (GeneB Marginal)"] = standard_error_all

    return df


def calculate_weights(df):
    """  
    Parameters: df is the data frame containing the single and double mutant fitnesses
                as well as the predicted double mutant fitnesses and their variances from rsd
    Purpose:    to compute the genetic interaction for each double knockout using the predicted
                double mutant fitness from each marginal distribution, as well as the number of
                standard deviations the double mutant fitness is from each marginal prediction
    Returns:    a pandas dataframe containing the DI weight and the standard deviations from predicted
    """
    df["DI Weight"] = df["DMF"] - gmean(
        df[["Predicted (GeneA Marginal)", "Predicted (GeneB Marginal)"]], axis=1
    )

    df["GeneA RSD"], df["GeneB RSD"] = np.sqrt(
        df["Variance (GeneA Marginal)"]
    ), np.sqrt(df["Variance (GeneB Marginal)"])
    df["SD from PredictA"] = (
        (df["DMF"] - df["Predicted (GeneA Marginal)"]) / df["GeneA RSD"]
    ).abs()
    df["SD from PredictB"] = (
        (df["DMF"] - df["Predicted (GeneB Marginal)"]) / df["GeneB RSD"]
    ).abs()

    df["Min SD from Prediction"] = df[["SD from PredictA", "SD from PredictB"]].min(
        axis=1
    )

    return df


def generate_weight_file(df, cutoff, file_path, no_filtering):
    """  
    Parameters: df is the data frame containing the single and double mutant fitnesses, their
                    interaction scores, and the minimum rsd from predicted
                cutoff is the minimum number of rsds from both predictions that the interaction must be 
                    to be kept in the network
                file_path is the location of the resulting genetic interaction network file
                no_filtering is a boolean for when the network should not be filtered at all using a cutoff
    Purpose:    to filter the data based on the minimum rsd from predicted and output the
                genetic interaction file
    Returns:    n/a
    """

    if not no_filtering:
        df = df[df["Min SD from Prediction"] > cutoff]

    gene_gene_weight = df[["GeneA", "GeneB", "DI Weight"]]
    gene_gene_weight.to_csv(file_path, sep="\t", header=False, index=False)


def argument_parser():
    parser = ArgumentParser()

    parser.add_argument(
        "-gifp",
        "--genetic-interaction-filepath",
        required=True,
        help="the csv file containing the genetic interactions, with the columns 'GeneA','GeneB','SMFA','SMFB', and 'DMF' representing the single and double mutant fitnesses of genes A and B respectively",
        type=str,
    )
    parser.add_argument(
        "-outd",
        "--output-file",
        required=True,
        help="the desired location and name of the genetic interaction file with the computed weights",
        type=str,
    )
    parser.add_argument(
        "-cut",
        "--rsd-cutoff-value",
        required=False,
        help="the distance with respect to residual standard deviation (rsd) of both component genes required for an interaction to be kept. Default is 2.",
        type=float,
        default=2,
    )
    parser.add_argument(
        "-nf",
        "--no-filtering",
        required=False,
        action="store_true",
        help="Whether the genetic interaction network is filtered. Default is False, meaning that the genetic interaction network is filtered based on the rsd cutoff value",
    )
    return parser


## Assuming columns are: GeneA, GeneB, SMFA, SMFB, DMF

if __name__ == "__main__":
    args = argument_parser().parse_args()

    input_file_path = args.genetic_interaction_filepath
    output_file_path = args.output_file
    rsd_cutoff = float(args.rsd_cutoff_value)

    print("Reading data...", flush=True)
    inputData = read_data(input_file_path)

    print("Starting predictions...", flush=True)
    dmf_predictions = predict_using_marginals(inputData)

    print("Calculating and outputting weights...", flush=True)
    weights = calculate_weights(dmf_predictions)
    generate_weight_file(weights, rsd_cutoff, output_file_path, args.no_filtering)
