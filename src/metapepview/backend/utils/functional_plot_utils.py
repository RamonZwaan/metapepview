from typing import List

import pandas as pd


def calculate_frac_abundance(peptide_df: pd.DataFrame,
                             abundance_col: str) -> pd.DataFrame:
    """Convert absolute abundances to fractional abundance per sample name.

    Args:
        peptide_df (pd.DataFrame): MetaPepTable peptide dataset.
        abundance_col (str): Abundance column in dataset.

    Returns:
        pd.DataFrame: MetaPepTable peptide dataset with fractional abundances.
    """
    for sample in peptide_df["Sample Name"].unique():
        sample_idx = peptide_df[peptide_df["Sample Name"] == sample].index.to_list()
        peptide_df.loc[sample_idx, abundance_col] /= peptide_df.loc[sample_idx, abundance_col].sum() # type:ignore

    return peptide_df


def configure_plot_title(kegg_group_format: str,
                         filter_clade: bool,
                         clade_rank: str,
                         func_abundances: pd.DataFrame):
    # title prefix
    title = "functional abundance"

    # add clade filter to title
    if filter_clade and clade_rank and clade_rank != 'Root':
        title += f" from {filter_clade} clade"

    # filter dataset with annotations belonging to nitrogen cycle
    if kegg_group_format == "Protein/Gene name":
        title = "Protein " + title
    elif kegg_group_format == "EC":
        title = "Enzyme " + title
    elif kegg_group_format == "Module":
        title = "KEGG module " + title
    else:
        title = "KEGG orthology " + title

    # check if elements exceeds limit, set in plot title
    if func_abundances.shape[0] > 20:
        title += " (top 20 abundant)"
    
    return title


def get_top_functions(peptide_df: pd.DataFrame,
                      func_abundances: pd.DataFrame):
    func_abundance_sorted = func_abundances.sort_values(ascending=False)

    # check if elements exceeds limit, set in plot title
    if func_abundance_sorted.shape[0] > 20:
        top_n_func = func_abundance_sorted.head(20)
    else:
        top_n_func = func_abundance_sorted
    
    peptide_df = peptide_df[peptide_df["Protein Name"].isin(set(top_n_func.index))]

    return peptide_df


def process_kegg_groups(
        peptide_df: pd.DataFrame,
        kegg_db: "KeggDatabase",
        ycol: str,
        include_taxa: bool,
        kegg_group_method: str,
        kegg_group_format: str,
        custom_prot: List[str] | None,
        predifined_pathway: str | None = None,
        predifined_module: str | None = None,
        predifined_brite_group: str | None = None,
        combine_annotations: bool = True,

):
    # define column lists required during processing when taxonomy is included and when not 
    if include_taxa is True:
        cols_filter_df = [ycol, "Sample Name", "Taxonomy Name", "KEGG_ko"]              # use to filter unrelevant columns
        cols_group_peptides = ["Peptide Index", ycol, "Taxonomy Name", "Sample Name"]   # use to merge protein names by peptide sequence
        cols_group_names = ["Protein Name", "Taxonomy Name", "Sample Name"]                                   # use to group and count psm's by function and taxonomy
    else:
        cols_filter_df = [ycol, "Sample Name", "KEGG_ko"]                             # use to filter unrelevant columns
        cols_group_peptides = ["Peptide Index", ycol, "Sample Name"]                  # use to merge protein names by peptide sequence
        cols_group_names = ["Protein Name", "Sample Name"]                                                  # use to group and count psm's by function
        
    # filter unrelevant columns out of the dataset    
    peptide_df = peptide_df[cols_filter_df]
    
    # split cells with multiple kegg annotations into separate rows
    peptide_df.loc[:, "KEGG_ko"] = peptide_df["KEGG_ko"].str.split(",") 

    # ensure unique index prior to explode, so that duplicate id's correspond to same original row
    peptide_df = peptide_df.reset_index(drop=True) 
    # put peptides annotated towards multiple kegg ko's in separate rows
    peptide_df = peptide_df.explode("KEGG_ko")
    
    # filter peptide dataset outside of predifined pathways/modules/brite
    if kegg_group_method == "Pathway":
        valid_ko = kegg_db.list_ko(pathway=predifined_pathway,
                                   module=predifined_module)
    elif kegg_group_method == "BRITE":
        valid_ko = kegg_db.list_ko(brite=predifined_brite_group)
    else:
        valid_ko = custom_prot
    peptide_df = peptide_df[peptide_df["KEGG_ko"].isin(valid_ko)]
    
    # filter dataset with annotations belonging to nitrogen cycle
    if kegg_group_format == "Protein/Gene name":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_symbol)
    elif kegg_group_format == "EC":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_ec)
    elif kegg_group_format == "Module":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_module)
    else:
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"]
        
    # one ko may correspond to multiple modules or ec, explode to separate rows
    peptide_df = peptide_df.explode("Protein Name")
    # ko under pathway may also match towards modules outside of pathway, filter these
    if kegg_group_format == "Module" and predifined_pathway is not None:
        peptide_df[peptide_df["Protein Name"].isin(kegg_db.list_modules(predifined_pathway))]
    
    peptide_df.dropna(subset="Protein Name", inplace=True)
    
    # drop kegg ko column
    peptide_df.drop("KEGG_ko", axis=1, inplace=True)
    
    # merge duplicate peptides whose protein name are the same
    peptide_df["Peptide Index"] = peptide_df.index
    peptide_df.drop_duplicates(subset=["Peptide Index", "Protein Name"], inplace=True)
    
    # for peptide duplicates with different protein names there are two options:
    #   Count psm towards both names: this does result in inflated PSM numbers
    #   Combine names (alphabetically): Most realistic data, but results in more categories that may not show in the plot.
    if combine_annotations is True:
        peptide_df = peptide_df.groupby(by=cols_group_peptides)["Protein Name"]\
            .agg(lambda x: ",".join(sorted(x)))
        peptide_df = peptide_df.to_frame().reset_index(names=cols_group_peptides)
        
    # with combined protein name categories, group by these categories to obtain psm counts
    peptide_df = peptide_df.groupby(by=cols_group_names)[ycol].agg('sum')
    peptide_df = peptide_df.to_frame().reset_index(names=cols_group_names)

    return peptide_df

