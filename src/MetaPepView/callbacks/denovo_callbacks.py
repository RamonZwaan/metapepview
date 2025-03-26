from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from MetaPepView.server import app

from backend import *
from backend.plots import taxonomic_abundance_barplot, taxonomic_abundance_heatmap, de_novo_fraction_barplot

import base64
import io
import pandas as pd


@app.callback(
    Output('denovo_fraction_barplot_graph', 'children'),
    Input('peptides', 'data'),
    Input('de_novo_fraction_alc_cutoff', 'value'),
    Input('de_novo_fraction_radio', 'value')
)
def denovo_only_fraction(peptide_json,
                         score_cutoff=80,
                         show_fraction=False):
    """Obtain for each sample the fraction of sequences found through
    de novo processing not identified by db search.

    """
    
    score_cutoff = float(score_cutoff)
    
    # get pandas df from json
    if peptide_json is None:
        return html.P("Import peptide dataset...")
    
    peptide_dataset = pd.read_json(peptide_json)
    
    # filter out unrelevant columns
    peptide_dataset = peptide_dataset[["Sample Name", "ALC (%)", "peptide_spectrum_matches"]]
    
    # score to float
    peptide_dataset["ALC (%)"] = peptide_dataset["ALC (%)"].astype(float)
    
    # assign rows db search identified and non identified values
    peptide_dataset["db search identified"] = ~peptide_dataset["peptide_spectrum_matches"].isna()
    
    
    # filter out de novo only (no db search) rows with score below cutoff
    peptide_dataset = peptide_dataset[(peptide_dataset["ALC (%)"] >= score_cutoff) | (peptide_dataset["db search identified"] == True)]

    # group samples and count db search identified and non identified
    group_df = peptide_dataset.groupby(by=["Sample Name", "db search identified"])
    
    # obtain group sizes
    group_sizes = group_df.size()
    fraction_df = group_sizes.reset_index(name="peptides")

    # if fractions shown, divide db search identified values by sum of groups for each sample
    if show_fraction is True:
        sample_sums = fraction_df.groupby("Sample Name")["peptides"].agg("sum")
        division_col = fraction_df["Sample Name"].apply(lambda x: sample_sums.loc[x])
        
        # divide peptide counts by the sum of the sample name counts
        fraction_df["peptides"] = fraction_df["peptides"].divide(division_col)

    # create the figure
    plot = de_novo_fraction_barplot(fraction_df)
    return dcc.Graph(figure=plot, id="denovo_fraction_figure")
    