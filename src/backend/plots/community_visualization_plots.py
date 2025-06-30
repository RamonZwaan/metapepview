import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# import dash_bio

from typing import List
from copy import deepcopy

import pandas as pd
import numpy as np


from backend.types import MetaPepDbSearch, MetaPepDeNovo
from backend.type_operations import assert_db_search_compatibility, assert_de_novo_compatibility
from backend.post_processing import reference_score_distribution, reference_score_dist_peaks
from backend.spectral_ref_builder import *
from backend.utils.graph_utils import *
from backend.utils.pd_utils import *
from backend.io.import_spectra import *
from constants import *

# assign template theme for plotly figures
pio.templates.default = GraphConstants.default_template


# taxonomy abundance composition
def taxonomic_abundance_barplot(peptide_dataset: pd.DataFrame,
                                rank: RankType='Phylum',
                                abundance_metric: AbundanceMetricType='Match Count',
                                fetch_names: bool=True,
                                include_undefined: bool=False,
                                denovo_threshold: float=90,
                                fractional_abundance: bool=False,
                                topn: int=0):
    # Set column names to process for the figure
    xcol = 'Sample Name'
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'

    # configure taxonomy columns to process based on rank and selected display format
    (display_suffix, hidden_suffix) = (' Name', ' Id') if fetch_names is True \
        else (' Id', ' Name')
    rank_display_col = rank + display_suffix    
    rank_hidden_col = rank + hidden_suffix
    
    # generate lookup table to match display format (id/name) to its corresponding hidden format
    rank_cols = [rank + sfx for sfx in (display_suffix, hidden_suffix)]
    id_name_match = {displ_name: hidden_name for displ_name, hidden_name \
        in peptide_dataset[rank_cols].drop_duplicates().dropna().values}

    # set correct discrete color scale depending on number of taxa
    color_scale = deepcopy(GraphConstants.wide_color_palette) if topn > 10 \
        else deepcopy(GraphConstants.color_palette)
    
    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_display_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0
    other_val = 'Other'

    # convert dataset to rank name and their aggregated sum
    comp_df = peptide_dataset[[ycol, rank_display_col, xcol]]
    
    # give value to all nan cells to represent undefined or unknown
    if include_undefined is True:
        comp_df[rank_display_col].fillna(undefined_val, inplace=True)
            
    # group all peptides belonging to rank, sum separate abundance scores
    comp_df = comp_df.groupby(by=[xcol, rank_display_col])[[ycol]].sum().reset_index()
    
    # divide psm's by sum of psm's per sample name
    if fractional_abundance is True:
        comp_df.loc[:, ycol] /= comp_df.groupby(by=[xcol])[ycol].transform('sum')
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        sorted_taxa_counts = comp_df[comp_df[rank_display_col] != undefined_val]\
            .groupby(by=rank_display_col)[ycol]\
            .sum()\
            .sort_values(ascending=False)
        
        top_taxa = sorted_taxa_counts.head(topn)
        top_taxa = top_taxa.index.to_list()
        
        if include_undefined == True:
            top_taxa.append(undefined_val)
        

        # summ smaller taxa into 'other' group
        top_taxa_indices = comp_df[rank_display_col].isin(top_taxa + [undefined_val])
        other_sum = comp_df[~top_taxa_indices].groupby(by=xcol)[ycol].sum()
        other_df = other_sum.reset_index()
        other_df[rank_display_col] = other_val
        
        comp_df = comp_df[comp_df[rank_display_col].isin(top_taxa)]
        comp_df = comp_df.reset_index(drop=True)
        comp_df = comp_df.sort_values(by=rank_display_col)
        comp_df[rank_display_col] = comp_df[rank_display_col].astype(str)
        
        comp_df = pd.concat([comp_df, other_df]) 
        
    # wether taxonomy id or name is displayed, include the other as hidden data
    comp_df.loc[:, rank_hidden_col] = comp_df[rank_display_col].apply(
        lambda x: id_name_match.get(x, np.nan)
    )
    
    # order categories to have 'Other' and 'Undefined' at the end
    for custom_cat in ['Other', 'Undefined']:
        if custom_cat in comp_df[rank_display_col].values:
            other_df = comp_df[comp_df[rank_display_col] != custom_cat]
            custom_cat_df = comp_df[comp_df[rank_display_col] == custom_cat]
            comp_df = pd.concat([other_df, custom_cat_df])
        
    # assign clear distinct color for 'Undefined group
    if 'Undefined' in comp_df[rank_display_col].values:
        ncats = len(comp_df[rank_display_col].unique())
        # if too more colors than scale + undef color, do nothing
        if ncats <= len(color_scale):
            color_scale[ncats - 1] = GraphConstants.undefined_color
        elif ncats == len(color_scale) + 1:
            color_scale.append(GraphConstants.undefined_color)
            

    fig = px.bar(comp_df,
                 title=f"Distribution PSM over taxa ({rank})",
                 x=xcol,
                 y=ycol,
                 color=rank_display_col,
                 custom_data=rank_hidden_col,
                 color_discrete_sequence=color_scale)
    
    fig.update_layout(GraphConstants.default_layout)
    fig.update_layout(title="")
     
    # configure y-axis title
    ytitle = 'Peptide spectrum matches' if abundance_metric == 'Match Count' else 'Area'
    if fractional_abundance is True:
        ytitle = "Fraction " + ytitle.lower()
    
    fig.update_xaxes(showline=True, 
                     linecolor="Black", 
                     categoryorder="category ascending")
    
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     nticks=5,
                     title=ytitle)
    
    # fig.update_layout(width=900, height=500)
    return fig


# taxonomy abundance composition
def taxonomic_abundance_heatmap(peptide_dataset: pd.DataFrame,
                                rank: RankType="Phylum",
                                abundance_metric: AbundanceMetricType="Match Count",
                                fetch_names: bool=True,
                                include_undefined: bool =False,
                                include_denovo_only: bool=False,
                                denovo_threshold: float=90,
                                fractional_abundance: bool=False,
                                topn: int=0):
    # column name for abundances
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'
    xcol = 'Sample Name'
    # get correct suffix from lineage columns
    suffix = ' Name' if fetch_names is True else ' Id'
    rank_col = rank + suffix            
    
    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0
    other_val = 'Other'
    tax_comp_list = list()
    
    # get composition data for each sample
    for i, sample_name in enumerate(peptide_dataset[xcol].unique()):
        dataset = peptide_dataset[peptide_dataset[xcol] == sample_name]
            
        # convert dataset to rank name and their aggregated sum
        rank_composition = dataset[[ycol, rank_col]]
        
        # give value to all nan cells to represent undefined or unknown
        if include_undefined is True:
            rank_composition[rank_col].fillna(undefined_val, inplace=True)
                
        # group all peptides belonging to rank, sum separate abundance scores
        rank_composition = rank_composition.groupby(by=rank_col).sum().reset_index()
        
        if fractional_abundance is True:
            rank_composition[ycol] /= rank_composition[ycol].sum()
        
        rank_composition.loc[:, xcol] = sample_name
        
        tax_comp_list.append(rank_composition)
        
    # combine samples into one datast
    comp_df = pd.concat(tax_comp_list)
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = comp_df[comp_df[rank + suffix] != undefined_val]\
            .groupby(by=rank + suffix)\
            .sum()[ycol]\
            .sort_values(ascending=False)\
            .head(topn)
        top_taxa = top_taxa.index.to_list()
        
        if include_undefined == True:
            top_taxa.append(undefined_val)
        
        comp_df = comp_df[comp_df[rank_col].isin(top_taxa)]
        comp_df[rank_col] = comp_df[rank_col].astype(str)
        
    # if include_denovo_only is True and denovo_dataset is not None:
    #     de_novo_psm = denovo_dataset[denovo_dataset["ALC (%)"] > denovo_threshold]["peptide_spectrum_matches"].sum()
    #     comp_df.loc[len(comp_df.index)] = [f"de novo only (>{denovo_threshold})",
    #                                        de_novo_psm,
    #                                        denovo_dataset_name]
    
    value_matrix = comp_df[[ycol, xcol, rank_col]].pivot(index=xcol,
                                                         columns=rank_col,
                                                         values=ycol)
    
    # create simple annotated heatmap
    fig = px.imshow(value_matrix,
                    color_continuous_scale=GraphConstants.continuous_color_scale)
    
    # configure x-axis title
    xtitle = 'Peptide spectrum matches' if abundance_metric == 'Match Count' else 'Area'
    if fractional_abundance is True:
        xtitle = "Fraction " + xtitle.lower()
    
    fig.update_xaxes(side="top",
                     title=xtitle)
    fig.update_layout(GraphConstants.default_layout)
    
    return fig


def taxonomy_dropoff_scatter(peptide_df: pd.DataFrame,
                             lineage_counts: Tuple[str, int | float],
                             lineage_dropoff: List[Tuple[str, Tuple[int | float]]],
                             normalize_bars: bool = False):
    
    # fetch allocation categories into separate arrays
    lin_names = [x[0] for x in lineage_dropoff]
    annotation_dropoff = np.array([x[0] for i, x in lineage_dropoff], dtype=np.float64)
    branching_dropoff = np.array([x[1] for i, x in lineage_dropoff], dtype=np.float64)
    valid_counts = np.array([x[2] for i, x in lineage_dropoff], dtype=np.float64)

    other_branches_names = []
    other_branches_values =[]
    for i, x in lineage_dropoff:
        # check nan
        if x[4] == x[4]:
            other_branches_names += [list(x[4].keys())]
            other_branches_values += [np.array(list(x[4].values()), dtype=np.float64)]

    
    # normalize peptide annotation allocation to 100%
    sum_array = annotation_dropoff + branching_dropoff + valid_counts
    if normalize_bars is True:
        annotation_dropoff /= sum_array / 100
        branching_dropoff /= sum_array / 100
        valid_counts /= sum_array / 100
        other_branches_values = [rank / (sum_val / 100) for rank, sum_val in\
            zip(other_branches_values, sum_array)]

        secondary_y = True
    else:
        secondary_y = False
    
    
    # create new figure
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # add lineage counts as scatter trace
    fig.add_trace(go.Scatter(x=[str(row[0]) for row in lineage_counts],
                             y=[row[1] for row in lineage_counts],
                             mode='lines+markers',
                             marker=dict(
                                size=8,
                                color="Black"
                             ),
                             line=dict(width=3),
                             connectgaps=True,
                             name="lineage abundance"),
                  secondary_y=secondary_y,
                  )
    # add lineage dropoff as bar trace
    fig.add_trace(go.Bar(x=lin_names,
                         y=valid_counts,
                         marker=dict(
                                color=GraphConstants.color_palette[0],
                                line_width=1
                         ),
                         name="lineage clade")
                  )
    
    x_data, y_data, custom_data = [], [], []
    for c, lin_elem in enumerate(lin_names):
        # stop if end of list
        if len(other_branches_names) > c:
            x_data += [lin_elem] * len(other_branches_names[c])
            y_data += list(other_branches_values[c])
            custom_data += other_branches_names[c]
    fig.add_trace(go.Bar(x=x_data,
                         y=y_data,
                         customdata=custom_data,
                         hovertemplate="(%{customdata}, %{y})",
                         marker=dict(
                            color=GraphConstants.color_palette[1],
                            line_width=1
                            ),
                        name="diverging clades")
    )
    # fig.add_trace(go.Bar(x=lin_names,
    #                      y=branching_dropoff,
    #                      marker=dict(
    #                             color=GraphConstants.color_palette[1]
    #                          ),
    #                      name="diverging clades")
    #              )
    fig.add_trace(go.Bar(x=lin_names,
                         y=annotation_dropoff,
                         marker=dict(
                                color=GraphConstants.color_palette[2],
                                line_width=1
                             ),
                         name="annotation drop")
                  )
    
    fig.update_layout(GraphConstants.default_layout)
    fig.update_layout(barmode='stack',
                      legend=dict(x=1.1))
    
    if normalize_bars is True:
        fig.update_yaxes(title="abundance alloc lower clade (%)",
                         tickmode="sync",
                         gridcolor=GraphConstants.gridcolor,
                         gridwidth=GraphConstants.gridwidth,
                         nticks=6)
    fig.update_yaxes(title="Peptide spectrum matches",
                     secondary_y=secondary_y,
                     gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     nticks=6,
                     rangemode="tozero")
    fig.update_xaxes(showline=True, linecolor="Black")
    
    return fig


# taxonomy abundance composition
def tax_differential_barplot(peptide_dataset: pd.DataFrame,
                             ratio_numerator_col: str,
                             ratio_denominator_col: str,
                             rank: RankType='Phylum',
                             abundance_metric: AbundanceMetricType='Match Count',
                             fetch_names: bool=True,
                             topn: int=0,
                             fractional_abundance_threshold: float = 0.001,
                             show_legend: bool = True):
    # Set column names to process for the figure
    sample_name = 'Sample Name'
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'

    # get correct suffix from lineage columns
    suffix = ' Name' if fetch_names is True else ' Id'
    rank_col = rank + suffix
    
    # set color scale based on topn count
    color_scale = GraphConstants.wide_color_palette if topn > 10 else GraphConstants.color_palette

    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0

    # convert dataset to rank name and their aggregated sum
    comp_df = peptide_dataset[[ycol, rank_col, sample_name]]
            
    # group all peptides belonging to rank, sum separate abundance scores
    comp_df = comp_df.groupby(by=[sample_name, rank_col]).sum().reset_index()
    
    # divide psm's by sum of psm's per sample name
    comp_df.loc[:, ycol] /= comp_df.groupby(by=[sample_name])[ycol].transform('sum')
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        sorted_taxa_counts = comp_df[comp_df[rank_col] != undefined_val]\
            .groupby(by=rank_col)[ycol]\
            .sum()\
            .sort_values(ascending=False)
        
        # show one less, as this is given in barplot to "Other"
        top_taxa = sorted_taxa_counts.head(topn)
        top_taxa = top_taxa.index.to_list()
        
        comp_df = comp_df[comp_df[rank_col].isin(top_taxa)].reset_index(drop=True)
        comp_df[rank_col] = comp_df[rank_col].astype(str)
    
    # for each taxa, compute ratio between samples
    comp_pivot_df = comp_df.pivot(index=rank_col,
                                  columns=sample_name,
                                  values=ycol).fillna(0.0)
    
    # compute ratio of fractions by dividing the numerator sample with the denominator sample
    comp_pivot_df['Ratio'] = comp_pivot_df[ratio_numerator_col] / comp_pivot_df[ratio_denominator_col]
    comp_pivot_df['Size'] = comp_pivot_df[ratio_numerator_col] + comp_pivot_df[ratio_denominator_col]
    comp_pivot_df['yaxis scatter'] = "abundance"
    
    # rescale ratio to (-inf, inf) to (-1, 1) with log10 and tanh function
    comp_pivot_df['Ratio'] = comp_pivot_df['Ratio'].apply(np.log10).apply(np.tanh)
    comp_pivot_df = comp_pivot_df.reset_index()
    
    # build figure
    fig = make_subplots(rows=2,
                        cols=1,
                        row_heights=[0.7, 0.3],
                        shared_xaxes=True,
                        vertical_spacing=0.01)
    
    fig.add_trace(
        go.Bar(x=comp_pivot_df[rank_col],
               y=comp_pivot_df['Ratio'],
               marker_color=color_scale,
               #marker_color=comp_pivot_df[rank_col],
               #color_discrete_sequence=GraphConstants.color_palette
        ),
        row=1,
        col=1
    )
    
    # take root to let marker area correspond with fraction, multiply to scale to plot
    marker_size_values = (comp_pivot_df['Size'] ** 0.5) * 30
    
    fig.add_trace(
        go.Scatter(x=comp_pivot_df[rank_col],
                   y=comp_pivot_df['yaxis scatter'],
                   mode='markers',
                   marker=dict(
                       size=marker_size_values,
                       color=color_scale,
                   )
        ),
        row=2,
        col=1
    ) 
    
    fig.update_layout(GraphConstants.default_layout)
    fig.update_layout(showlegend=show_legend,
                      margin=dict(l=20, r=20, t=35, b=20),
                      title="Over/under representation de novo annotation over db matching")
    fig.update_traces(width=0.5, row=1, col=1)
    fig.update_layout(title=f"Abundance ratio {ratio_numerator_col}/{ratio_denominator_col}")
    
    fig.update_xaxes(title=None)
    
    fig.update_layout(
        yaxis=dict(gridcolor=GraphConstants.gridcolor,
                   gridwidth=GraphConstants.gridwidth,
                   nticks=5,
                   zeroline=True, 
                   zerolinecolor="Black",
                   showline=True, linecolor="Black", 
                   title='normalized log ratio'),
        yaxis2=dict(gridcolor='rgba(0,0,0,0)',
                    gridwidth=0),
        xaxis2=dict(
            gridcolor='rgba(0,0,0,0)',
            tickvals=comp_pivot_df[rank_col],
            ticktext=comp_pivot_df[rank_col].apply(lambda x: truncate_end(x, 25)),
            gridwidth=0,
            tickangle=45)
    )

    return fig


def pathway_abundance_barplot(dataset: pd.DataFrame, # cols: Protein Name, (Taxonomy Id,) PSM Matches
                              xcol: str,
                              ycol: str,
                              prot_col: str,
                              tax_col: str | None=None,
                              filter_taxa: List[str] | str | None=None,
                              custom_title: str | None=None,
                              custom_xname: str | None=None,
                              custom_yname: str | None=None):
    for colname in [xcol, ycol, prot_col]:
        if colname not in dataset.columns:
            raise ValueError(f"Column {colname} not in dataset.")
        
    if tax_col is not None and tax_col not in dataset.columns:
            raise ValueError(f"Column {tax_col} not in dataset.")
    
    # if taxonomy column represented as numerical, convert to string to prevent color scale
    if tax_col is not None and dataset[tax_col].dtype != "O":
        dataset[tax_col] = dataset[tax_col].astype(str)

    # only show functional data for given taxa
    if filter_taxa is not None and tax_col is not None:
        # make sure format of parameter is list of strings
        if isinstance(filter_taxa, str):
            filter_taxa = [filter_taxa]
        
        # remove unrelated taxa from filter list
        dataset = dataset[dataset[tax_col].isin(filter_taxa)]

    # format title of plot
    if custom_title is None:
        title = "Abundance metabolic enzymes"
    else:
        title = custom_title

    fig = px.bar(dataset,
                 title=title,
                 x=prot_col,
                 y=ycol,
                 hover_name=tax_col,
                 color=xcol,
                 color_discrete_sequence=GraphConstants.color_palette,
                 barmode='group')
    
    fig.update_layout(GraphConstants.default_layout,
                      margin=dict(t=30))
    
    fig.update_xaxes(showline=True,
                     linecolor="Black",
                     categoryorder='total descending')
    
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     nticks=5)
    
    if custom_yname is not None:
        fig.update_yaxes(title = custom_yname)
    if custom_xname is not None:
        fig.update_xaxes(title = custom_xname)
        
    return fig


def de_novo_fraction_barplot(dataset,
                             xcol="Sample Name",
                             ycol="peptides",
                             color_col="db search identified",
                             barmode="relative"):
    title = "Fraction of peptides only de novo identified"
    
    fig = px.bar(dataset,
                 title=title,
                 x=xcol,
                 y=ycol,
                 color=color_col,
                 barmode=barmode)
    return fig