import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from typing import Any, Dict, Tuple, List
from copy import deepcopy

import pandas as pd
import numpy as np

from metapepview.constants import GraphConstants
from metapepview.backend.utils.string_utils import truncate_end



def process_tax_abundance(peptide_dataset: pd.DataFrame,
                          abundance_metric: str,
                          tax_col: str,
                          include_undefined: bool,
                          undefined_val: Any,
                          fractional_abundance: bool) -> Tuple[pd.DataFrame, str]:
    # Set column names to process for the figure
    xcol = 'Sample Name'
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'
    
    # convert dataset to rank name and their aggregated sum
    comp_df = peptide_dataset[[ycol, tax_col, xcol]]
    
    # give value to all nan cells to represent undefined or unknown
    if include_undefined is True:
        comp_df.loc[:, tax_col] = comp_df[tax_col].fillna(undefined_val)
            
    # group all peptides belonging to rank, sum separate abundance scores
    comp_df = comp_df.groupby(by=[xcol, tax_col])[[ycol]].sum().reset_index()
    
    # divide psm's by sum of psm's per sample name
    if fractional_abundance is True:
        comp_df.loc[:, ycol] /= comp_df.groupby(by=[xcol])[ycol].transform('sum')

    return (comp_df, ycol)


def determine_top_taxa(
        comp_df: pd.DataFrame,
        facet_df: pd.DataFrame | None,
        topn: int,
        tax_col: str,
        ycol: str,
        facet_ycol: str,
        undefined_val: Any,
        include_undefined: bool) -> List[str]:
    taxa_counts = comp_df[comp_df[tax_col] != undefined_val]\
        .groupby(by=tax_col)[ycol]\
        .sum()
    if facet_df is None:
        # sort counts and get top n
        sorted_taxa_counts = taxa_counts.sort_values(ascending=False)
        top_taxa = sorted_taxa_counts.head(topn)
        top_taxa = top_taxa.index.to_list()
    else:
        facet_counts = facet_df[facet_df[tax_col] != undefined_val]\
            .groupby(by=tax_col)[facet_ycol]\
            .sum()
        # if facet plot present, normalize counts for both and add them
        taxa_counts /= taxa_counts.max()
        facet_counts /= facet_counts.max()

        combined_counts = taxa_counts\
            .add(facet_counts, fill_value=0)\
            .sort_values(ascending=False)

        top_taxa = combined_counts.head(topn)
        top_taxa = top_taxa.index.to_list()
        
    if include_undefined == True:
        top_taxa.append(undefined_val)
    
    return top_taxa

def add_other_group(
        comp_df: pd.DataFrame,
        tax_col: str,
        xcol: str,
        ycol: str,
        top_taxa,
        undefined_val) -> pd.DataFrame:
    # sum smaller taxa into 'other' group
    top_taxa_indices = comp_df[tax_col].isin(top_taxa + [undefined_val])
    other_sum = comp_df[~top_taxa_indices].groupby(by=xcol)[ycol].sum()
    other_df = other_sum.reset_index()
    other_df[tax_col] = 'Other'
    
    comp_df = comp_df[comp_df[tax_col].isin(top_taxa)]
    comp_df = comp_df.reset_index(drop=True)
    comp_df = comp_df.sort_values(by=tax_col)
    comp_df[tax_col] = comp_df[tax_col].astype(str)
    
    comp_df = pd.concat([comp_df, other_df]) 

    return comp_df


def add_missing_cats(comp_df: pd.DataFrame, 
                     facet_comp_df: pd.DataFrame,
                     rank_display_col: str,
                     rank_hidden_col: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # ensure both dataframes have all categories (taxa and samples), even if they are 0
    cats = pd.MultiIndex.from_frame(
        pd.concat([comp_df[["Sample Name", rank_display_col, rank_hidden_col]], 
                   facet_comp_df[["Sample Name", rank_display_col, rank_hidden_col]]])
        .drop_duplicates()
    )
    comp_df = (comp_df
        .set_index(["Sample Name", rank_display_col, rank_hidden_col])
        .reindex(cats, fill_value=0)
        .reset_index()
    )
    facet_comp_df = (facet_comp_df
        .set_index(["Sample Name", rank_display_col, rank_hidden_col])
        .reindex(cats, fill_value=0)
        .reset_index()
    )
    return (comp_df, facet_comp_df)


def order_comp_categories(
    comp_df: pd.DataFrame,
    tax_col: str,
    undefined_val: Any) -> pd.DataFrame:
    # order categories to have 'Other' and 'Undefined' at the end
    for custom_cat in ['Other', undefined_val]:
        if custom_cat in comp_df[tax_col].values:
            other_df = comp_df[comp_df[tax_col] != custom_cat]
            custom_cat_df = comp_df[comp_df[tax_col] == custom_cat]
            comp_df = pd.concat([other_df, custom_cat_df])
    return comp_df


def add_hidden_data_col(
    comp_df: pd.DataFrame,
    tax_col: str,
    hidden_col: str,
    match_dict: Dict[str, str]) -> pd.DataFrame:
    comp_df.loc[:, hidden_col] = comp_df[tax_col].apply(
        lambda x: match_dict.get(x, np.nan)
    )
    return comp_df


def set_color_scale(
    topn: int,
    ncats: int,
    has_undefined: bool) -> List[str]:
    # set correct discrete color scale depending on number of taxa
    color_scale = deepcopy(GraphConstants.wide_color_palette) if topn > 10 \
        else deepcopy(GraphConstants.color_palette)
    
    # assign clear distinct color for 'Undefined group
    if has_undefined is True:
        # if more categories than scale + undef color, do nothing
        if ncats <= len(color_scale):
            color_scale[ncats - 1] = GraphConstants.undefined_color
        elif ncats == len(color_scale) + 1:
            color_scale.append(GraphConstants.undefined_color)

    return color_scale


def create_facet_barplot(
    comp_df: pd.DataFrame,
    facet_comp_df: pd.DataFrame,
    rank: str,
    rank_displ_col: str,
    rank_hid_col: str,
    xcol: str,
    ycol: str,
    color_scale: List[str],
    abundance_metric: str,
    facet_ycol: str,
    fractional_abundance: bool,
    facet_abundance_metric: str,
    facet_fractional_abundance: bool):
    share_y = (abundance_metric == facet_abundance_metric and 
               fractional_abundance == facet_fractional_abundance)
    
    spacing = 0.02 if share_y is True else 0.15

    # Create figure with two subplots, if y scaling is same, share axis
    fig = make_subplots(rows=1, cols=2, subplot_titles=["", "Facet plot"],
                        shared_yaxes=share_y, horizontal_spacing=spacing)

    # build the barplot traces for both facets, the facet plot is on the second column.
    for i, df in enumerate([comp_df, facet_comp_df]):
        show_legend = True if i == 0 else False
        sub_ycol = ycol if i == 0 else facet_ycol
        # add column to show if it is a facet df for customdata
        df["Is facet"] = False if i == 0 else True

        for j, tax_value in enumerate(df[rank_displ_col].unique()):
            cat_df = df[df[rank_displ_col] == tax_value]
            legend_name = truncate_end(tax_value, 30)
            fig.add_trace(
                go.Bar(
                    x=cat_df[xcol],
                    y=cat_df[sub_ycol],
                    name=legend_name,
                    customdata=cat_df[[rank_hid_col, rank_displ_col, 'Is facet']],
                    marker_color=color_scale[j],
                    showlegend=show_legend,
                    legendgroup=f"{tax_value}",
                    hovertemplate=(
                        f"{xcol}: %{{x}}<br>"
                        f"{ycol}: %{{y}}<br>"
                        f"{rank_displ_col}: %{{customdata[1]}}<br>"
                        f"{rank_hid_col}: %{{customdata[0]}}"
                    )
                ),
                row=1,
                col=i+1
            )

    fig.update_layout(GraphConstants.default_layout)
    fig.update_layout(margin=dict(l=20, r=20, t=30, b=10))
    fig.update_layout(
        barmode='stack', 
        #legend_groupclick='toggleitem',
        #legend_itemsizing='trace',
        legend_tracegroupgap=0,
        legend_title_text=rank_displ_col)
    
    fig.update_xaxes(showline=True, 
                     linecolor="Black", 
                     categoryorder="category ascending")
    
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     nticks=6)
    # update y axis based on abundance settings for each plot
    # configure y-axis title
    ytitle = 'Peptide spectrum matches' if abundance_metric == 'Match Count' else 'Area'
    if fractional_abundance is True:
        ytitle = "Fraction " + ytitle.lower()
    fig.update_yaxes(title=ytitle,
                     col=1)

    if share_y is False:
        facet_ytitle = 'Peptide spectrum matches' if facet_abundance_metric == 'Match Count' \
            else 'Area'
        if facet_fractional_abundance is True:
            facet_ytitle = "Fraction " + facet_ytitle.lower()
        fig.update_yaxes(title=facet_ytitle,
                         title_standoff=5,
                         col=2)

    return fig
