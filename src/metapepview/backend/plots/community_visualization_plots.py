import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# import dash_bio

from typing import List
from copy import deepcopy

import pandas as pd
import numpy as np

from metapepview.backend.spectral_ref_builder import truncate_end
from metapepview.backend.utils.graph_utils import *
from metapepview.backend.utils.pd_utils import *
from metapepview.backend.io.import_spectra import *
from metapepview.constants import *


from metapepview.backend.utils.taxonomy_plot_utils import *


# assign template theme for plotly figures
pio.templates.default = GraphConstants.default_template


# taxonomy abundance composition
def taxonomic_abundance_barplot(peptide_dataset: pd.DataFrame,
                                rank: RankType='Phylum',
                                abundance_metric: AbundanceMetricType='Match Count',
                                fetch_names: bool=True,
                                include_undefined: bool=False,
                                fractional_abundance: bool=False,
                                topn: int=0):
    # Set column names to process for the figure
    xcol = 'Sample Name'

    # configure taxonomy columns to process based on rank and selected display format
    (display_suffix, hidden_suffix) = (' Name', ' Id') if fetch_names is True \
        else (' Id', ' Name')
    rank_display_col = rank + display_suffix    
    rank_hidden_col = rank + hidden_suffix

    # generate lookup table to match display format (id/name) to its corresponding hidden format
    rank_cols = [rank + sfx for sfx in (display_suffix, hidden_suffix)]
    id_name_match = {displ_name: hidden_name for displ_name, hidden_name \
        in peptide_dataset[rank_cols].drop_duplicates().dropna().values}
    
    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_display_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0

    comp_df, ycol = process_tax_abundance(
        peptide_dataset,
        abundance_metric,
        rank_display_col,
        include_undefined,
        undefined_val,
        fractional_abundance
    )
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = determine_top_taxa(
            comp_df=comp_df,
            facet_df=None,
            topn=topn,
            tax_col=rank_display_col,
            ycol=ycol,
            facet_ycol=None,
            undefined_val=undefined_val,
            include_undefined=include_undefined)
        # summ smaller taxa into 'other' group
        comp_df = add_other_group(
            comp_df, 
            rank_display_col, 
            xcol, 
            ycol, 
            top_taxa, 
            undefined_val
        )  
        
    # whether taxonomy id or name is displayed, include the other as hidden data
    comp_df = add_hidden_data_col(comp_df, rank_display_col, rank_hidden_col, id_name_match)
    
    # order categories to have 'Other' and 'Undefined' at the end
    comp_df = order_comp_categories(comp_df, rank_display_col, undefined_val)
    
    # set correct discrete color scale depending on number of taxa
    color_scale = deepcopy(GraphConstants.wide_color_palette) if topn > 10 \
        else deepcopy(GraphConstants.color_palette)

    # configure colorscale to use
    ncats = len(comp_df[rank_display_col].unique())
    has_undefined = 'Undefined' in comp_df[rank_display_col].values
    color_scale = set_color_scale(topn, ncats, has_undefined)

    fig = px.bar(comp_df,
                 title=f"Distribution PSM over taxa ({rank})",
                 x=xcol,
                 y=ycol,
                 color=rank_display_col,
                 custom_data=[rank_hidden_col, rank_display_col],
                 color_discrete_sequence=color_scale)
    
    fig.update_layout(GraphConstants.default_layout)
    fig.update_layout(title="")
     
    # configure y-axis title
    ytitle = 'Peptide spectrum matches' if abundance_metric == 'Match Count' else 'Area'
    if fractional_abundance is True:
        ytitle = "Fraction " + ytitle.lower()

    # configure hovertemplate (Essential elements to be presented in clickdata)
    fig.update_traces(
        hovertemplate=f"{xcol}: %{{x}}<br>{ycol}: %{{y}}<br>{rank_display_col}: %{{customdata[1]}}<br>{rank_hidden_col}: %{{customdata[0]}}"
    )
    
    fig.update_xaxes(showline=True, 
                     linecolor="Black", 
                     categoryorder="category ascending")
    
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     nticks=5,
                     title=ytitle)
    
    # fig.update_layout(width=900, height=500)
    return fig



def facet_taxonomic_abundance_barplot(peptide_dataset: pd.DataFrame,
                                      facet_dataset: pd.DataFrame | None,
                                      rank: RankType='Phylum',
                                      abundance_metric: AbundanceMetricType='Match Count',
                                      fetch_names: bool=True,
                                      include_undefined: bool=False,
                                      fractional_abundance: bool=False,
                                      facet_abundance_metric: AbundanceMetricType='Match Count',
                                      facet_include_undefined: bool=False,
                                      facet_fractional_abundance: bool=False,
                                      topn: int=0):
    # Set column names to process for the figure
    xcol = 'Sample Name'

    # configure taxonomy columns to process based on rank and selected display format
    (display_suffix, hidden_suffix) = (' Name', ' Id') if fetch_names is True \
        else (' Id', ' Name')
    rank_display_col = rank + display_suffix    
    rank_hidden_col = rank + hidden_suffix

    # generate lookup table to match display format (id/name) to its corresponding hidden format
    rank_cols = [rank + sfx for sfx in (display_suffix, hidden_suffix)]
    id_name_match = {displ_name: hidden_name for displ_name, hidden_name \
        in peptide_dataset[rank_cols].drop_duplicates().dropna().values}

    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_display_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0
    
    # facet dataset may be already provided (because Unipept suppl. was performed)
    if facet_dataset is None:
        facet_dataset = deepcopy(peptide_dataset)

    comp_df, ycol = process_tax_abundance(
        peptide_dataset,
        abundance_metric,
        rank_display_col,
        include_undefined,
        undefined_val,
        fractional_abundance
    )
    facet_comp_df, facet_ycol = process_tax_abundance(
        facet_dataset,
        facet_abundance_metric,
        rank_display_col,
        facet_include_undefined,
        undefined_val,
        facet_fractional_abundance
    )
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = determine_top_taxa(comp_df,
                                      facet_comp_df,
                                      topn,
                                      rank_display_col,
                                      ycol,
                                      facet_ycol,
                                      undefined_val,
                                      any([include_undefined, facet_include_undefined]))
        comp_df = add_other_group(
            comp_df, rank_display_col, xcol, ycol, top_taxa, undefined_val
        ) 
        facet_comp_df = add_other_group(
            facet_comp_df, rank_display_col, xcol, facet_ycol, top_taxa, undefined_val
        ) 
        
    # whether taxonomy id or name is displayed, include the other as hidden data
    comp_df = add_hidden_data_col(comp_df, rank_display_col, rank_hidden_col, id_name_match)
    facet_comp_df = add_hidden_data_col(facet_comp_df, rank_display_col, rank_hidden_col, id_name_match)
    
    # ensure both dataframes have all categories (taxa and samples), even if they are 0
    comp_df, facet_comp_df = add_missing_cats(
        comp_df, facet_comp_df, rank_display_col, rank_hidden_col
    )

    # order categories to have 'Other' and 'Undefined' at the end
    comp_df = order_comp_categories(comp_df, rank_display_col, undefined_val)
    facet_comp_df = order_comp_categories(facet_comp_df, rank_display_col, undefined_val)
    

    # configure colorscale to use
    ncats = max(len(comp_df[rank_display_col].unique()),
                len(facet_comp_df[rank_display_col].unique()))
    has_undefined = ('Undefined' in comp_df[rank_display_col].values or 
                     'Undefined' in facet_comp_df[rank_display_col].values)
    color_scale = set_color_scale(topn, ncats, has_undefined)

    fig = create_facet_barplot(
        comp_df=comp_df,
        facet_comp_df=facet_comp_df,
        rank=rank,
        rank_displ_col=rank_display_col,
        rank_hid_col=rank_hidden_col,
        xcol=xcol,
        ycol=ycol,
        color_scale=color_scale,
        abundance_metric=abundance_metric,
        facet_ycol=facet_ycol,
        fractional_abundance=fractional_abundance,
        facet_abundance_metric=facet_abundance_metric,
        facet_fractional_abundance=facet_fractional_abundance
    )
    
    # fig.update_layout(width=900, height=500)
    return fig




# taxonomy abundance composition
def taxonomic_abundance_heatmap(peptide_dataset: pd.DataFrame,
                                rank: RankType="Phylum",
                                abundance_metric: AbundanceMetricType="Match Count",
                                fetch_names: bool=True,
                                include_undefined: bool =False,
                                fractional_abundance: bool=False,
                                topn: int=0):
    # column name for abundances
    xcol = 'Sample Name'
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'
    # get correct suffix from lineage columns
    suffix = ' Name' if fetch_names is True else ' Id'
    rank_col = rank + suffix            
    
    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0

    # tax_comp_list = list()
    
    # # get composition data for each sample
    # for i, sample_name in enumerate(peptide_dataset[xcol].unique()):
    #     dataset = peptide_dataset[peptide_dataset[xcol] == sample_name]
            
    #     # convert dataset to rank name and their aggregated sum
    #     rank_composition = dataset[[ycol, rank_col]]
        
    #     # give value to all nan cells to represent undefined or unknown
    #     if include_undefined is True:
    #         rank_composition[rank_col].fillna(undefined_val, inplace=True)
                
    #     # group all peptides belonging to rank, sum separate abundance scores
    #     rank_composition = rank_composition.groupby(by=rank_col).sum().reset_index()
        
    #     if fractional_abundance is True:
    #         rank_composition[ycol] /= rank_composition[ycol].sum()
        
    #     rank_composition.loc[:, xcol] = sample_name
        
    #     tax_comp_list.append(rank_composition)

    # # combine samples into one datast
    # comp_df = pd.concat(tax_comp_list)

    comp_df, ycol = process_tax_abundance(
        peptide_dataset=peptide_dataset,
        abundance_metric=abundance_metric,
        tax_col=rank_col,
        include_undefined=include_undefined,
        undefined_val=undefined_val,
        fractional_abundance=fractional_abundance
    )
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = determine_top_taxa(
            comp_df=comp_df,
            facet_df=None,
            topn=topn,
            tax_col=rank_col,
            ycol=ycol,
            facet_ycol=None,
            undefined_val=undefined_val,
            include_undefined=include_undefined
        )
        
        # in contrast to barplot, no 'Other' is added, just the top taxa (+ undefined)
        comp_df = comp_df[comp_df[rank_col].isin(top_taxa)]
        comp_df[rank_col] = comp_df[rank_col].astype(str)
    
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


# taxonomy abundance composition
def facet_taxonomic_abundance_heatmap(peptide_dataset: pd.DataFrame,
                                      facet_dataset: pd.DataFrame | None,
                                      rank: RankType="Phylum",
                                      abundance_metric: AbundanceMetricType="Match Count",
                                      fetch_names: bool=True,
                                      include_undefined: bool =False,
                                      fractional_abundance: bool=False,
                                      facet_abundance_metric: AbundanceMetricType='Match Count',
                                      facet_include_undefined: bool=False,
                                      facet_fractional_abundance: bool=False,
                                      topn: int=0):
    # column name for abundances
    xcol = 'Sample Name'
    ycol = 'PSM Count' if abundance_metric == 'Match Count' else 'Area'
    # get correct suffix from lineage columns
    suffix = ' Name' if fetch_names is True else ' Id'
    rank_col = rank + suffix            
    

    # Set replacement values for specific cases
    # value of undefined depends on dtype of column (id number or name string)
    if pd.api.types.is_string_dtype(peptide_dataset[rank_col].dtype):
        undefined_val = "Undefined"
    else:
        undefined_val = 0

    # facet dataset may be already provided (because Unipept suppl. was performed)
    if facet_dataset is None:
        facet_dataset = deepcopy(peptide_dataset)

    comp_df, ycol = process_tax_abundance(
        peptide_dataset=peptide_dataset,
        abundance_metric=abundance_metric,
        tax_col=rank_col,
        include_undefined=include_undefined,
        undefined_val=undefined_val,
        fractional_abundance=fractional_abundance
    )
    facet_comp_df, facet_ycol = process_tax_abundance(
        peptide_dataset=facet_dataset,
        abundance_metric=facet_abundance_metric,
        tax_col=rank_col,
        include_undefined=facet_include_undefined,
        undefined_val=undefined_val,
        fractional_abundance=facet_fractional_abundance
    )
    
    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = determine_top_taxa(
            comp_df=comp_df,
            facet_df=facet_comp_df,
            topn=topn,
            tax_col=rank_col,
            ycol=ycol,
            facet_ycol=facet_ycol,
            undefined_val=undefined_val,
            include_undefined=any([include_undefined, facet_include_undefined])
        )

        # in contrast to barplot, no 'Other' is added, just the top taxa (+ undefined)
        comp_df = comp_df[comp_df[rank_col].isin(top_taxa)]
        comp_df[rank_col] = comp_df[rank_col].astype(str)

        facet_comp_df = facet_comp_df[facet_comp_df[rank_col].isin(top_taxa)]
        facet_comp_df[rank_col] = facet_comp_df[rank_col].astype(str)
    
    # ensure both dataframes have all categories (taxa and samples), even if they are 0
    comp_df, facet_comp_df = add_missing_cats(comp_df, facet_comp_df, rank_col, fill_value=np.nan)

    value_matrix = comp_df[[ycol, xcol, rank_col]].pivot(index=xcol,
                                                         columns=rank_col,
                                                         values=ycol)
    facet_value_matrix = facet_comp_df[[facet_ycol, xcol, rank_col]].pivot(index=xcol,
                                                                           columns=rank_col,
                                                                           values=facet_ycol)
    
    fig = create_facet_heatmap(
        comp_df=value_matrix,
        facet_comp_df=facet_value_matrix,
        rank_col=rank_col,
        color_scale=GraphConstants.continuous_color_scale,
        abundance_metric=abundance_metric,
        fractional_abundance=fractional_abundance,
        facet_abundance_metric=facet_abundance_metric,
        facet_fractional_abundance=facet_fractional_abundance
    )

    # # create simple annotated heatmap
    # fig = px.imshow(value_matrix,
    #                 color_continuous_scale=GraphConstants.continuous_color_scale)
    
    # # configure x-axis title
    # xtitle = 'Peptide spectrum matches' if abundance_metric == 'Match Count' else 'Area'
    # if fractional_abundance is True:
    #     xtitle = "Fraction " + xtitle.lower()
    
    # fig.update_xaxes(side="top",
    #                  title=xtitle)
    # fig.update_layout(GraphConstants.default_layout)
    
    return fig


def taxonomy_dropoff_scatter(lineage_counts: List[Tuple[str, int | float]],
                             lineage_dropoff: List[Tuple[str, 
                                                         Tuple[float,
                                                               float,
                                                               float,
                                                               float,
                                                               Dict[str, 
                                                                    float]]]],
                             normalize_bars: bool = False,
                             ytitle: str = "Peptide spectrum matches"):
    
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
        else:
            other_branches_names += [[]]
            other_branches_values += [np.array([])]
            
    
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
    
    

    fig.update_yaxes(title=ytitle,
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
                             fractional_abundance: bool=False,
                             fractional_abundance_threshold: float=0.001,
                             show_legend: bool=True):
    # Set column names to process for the figure
    sample_name = 'Sample Name'

    # get correct suffix from lineage columns
    suffix = ' Name' if fetch_names is True else ' Id'
    rank_col = rank + suffix
    
    # set color scale based on topn count
    color_scale = GraphConstants.wide_color_palette if topn > 10 else GraphConstants.color_palette
    
    # For top taxa determination, calculate tax abundances with same settings as barplot
    # TODO: export top taxa from barplot to dcc.Store, then import them here
    top_tax_comp_df, _ = process_tax_abundance(
        peptide_dataset=peptide_dataset,
        abundance_metric=abundance_metric,
        tax_col=rank_col,
        include_undefined=False,
        undefined_val=None,
        fractional_abundance=fractional_abundance
    )
    comp_df, ycol = process_tax_abundance(
        peptide_dataset=peptide_dataset,
        abundance_metric=abundance_metric,
        tax_col=rank_col,
        include_undefined=False,
        undefined_val=None,
        fractional_abundance=True
    )
    

    # only keep the most abundant taxa over both samples, optionally, max is taken instead of sum
    if topn > 0:
        top_taxa = determine_top_taxa(
            comp_df=top_tax_comp_df,
            facet_df=None,
            topn=topn,
            tax_col=rank_col,
            ycol=ycol,
            facet_ycol=None,
            undefined_val=None,
            include_undefined=False
        )
        
        comp_df = filter_dataset_top_taxa(comp_df=comp_df,
                                          tax_col=rank_col,
                                          top_taxa=top_taxa)
    
    # for each taxa, compute ratio between samples
    comp_pivot_df = comp_df.pivot(index=rank_col,
                                  columns=sample_name,
                                  values=ycol).fillna(0.0)
    
    # compute ratio of fractions by dividing the numerator sample with the denominator sample
    comp_pivot_df['Ratio'] = comp_pivot_df[ratio_numerator_col] / comp_pivot_df[ratio_denominator_col]
    comp_pivot_df['Size'] = comp_pivot_df[ratio_numerator_col] + comp_pivot_df[ratio_denominator_col]
    comp_pivot_df['yaxis scatter'] = "abundance (%)"
    
    # rescale ratio to (-inf, inf) to (-1, 1) with log10 and tanh function
    comp_pivot_df['Norm ratio'] = comp_pivot_df['Ratio'].apply(np.log10).apply(np.tanh)
    comp_pivot_df = comp_pivot_df.reset_index()
    
    # build figure
    fig = make_subplots(rows=2,
                        cols=1,
                        row_heights=[0.5, 0.5],
                        shared_xaxes=True,
                        vertical_spacing=0.01)
    
    fig.add_trace(
        go.Bar(x=comp_pivot_df[rank_col],
               y=comp_pivot_df['Norm ratio'],
               marker_color=color_scale,
               #marker_color=comp_pivot_df[rank_col],
               #color_discrete_sequence=GraphConstants.color_palette
        ),
        row=1,
        col=1
    )
    
    # take root to let marker area correspond with fraction, multiply to scale to plot
    marker_size_values = (comp_pivot_df['Size'] ** 0.5) * 30
    # give size fraction percentage as string value
    marker_size_str = (comp_pivot_df['Size'] * 100 / 2).apply(lambda x: '{0:.1f}'.format(x))
    
    fig.add_trace(
        go.Scatter(x=comp_pivot_df[rank_col],
                   y=comp_pivot_df['yaxis scatter'],
                   text=marker_size_str,
                   textposition="top center",
                   mode='markers+text',
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
                      margin=dict(l=20, r=20, t=35, b=20))
    fig.update_traces(width=0.5, row=1, col=1)
    fig.update_layout(title=dict(
                        text=f"Abundance ratio: {ratio_numerator_col} / {ratio_denominator_col}",
                        font_size=16)
                      )
    
    fig.update_xaxes(title=None)
    
    fig.update_layout(
        yaxis=dict(gridcolor=GraphConstants.gridcolor,
                   gridwidth=GraphConstants.gridwidth,
                   nticks=5,
                   zeroline=True, 
                   zerolinecolor="Black",
                   showline=True, linecolor="Black", 
                   title='norm log ratio'),
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
                      margin=dict(t=30),
                      legend_title_text=None)
    
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