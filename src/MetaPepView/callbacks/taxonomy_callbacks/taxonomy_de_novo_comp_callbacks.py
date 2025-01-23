from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from server import app

# import layout elements
from layout.style_constants import *
from layout.validation_page import ms_performance
from layout.taxonomy_page import taxonomy_sample_analysis
from layout.func_annot_page import *
from layout.data_page import data_visual
from layout.header import content_header

from backend import *
from backend.plots import taxonomic_abundance_barplot,\
    taxonomic_abundance_heatmap,\
    tax_differential_barplot
from constants import *

import pandas as pd
    
    
@app.callback(
    Output('taxa_barplot_de_novo_graph', 'children'),
    Output('taxa_barplot_de_novo_graph', 'style'),
    Output('taxonomy_de_novo_figure_title', 'children'),
    Input('sidebar_taxonomy_de_novo_button', 'active'),
    Input('peptides', 'data'),
    Input("barplot_de_novo_sample_items", "value"),
    Input("taxonomy_stacked_barplot_button", "active"),
    Input('barplot_custom_taxa_items', 'value'),
    Input('tax_barplot_clade_selection_taxa', 'value'),
    Input('tax_barplot_clade_selection_rank', 'value'),
    Input('barplot_taxa_selector_radio', 'value'),
    Input('barplot_taxa_rank_items', 'value'),
    Input('barplot_taxa_fraction_checkbox', 'value'),
    Input('barplot_taxa_unannotated_checkbox', 'value'),
    Input('global_annot_de_novo_only_checkbox', 'value')
)
def update_de_novo_taxa_graph(page_active,
                              peptide_json,
                              sample_name,
                              bar_graph,
                              tax_ids,
                              filter_clade,
                              clade_rank,
                              top_taxa,
                              tax_rank,
                              fractional,
                              unannotated,
                              glob_annot_de_novo_only):
    if page_active is False:
        raise PreventUpdate
    
    if peptide_json is None:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Import PSM and protein db datasets...")
        return block_element, dict(), 'Figure'
    
    if sample_name is None:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select sample name from the dropdown menu...")
        return block_element, dict(), 'Figure'

    if tax_rank is None:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select taxonomy rank in dropdown menu...")
        return block_element, dict(), 'Figure'
    
    # fetch data from sample
    peptide_df = MetaPepTable.read_json(peptide_json).data
    peptide_df = peptide_df[peptide_df['Sample Name'] == sample_name]
    
    # reshape data: change sample name to metagenome annotation or unipept annotation
    peptide_df = reshape_taxonomy_df_to_denovo(peptide_df, glob_annot_de_novo_only)
    
    if filter_clade and clade_rank and clade_rank != 'Root':
        # filter the dataset based on taxa at corresponding rank
        peptide_df = peptide_df[peptide_df[clade_rank + " Name"] == filter_clade]
        plot_title = f"{sample_name} Taxonomic abundances from {filter_clade} clade, {tax_rank} rank"
    else:
        plot_title = f"{sample_name} Taxonomic abundances, {tax_rank} rank"
    
    if bar_graph is True:
        plot_method = taxonomic_abundance_barplot
    else:
        plot_method = taxonomic_abundance_heatmap
        
    
    if tax_ids != [] and top_taxa == 2:
        peptide_df = peptide_df[peptide_df[tax_rank + ' Name'].isin(tax_ids)]
        comp_plot = plot_method(peptide_df, 
                                rank=tax_rank,
                                fractional_abundance=fractional)
        dif_plot = tax_differential_barplot(peptide_df,
                                            'Global annotation',
                                            'Metagenome annotation',
                                            tax_rank)
    elif top_taxa == 2:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select custom tax id's...")
        return block_element, dict(), 'Figure'
    else:
        n_taxa = 9 if top_taxa == 1 else 24
        comp_plot = plot_method(peptide_df,
                                topn=n_taxa,
                                rank=tax_rank,
                                fractional_abundance=fractional,
                                include_undefined=unannotated)
        dif_plot = tax_differential_barplot(peptide_df,
                                            'Global annotation',
                                            'Metagenome annotation',
                                            tax_rank,
                                            topn=n_taxa,
                                            show_legend=False)

    
    graphs = [
        dcc.Graph(figure=comp_plot,
                  id="taxonomy_barplot_de_novo_figure",
                  style={"height": "28rem"}),
        html.Hr(),
        dcc.Graph(figure=dif_plot,
                  id="taxonomy_dif_barplot_de_novo_figure",
                  className="mt-2",
                  style={'height': '28rem'}),
    ]
    
    return graphs, dict(), plot_title
