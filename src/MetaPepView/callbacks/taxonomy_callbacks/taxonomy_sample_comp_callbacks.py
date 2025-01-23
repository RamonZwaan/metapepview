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
from backend.plots import taxonomic_abundance_barplot, taxonomic_abundance_heatmap
from constants import *
    
    
@app.callback(
    Output('taxa_barplot_graph', 'children'),
    Output('taxa_barplot_graph', 'style'),
    Output('taxonomy_figure_title', 'children'),
    Input('sidebar_taxonomy_button', 'active'),
    Input('peptides', 'data'),
    Input("taxonomy_stacked_barplot_button", "active"),
    Input('barplot_custom_taxa_items', 'value'),
    Input('tax_barplot_clade_selection_taxa', 'value'),
    Input('tax_barplot_clade_selection_rank', 'value'),
    Input('barplot_taxa_selector_radio', 'value'),
    Input('barplot_taxa_rank_items', 'value'),
    Input('barplot_taxa_fraction_checkbox', 'value'),
    Input('barplot_taxa_unannotated_checkbox', 'value'),
    Input('barplot_taxa_allow_global_annot_checkbox', 'value')
)
def update_taxa_graph(page_active,
                      peptide_json,
                      bar_graph,
                      tax_ids,
                      filter_clade,
                      clade_rank,
                      top_n,
                      tax_rank,
                      fractional,
                      unannotated,
                      global_annot_fallback):
    if page_active is False:
        raise PreventUpdate
    
    if peptide_json is None:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Import PSM and protein db datasets...")
        return block_element, dict(), 'Figure'
    
    if tax_rank is None:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select taxonomy rank in dropdown menu...")
        return block_element, dict(), 'Figure'
    
    peptide_df = MetaPepTable.read_json(peptide_json).data
    
    # substitute missing taxonomy annotation with global annotation if specified
    glob_tax_fields = GlobalConstants.metapep_table_global_taxonomy_lineage
    if global_annot_fallback is True and \
        all(i in peptide_df.columns for i in glob_tax_fields):
        peptide_df = substitute_lineage_with_global_lineage(peptide_df)        
        
    
    if filter_clade and clade_rank and clade_rank != 'Root':
        # filter the dataset based on taxa at corresponding rank
        peptide_df = peptide_df[peptide_df[clade_rank + " Name"] == filter_clade]
        plot_title = f"Taxonomic abundances from {filter_clade} clade, {tax_rank} rank"
    else:
        plot_title = f"Taxonomic abundances, {tax_rank} rank"
    
    if bar_graph is True:
        plot_method = taxonomic_abundance_barplot
    else:
        plot_method = taxonomic_abundance_heatmap
        
    
    if tax_ids != [] and top_n == 2:
        peptide_df = peptide_df[peptide_df[tax_rank + ' Name'].isin(tax_ids)]
        plot = plot_method(peptide_df, rank=tax_rank, fractional_abundance=fractional)
    elif top_n == 2:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select custom tax id's...")
        return block_element, dict(), 'Figure'
    else:
        n_taxa = 9 if top_n == 1 else 24
        plot = plot_method(peptide_df, topn=n_taxa, rank=tax_rank, fractional_abundance=fractional, include_undefined=unannotated)
    # plot.update_layout(height=500)
    
    return dcc.Graph(figure=plot,
                     id="taxonomy_barplot_figure",
                     style={'height': "45rem"}), dict(), plot_title


