from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

from metapepview.server import app

# import layout elements
from metapepview.html_templates import hidden_graph_with_text
from metapepview.backend import *
from metapepview.backend.plots import \
    taxonomic_abundance_barplot, \
    taxonomic_abundance_heatmap, \
    taxonomy_dropoff_scatter
from metapepview.constants import GlobalConstants as gc


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
    Input('barplot_taxa_quantification_column', 'value'),
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
                      quant_method,
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
        plot = plot_method(peptide_df, 
                           rank=tax_rank,
                           abundance_metric=quant_method,
                           fractional_abundance=fractional)
    elif top_n == 2:
        block_element = hidden_graph_with_text("taxonomy_barplot_figure",
                                               "Select custom tax id's...")
        return block_element, dict(), 'Figure'
    else:
        n_taxa = 9 if top_n == 1 else 20
        plot = plot_method(peptide_df, 
                           topn=n_taxa, 
                           rank=tax_rank,
                           abundance_metric=quant_method,
                           fractional_abundance=fractional,
                           include_undefined=unannotated)
    # plot.update_layout(height=500)

    return dcc.Graph(figure=plot,
                     id="taxonomy_barplot_figure",
                     style={'height': "45rem"}), dict(), plot_title


@app.callback(
    Output('taxonomic_dropoff_modal', 'is_open'),
    Output('taxonomic_dropoff_graph', 'style'),
    Output('taxonomic_dropoff_figure', 'figure'),
    Output('taxonomic_dropoff_title', 'children'),
    Output('lineage_tax_dropoff_text', 'children'),
    Output('global_av_dropoff_text', 'children'),
    Output('taxonomy_cumulative_dropoff_rank', 'options'),
    Output('taxonomy_cumulative_dropoff_rank', 'value'),
    Input('taxonomy_barplot_figure', 'clickData'),
    State('sidebar_taxonomy_button', 'active'),
    State('peptides', 'data'),
    State('barplot_taxa_rank_items', 'value'),
    State('barplot_taxa_quantification_column', 'value'),
    Input('taxonomic_dropoff_normalize', 'value'),
    Input('taxonomy_cumulative_dropoff_rank', 'value'),
    State('barplot_taxa_allow_global_annot_checkbox', 'value'),
    prevent_initial_call=True
)
def update_taxonomy_dropoff_graph(clickData,
                                  page_active,
                                  peptide_json,
                                  tax_rank,
                                  quant_method,
                                  normalize_bars,
                                  dropoff_root_rank,
                                  global_annot_fallback):
    if page_active is False \
        or clickData is None \
        or tax_rank is None \
        or peptide_json is None \
        or 'customdata' not in clickData['points'][0] \
        or gc.show_advanced_settings is False:
        return (False,
                {'display': 'None'},
                go.Figure(),
                None,
                "",
                "",
                [],
                "Root")

    # extract data from clickdata
    datapoint = clickData['points'][0]
    sample_name = datapoint['x']
    tax_id = datapoint['customdata'][0]

    # can only construct figure from valid tax id's.
    if tax_id is None:
        return (False,
                {'display': 'None'},
                go.Figure(),
                None,
                "",
                "",
                [],
                "Root")

    peptide_df = MetaPepTable.read_json(peptide_json).data
    peptide_df = peptide_df[peptide_df["Sample Name"] == sample_name]

    # substitute missing taxonomy annotation with global annotation if specified
    glob_tax_fields = GlobalConstants.metapep_table_global_taxonomy_lineage
    if global_annot_fallback is True and \
        all(i in peptide_df.columns for i in glob_tax_fields):
        peptide_df = substitute_lineage_with_global_lineage(peptide_df)

    # select column to sum, match count or total signal
    if quant_method == "Match Count":
        quant_col = "PSM Count"
    else:
        quant_col = "Area"

    # fetch lineage directly from peptide json
    rank_names = GlobalConstants.standard_lineage_ranks
    lineage_cols = [rank + ' Name' for rank in rank_names]
    lineage = peptide_df[peptide_df[tax_rank + ' Id'] == tax_id]

    # stop execution if no valid lineage found, else, retrieve first row and convert to list
    if lineage.shape[0] == 0:
        raise PreventUpdate
    else:
        tax_name = lineage.iloc[0][tax_rank + ' Name']
        lineage = lineage.iloc[0][lineage_cols].fillna("-").to_list()

        # cut off any rand definitions beyond rank name
        rank_idx = GlobalConstants.standard_lineage_ranks.index(tax_rank)
        for i in range(rank_idx + 1, len(lineage)):
            lineage[i] = "-"

    lineage_counts, pept_allocation = peptide_allocation_across_lineage(
        peptide_df,
        lineage,
        quant_col)

    # configure dropoff dropdown menu and check that current value is part of dropdown menu
    dropoff_menu_options = [
        {'label': 'Root', 'value':"Root"}
    ]

    for rank, rank_short, clade in zip(GlobalConstants.standard_lineage_ranks,
                                       GlobalConstants.lineage_ranks_short,
                                       lineage):
        if clade != "-":
            clade = clade
            dropoff_menu_options.append({'label': clade, 'value': rank})
    if dropoff_root_rank not in [x['value'] for x in dropoff_menu_options]:
        dropoff_root_rank = "Root"

    # get clade from root rank
    if dropoff_root_rank == "Root":
        dropoff_root_clade = "Root"
    else:
        root_rank_idx = GlobalConstants.standard_lineage_ranks.index(dropoff_root_rank)
        dropoff_root_clade = lineage[root_rank_idx]

    # compute cumulative dropoffs
    if dropoff_root_clade == dropoff_root_clade and dropoff_root_clade != "-":
        combined_lin_loss = compute_lineage_cumulative_annotation_drop(
            pept_allocation,
            dropoff_root_rank
        )
        if combined_lin_loss != combined_lin_loss:
            lin_loss_text = "-".format(combined_lin_loss)
        else:
            lin_loss_text = "{:.1f}%".format(combined_lin_loss)

        combined_glob_loss = compute_global_cumulative_annotation_drop(
            peptide_df,
            dropoff_root_rank,
            dropoff_root_clade,
            tax_rank,
            quant_col
        )
        if combined_glob_loss != combined_glob_loss:
            glob_loss_text = " -".format(combined_glob_loss)
        else:
            glob_loss_text = "{:.1f}%".format(combined_glob_loss)
    else:
        lin_loss_text, glob_loss_text = "-", "-"


    fig = taxonomy_dropoff_scatter(
        peptide_df,
        lineage_counts,
        pept_allocation,
        normalize_bars)

    title = "Taxonomic dropoff: {} from {}".format(tax_name, sample_name)

    return (True,
            dict(),
            fig,
            title,
            lin_loss_text,
            glob_loss_text,
            dropoff_menu_options,
            dropoff_root_rank)
