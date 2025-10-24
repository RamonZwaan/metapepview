from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import pandas as pd

from metapepview.server import app

from metapepview.backend import *
from metapepview.constants import *

@app.callback(
    Output("barplot_de_novo_sample_items", "options"),
    Input('peptides', 'data')
)
def update_barplot_sample_selector(peptide_json):
    """Update sample to visualize in the db search vs de novo
    taxonomy comparison.
    """
    if peptide_json is None:
        return []
    
    # TODO: This data should be in peptide metadata table
    metapep_obj = MetaPepTable.read_json(peptide_json)
    metapep_df = metapep_obj.data

    # get sample names and if they have de novo taxonomy annotation
    glob_annot_samples = metapep_df[['Sample Name', 
                                     'Global Taxonomy Annotation', 
                                     'Taxonomy Annotation']]\
        .drop_duplicates(subset='Sample Name')
    
    # only samples with both taxonomy annotations should be included
    annot_samples = glob_annot_samples[
        (glob_annot_samples['Global Taxonomy Annotation'] == True) & \
        (glob_annot_samples['Taxonomy Annotation'] == True)]['Sample Name']
    return annot_samples.to_list()
     
        


@app.callback(
    Output(component_id="barplot_custom_taxa_items", component_property="disabled"),
    Input(component_id="barplot_taxa_selector_radio", component_property="value")
)
def update_barplot_taxa_selector(selector_radio_value):
    """Update the element used to select the given metabolic pathway.
    Either a list of predifined pathways can be given, or a multi choose
    menu to select your own genes.
    """
    if selector_radio_value == 1:
        return True
    else:
        return False
    
    
@app.callback(
    Output(component_id="heatmap_custom_taxa_items", component_property="disabled"),
    Input(component_id="heatmap_taxa_selector_radio", component_property="value")
)
def update_heatmap_taxa_selector(selector_radio_value):
    """Update the element used to select the given metabolic pathway.
    Either a list of predifined pathways can be given, or a multi choose
    menu to select your own genes.
    """
    if selector_radio_value == 1:
        return True
    else:
        return False
    
    
@app.callback(
    Output('barplot_custom_taxa_items', 'options'),
    Output('barplot_custom_taxa_items', 'value'),
    Input('peptides', 'data'),
    Input('barplot_taxa_rank_items', 'value'),
    Input('tax_barplot_clade_selection_taxa', 'value'),
    Input('tax_barplot_clade_selection_rank', 'value'),
)
def update_barplot_custom_taxa_items_options(peptide_json, tax_rank, filter_clade, clade_rank):
    if peptide_json is None:
        return ([], [])
    
    if tax_rank is None:
        return ([], [])
    
    metapep_obj = MetaPepTable.read_json(peptide_json)
    
    # return a list of taxonomy id's, sorted and without nan's
    return (unique_taxa_from_rank(metapep_obj, tax_rank, filter_clade, clade_rank, True),
        [])
    

@app.callback(
    Output('tax_barplot_clade_selection_taxa', 'options'),
    Output('tax_barplot_clade_selection_taxa', 'value'),
    Output('tax_barplot_clade_selection_taxa', 'disabled'),
    Output('tax_barplot_clade_selection_rank', 'disabled'),
    Input('peptides', 'data'),
    Input('tax_barplot_clade_selection_rank', 'value')
)
def update_barplot_clade_selection_options(peptide_json, tax_rank):
    if peptide_json is None:
        return ([], [], True, True)
    
    if not tax_rank or tax_rank == "Root":
        return ([], [], True, False)
    
    metapep_obj = MetaPepTable.read_json(peptide_json)
    
    # return a list of taxonomy id's, sorted and without nan's
    return (unique_taxa_from_rank(metapep_obj, tax_rank), [], False, False)
    
    
@app.callback(
    Output('heatmap_custom_taxa_items', 'options'),
    Output('heatmap_custom_taxa_items', 'value'),
    Input('peptides', 'data'),
    Input('heatmap_taxa_rank_items', 'value')
)
def update_heatmap_custom_taxa_items_options(peptide_json, tax_rank):
    if peptide_json is None:
        return ([], [])
    peptide_df = MetaPepTable.read_json(peptide_json).data
    
    # return a list of taxonomy id's, sorted and without nan's
    return (peptide_df[tax_rank + ' Name']\
        .dropna()\
        .drop_duplicates()\
        .sort_values()\
        .to_list(),
        [])
    
    
@app.callback(
    Output("taxonomy_stacked_barplot_button", "active"),
    Output("taxonomy_heatmap_button", "active"),
    Input("taxonomy_stacked_barplot_button", "n_clicks"), 
    Input("taxonomy_heatmap_button", "n_clicks"),
    prevent_initial_call=True
)
def toggle_bar_graph(bar_click, heat_click):
    if ctx.triggered_id == "taxonomy_stacked_barplot_button":
        return (True, False)
    else:
        return (False, True)


@app.callback(
    Output('export_taxonomy_button', 'disabled'),
    Input('peptides', 'data')
)
def toggle_export_button(peptide_json):
    if peptide_json is None:
        return True
    else:
        return False
    

@app.callback(
    Output('download_taxonomy_composition_csv', 'data'),
    Input('export_taxonomy_button', 'n_clicks'),
    State('peptides', 'data'),
    State('tax_barplot_clade_selection_taxa', 'value'),
    State('tax_barplot_clade_selection_rank', 'value'),
    State('barplot_taxa_allow_global_annot_checkbox', 'value'),
    prevent_initial_call=True
)
def export_tax_composition(button_click, 
                           peptide_json,
                           filter_clade,
                           clade_rank,
                           global_annot_fallback):
    if peptide_json is None:
        raise PreventUpdate
    
    # import peptide dataset
    peptide_df = MetaPepTable.read_json(peptide_json).data
    
    # define metadata for export
    metadata = {
        "Root Clade": "global root",
        "Root Rank": "global",
    }
    
    # substitute missing taxonomy annotation with global annotation if specified
    glob_tax_fields = GlobalConstants.metapep_table_global_taxonomy_lineage
    if global_annot_fallback is True and \
        all(i in peptide_df.columns for i in glob_tax_fields):
        peptide_df = substitute_lineage_with_global_lineage(peptide_df)    
    
    # filter dataset by selected root clade    
    if filter_clade is not None and clade_rank is not None and clade_rank != 'Root':
        # filter the dataset based on taxa at corresponding rank
        peptide_df = peptide_df[peptide_df[clade_rank + " Name"] == filter_clade]
        
        # add clade info to metadata
        metadata["Root Clade"] = filter_clade
        metadata["Root Rank"] = clade_rank
    
    # create dataframe that stores abundances for each taxonomic group
    if clade_rank != 'Root' and clade_rank is not None:
        lin_idx = GlobalConstants.standard_lineage_ranks.index(clade_rank)
    else:
        lin_idx = 0
        
    tax_comp_df = []
    for current_rank in GlobalConstants.standard_lineage_ranks[lin_idx:]:
        rank_groups = peptide_df[[current_rank + ' Name', 'Area', 'PSM Count', 'Sample Name']].groupby([current_rank + ' Name', 'Sample Name'])
        group_aggs = rank_groups[['Area', 'PSM Count']].agg('sum')
        group_aggs = group_aggs.rename_axis(['Taxonomy Name', 'Sample Name']).reset_index()
        group_aggs['Taxonomy Rank'] = current_rank
        tax_comp_df.append(group_aggs) 
    tax_comp_df = pd.concat(tax_comp_df).reset_index(drop=True)
    
    return dcc.send_data_frame(tax_comp_df.to_csv, "tax_composition.tsv", sep="\t")


@app.callback(
    Output('download_taxonomy_composition_de_novo_csv', 'data'),
    Input('export_taxonomy_eval_button', 'n_clicks'),
    State('peptides', 'data'),
    State("barplot_de_novo_sample_items", "value"),
    State('barplot_taxa_quantification_column', 'value'),
    State('global_annot_de_novo_only_checkbox', 'value'),
    prevent_initial_call=True
)
def export_de_novo_tax_composition(button_click, 
                                   peptide_json,
                                   sample_name,
                                   quant_method,
                                   glob_annot_de_novo_only):
    if peptide_json is None:
        raise PreventUpdate
    
    quant_col = 'PSM Count' if quant_method == 'Match Count' else 'Area'
    
    # import peptide dataset
    peptide_df = MetaPepTable.read_json(peptide_json).data
    peptide_df = peptide_df[peptide_df['Sample Name'] == sample_name]

    # reshape data: change sample name to metagenome annotation or unipept annotation
    peptide_df, db_search_col, unipept_col = reshape_taxonomy_df_to_denovo(
        peptide_df, 
        glob_annot_de_novo_only
    )
        
    tax_comp_df = []
    for current_rank in GlobalConstants.standard_lineage_ranks:
        rank_groups = (peptide_df[[current_rank + ' Name', 
                                  'Area', 
                                  'PSM Count', 
                                  'Sample Name']]
                       .groupby([current_rank + ' Name', 'Sample Name'])
        )
        group_aggs = rank_groups[['Area', 'PSM Count']].agg('sum')
        group_aggs = group_aggs.rename_axis(['Taxonomy Name', 'Sample Name']).reset_index()
        group_aggs['Taxonomy Rank'] = current_rank
        tax_comp_df.append(group_aggs) 
    tax_comp_df = pd.concat(tax_comp_df).reset_index(drop=True)

    # arrange two sample data next to each other
    pivot_df = tax_comp_df.pivot(index=["Taxonomy Name", "Taxonomy Rank"],
                                 columns="Sample Name",
                                 values=quant_col)
    pivot_df = pivot_df.reset_index()

    # if no values present in any sample name, replace with empty values    
    if db_search_col not in pivot_df.columns: pivot_df[db_search_col] = np.nan
    if unipept_col not in pivot_df.columns: pivot_df[unipept_col] = np.nan

    # calculate ratio in final column
    pivot_df[f"Ratio {db_search_col} / {unipept_col}"] = (pivot_df[db_search_col] / pivot_df[unipept_col]).round(decimals=3)
    
    # order data by rank and by total PSM (fractions) between samples
    quant_arg_sort = (pivot_df[db_search_col].fillna(0) + 
                      pivot_df[unipept_col].fillna(0)).argsort()[::-1]
    pivot_df = pivot_df.loc[quant_arg_sort, :]
    pivot_df = pivot_df.sort_values(by="Taxonomy Rank", 
                                    key=lambda x: x.apply(GlobalConstants.standard_lineage_ranks.index)
    )
    pivot_df = pivot_df.reset_index(drop=True)

    return dcc.send_data_frame(pivot_df.to_csv, "composition_local_vs_unipept.tsv", sep="\t")
