from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

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
from backend.utils import truncate_end
        
        
@app.callback(
    Output('experiment_sample_table', 'data'),
    Output('peptides_db_search_format', 'children'),
    Output('peptides_de_novo_format', 'children'),
    Output('peptides_taxonomy_db_format', 'children'),
    Output('peptides_function_db_format', 'children'),
    Input('peptides', 'data')
)
def show_samples_data(peptides_json):
    """Display sample of protein db data into datatable and display format
    information.
    """
    
    table_cols = GlobalConstants.experiment_sample_table_cols
    
    # display message to import data if no peptides dataset is present
    if peptides_json is None:
        return (None,
                "-",
                "-",
                "-",
                "-")
        # return ([html.H4("Import peptides dataset or annotate manually...",
        #                  className="px-2 py-2")],
        #         "-",
        #         "-",
        #         "-",
        #         "-")
    
    peptides_obj = MetaPepTable.read_json(peptides_json)
    peptides_df = peptides_obj.data

    db_search_format = peptides_obj.db_search_format
    de_novo_format = peptides_obj.de_novo_format
    tax_db_format = peptides_obj.taxonomy_db_format
    func_db_format = peptides_obj.functional_db_format

    # fetch sample names and annotation db name + formats
    sample_df = peptides_df.drop_duplicates(subset=['Sample Name'], keep="first")
    sample_df = sample_df[table_cols]
    
    # limit length of name columns in dataset
    for col in ['Sample Name', 'Taxonomy DB Name', 'Functional Annotation DB Name']:
        def text_processing(cell) -> str:
            trunc = truncate_end(cell, 40)
            if not isinstance(trunc, str):
                trunc = "None"
            return trunc
        sample_df[col] = sample_df[col].apply(text_processing)
    
    
    return (sample_df.to_dict('records'),
            db_search_format, 
            de_novo_format, 
            tax_db_format, 
            func_db_format)
    
    # return (dash_table.DataTable(data=sample_df.to_dict('records'),
    #             columns=[{'id': c, 'name': c} for c in sample_df.columns],
    #             style_data={'table-layout': 'fixed'},
    #             style_header={'backgroundColor': 'rgb(210, 210, 210)',
    #                             'color': 'black',
    #                             'fontWeight': 'bold'},
    #             style_cell={'textAlign': 'left', 'font-family': 'Arial'},
    #             style_cell_conditional=[
    #                 {'if': {'column_id': 'DB Search Imported'},
    #                 'width': '12%'},
    #                 {'if': {'column_id': 'De Novo Imported'},
    #                 'width': '12%'}],
    #             style_as_list_view=True,
    #             row_deletable=True,
    #             page_size=10), 
    #         db_search_format, 
    #         de_novo_format, 
    #         tax_db_format, 
    #         func_db_format)


@app.callback(
    Output('peptides', 'data', allow_duplicate=True),
    Input('experiment_sample_table', 'data'),
    State('peptides', 'data'),
    prevent_initial_call=True
)
def remove_peptide_data(datatable_data,
                        peptides_json):
    if peptides_json is None:
        raise PreventUpdate
        
    metapep_obj = MetaPepTable.read_json(peptides_json)
    metapep_samples = metapep_obj.sample_names
    
    samples_datatable = [row.get('Sample Name') for row in datatable_data]
    filter_samples = set(metapep_samples) - set(samples_datatable)
    if len(filter_samples) == 0:
        raise PreventUpdate
    
    metapep_obj = metapep_obj.remove_samples(filter_samples)
    
    return metapep_obj.to_json()


@app.callback(
    Output("experiment_name_field", "value", allow_duplicate=True),
    Output('experiment_name', 'data'),
    Input("experiment_name_field", "value"),
    State('experiment_name', 'data'),
    prevent_initial_call=True
)
def update_experiment_name(current_field_name, stored_data):
    if current_field_name is None:
        return stored_data, stored_data
    else:
        return current_field_name, current_field_name


@app.callback(
    Output("experiment_name_field", "value", allow_duplicate=True),
    Output('peptides', 'data', allow_duplicate=True),
    Input('clear_peptides_data', 'n_clicks'),
    prevent_initial_call=True
)
def clear_peptide_data(_):
    return None , None
 
            
@app.callback(
    Output('psm_table_selector', 'options'),
    Input('db_search_psm_upload', 'filename'),
)
def update_psm_sample_menu(names):
    """Update dropdown menu with imported sample names
    """
    if names is None:
        return []
    else:
        return names
    
    
