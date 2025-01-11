from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
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
from backend.utils import truncate_end
        
        
@app.callback(
    Output('sample_table', 'children'),
    Input('peptides', 'data'),
)
def show_samples_table(peptides_json):
    """Display sample of protein db data into datatable
    """
    # display message to import data if no peptides dataset is present
    if peptides_json is None:
        return [html.H4("Import peptides dataset or annotate manually...",
                        className="px-2 py-2")]
    
    peptides_df = MetaPepTable.read_json(peptides_json).data

    # fetch sample names and annotation db name + formats
    sample_df = peptides_df.drop_duplicates(subset=['Sample Name'], keep="first")
    sample_df = sample_df[['Sample Name',
                           'DB Search Imported',
                           'De Novo Imported',
                           'Taxonomy DB Name',
                           'Functional Annotation DB Name']]
    
    # limit length of name columns in dataset
    for col in ['Sample Name', 'Taxonomy DB Name', 'Functional Annotation DB Name']:
        def text_processing(cell) -> str:
            trunc = truncate_end(cell, 40)
            if not isinstance(trunc, str):
                trunc = "None"
            return trunc
        sample_df[col] = sample_df[col].apply(text_processing)
    
    return dash_table.DataTable(data=sample_df.to_dict('records'),
        columns=[{'id': c, 'name': c} for c in sample_df.columns],
        style_data={'table-layout': 'fixed'},
        style_header={'backgroundColor': 'rgb(210, 210, 210)',
                        'color': 'black',
                        'fontWeight': 'bold'},
        style_cell={'textAlign': 'left', 'font-family': 'Arial'},
        style_cell_conditional=[
            {'if': {'column_id': 'DB Search Imported'},
            'width': '12%'},
            {'if': {'column_id': 'De Novo Imported'},
            'width': '12%'}],
        style_as_list_view=True,
        page_size=10)

            
@app.callback(
    Output('psm_table_selector', 'options'),
    Input('db_search_psm_upload', 'filename'),
)
def update_psm_sample_menu(names):
    """Update dropdown menu with imported sample names
    """
    print(names)
    if names is None:
        return []
    else:
        return names
    
    
