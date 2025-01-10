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

import base64
import io
import pandas as pd


# @app.callback(
#     Output('db_search_psm_table', 'children'),
#     Input('db_search_psm_upload', 'contents'),
#     Input('db_search_psm_upload', 'filename'),
#     Input('psm_table_selector', 'value')
# )
# def show_db_search_psm_table(content, names, sample):
#     """Display sample of db_search_psm data into datatable"""
#     if content is None:
#          return [html.H5("Import DB Search PSM...")]
#     elif sample is None and len(content) > 1:
#          return [html.H5("Select sample to display...")]
#     else:
#         sample_size = 500
        
#         content_idx = 0
#         if len(content) > 1:
#             content_idx = names.index(sample)      
        
#         df = import_csv_file(content[content_idx])
#         df = df[["Peptide", "-10lgP", "Mass", "Length", "ppm", "m/z", "Z", "RT", "Area", "Scan"]]
        
#         if df.shape[0] > sample_size:
#             df = df[:sample_size]
        
#         # limit peptide length for display
#         df['Peptide'] = df['Peptide'].apply(lambda x: x[:27] + "..." if len(x) > 30 else x)
        
#         return dash_table.DataTable(data=df.to_dict('records'),
#             columns=[{'id': c, 'name': c} for c in df.columns],
#             style_cell_conditional=[
#                 {'if': {'column_id': 'Peptide'},
#                 'width': '25%'},
#             ],
#             style_data={'table-layout': 'fixed'},
#             style_header={'backgroundColor': 'rgb(200, 200, 200)',
#                           'color': 'black',
#                           'fontWeight': 'bold'},
#             page_size=10)


# @app.callback(
#     Output('taxonomy_db_table', 'children'),
#     Input('taxonomy_db_upload', 'contents'),
#     Input('taxonomy_db_format_radio', 'value'),
# )
# def show_taxonomy_db_table(content, format):
#     """Display sample of protein db data into datatable
#     """
#     if content is None:
#          return [html.H5("Import Protein DB...")]
#     else:
#         if format == "accession-taxonomy map":
#             prot_obj = import_acc_tax_map(content)
#             df = prot_obj.accession_to_taxa
#         else:
#             return [html.H5("No fasta support...")]
#             df = import_fasta_file(content)
        
#         if df.shape[0] > 500:
#             df = df[:500]
        
#         return dash_table.DataTable(data=df.to_dict('records'),
#             columns=[{'id': c, 'name': c} for c in df.columns],
#             style_data={'table-layout': 'fixed'},
#             style_header={'backgroundColor': 'rgb(200, 200, 200)',
#                           'color': 'black',
#                           'fontWeight': 'bold'},
#             page_size=10)
        
        
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

# @app.callback(
#     Output('func_annot_table', 'children'),
#     Input('func_annot_db_upload', 'contents'),
#     Input('func_annot_db_format_radio', 'value'),
# )
# def show_functional_annot_table(content, format):
#     """Display sample of functional annotation data into datatable
#     """
#     if content is None:
#          return [html.H5("Import Functional Annotation Data...")]
#     else:
#         if format == 'eggnog':
#             df = import_func_map(content)
#             df = df.reset_index()
            
#             df.rename(columns={"Preferred_name": "Prot_name", "COG_category": "COG_cat"}, inplace=True)
            
#             # select specific columns
#             df = df[["query", "evalue", "Prot_name", "COG_cat", "EC", "KEGG_ko", "CAZy"]]
            
#             if df.shape[0] > 500:
#                 df = df[:500]
            
#             # limit lengths of certain columns that may contain much data
#             # df['seed_ortholog'] = df['seed_ortholog'].apply(lambda x: x[:25] + "..." if len(x) > 28 else x)
#             df['EC'] = df['EC'].astype(str).replace({'nan': ''})
#             df['EC'] = df['EC'].apply(lambda x: ",".join(x.split(",")[:2]) + ",..." if len(x.split(",")) > 2 else x)
#             df['KEGG_ko'] = df['KEGG_ko'].astype(str).replace({'nan': ''})
#             df['KEGG_ko'] = df['KEGG_ko'].apply(lambda x: ",".join(x.split(",")[:2]) + ",..." if len(x.split(",")) > 2 else x)
#             df['evalue'] = df['evalue'].apply(lambda x: "{res:.3E}".format(res=x))
            
            
#             return dash_table.DataTable(data=df.to_dict('records'),
#                 columns=[{'id': c, 'name': c} for c in df.columns],
#                 style_data={'table-layout': 'fixed'},
#                 style_header={'backgroundColor': 'rgb(200, 200, 200)',
#                             'color': 'black',
#                             'fontWeight': 'bold'},
#                 page_size=10)
#         else:
#             return [html.H5("No KEGG support...")]
#             df = import_fasta_file(content)
            
            
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
    
    
