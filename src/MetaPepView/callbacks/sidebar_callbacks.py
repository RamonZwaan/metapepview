from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from MetaPepView.server import app

from pathlib import Path

from MetaPepView.layout.sidebar import *
from backend import \
    check_ncbi_taxonomy_present, \
    check_kegg_mapping_present, \
    check_gtdb_taxonomy_present, \
    download_ncbi_taxonomy, \
    download_kegg_ko_map, \
    download_gtdb_taxonomy


import base64
import io
import pandas as pd
from textwrap import dedent


@app.callback(
    Output('database_present_status', 'data'),
    Output('ncbi_db_presence_check', 'children'),
    Output('gtdb_db_presence_check', 'children'),
    Output('kegg_map_presence_check', 'children'),
    Output('ncbi_taxonomy_db_loc', 'valid'),
    Output('gtdb_taxonomy_db_loc', 'valid'),
    Output('kegg_map_loc', 'valid'),
    Output('start_db_check', 'interval'),
    Input('start_db_check', 'n_intervals'),
    Input('data_imported_trigger', 'data'),
    Input('ncbi_taxonomy_db_loc', 'value'),
    Input('gtdb_taxonomy_db_loc', 'value'),
)
def validate_db_presence(n,
                         import_finish_trigger,
                         ncbi_path_change,
                         gtdb_path_change):
    db_status_dict = {
        "ncbi_taxonomy": False,
        "gtdb_taxonomy": False,
        "kegg_map": False
    }
    
    new_interval = 1e5
    success_icon = html.I(className="bi bi-check-circle-fill me-3 ms-3 fs-5 text-success")
    failed_icon = html.I(className="bi bi-x-circle-fill me-3 ms-3 fs-5 text-danger")
    
    if check_ncbi_taxonomy_present(Path(ncbi_path_change))[0] is True:
        ncbi_status = [success_icon, html.P("NCBI Taxonomy", className="fs-5")]
        ncbi_path_valid = True
        db_status_dict["ncbi_taxonomy"] = True
    else:
        new_interval = 1e5
        ncbi_status = [failed_icon, html.P("NCBI Taxonomy", className="fs-5")]
        ncbi_path_valid = False
    
    if check_gtdb_taxonomy_present(Path(gtdb_path_change))[0] is True:
        gtdb_status = [success_icon, html.P("GTDB Taxonomy", className="fs-5")]
        gtdb_path_valid = True
        db_status_dict["gtdb_taxonomy"] = True
    else:
        new_interval = 1e5
        gtdb_status = [failed_icon, html.P("GTDB Taxonomy", className="fs-5")]
        gtdb_path_valid = False
    
    if check_kegg_mapping_present() is True:
        kegg_status = [success_icon, html.P("KEGG Dataset", className="fs-5")]
        kegg_path_valid = True
        db_status_dict["kegg_map"] = True
    else:
        new_interval = 1e5
        kegg_status = [failed_icon, html.P("KEGG Dataset", className="fs-5")]
        kegg_path_valid = False

    return (
        db_status_dict,
        ncbi_status,
        gtdb_status,
        kegg_status,
        ncbi_path_valid,
        gtdb_path_valid,
        kegg_path_valid,
        new_interval
    )

@app.callback(
    Output('db_download_status_alert', 'children'),
    Output('db_download_status_alert', 'is_open'),
    Output('loader_fetch_db', 'children'),
    Output('data_imported_trigger', 'data'),
    Input('start_database_download', 'n_clicks'),
    State('fetch_ncbi_taxonomy_checkbox', 'value'),
    State('fetch_gtdb_taxonomy_checkbox', 'value'),
    State('fetch_kegg_map_checkbox', 'value'),
    State('ncbi_taxonomy_db_loc', 'value'),
    State('gtdb_taxonomy_db_loc', 'value'),
    State('ncbi_taxonomy_db_source_url', 'value'),
    State('gtdb_taxonomy_db_source_url', 'value'),
    State('ncbi_taxonomy_overwrite_old_checkbox', 'value'),
    State('gtdb_taxonomy_overwrite_old_checkbox', 'value'),
    State('kegg_map_overwrite_old_checkbox', 'value'),
    State('ncbi_taxonomy_create_parent_dirs_checkbox', 'value'),
    State('gtdb_taxonomy_create_parent_dirs_checkbox', 'value'),
    State('kegg_map_create_parent_dirs_checkbox', 'value'),
)
def fetch_databases(button_click,
                    ncbi_check,
                    gtdb_check,
                    kegg_check,
                    ncbi_loc,
                    gtdb_loc,
                    ncbi_source,
                    gtdb_source,
                    ncbi_overwrite,
                    gtdb_overwrite,
                    kegg_overwrite,
                    ncbi_create_dirs,
                    gtdb_create_dirs,
                    kegg_create_dirs):
    # format alert message
    alert_msg = [html.P("Failed to download the following:", className="mb-2")]
    fail_encountered = False

    # attempt to download datasets
    if ncbi_check is True:
        print("fetching ncbi...")
        ncbi_success, ncbi_msg = download_ncbi_taxonomy(ncbi_loc,
                                                        ncbi_source,
                                                        ncbi_overwrite,
                                                        ncbi_create_dirs)
        if ncbi_success is False:
            fail_encountered = True
            alert_msg+= [html.P(f"\nNCBI: {ncbi_msg}")]
        # attempt to download datasets
    if gtdb_check is True:
        print("fetching gtdb...")
        gtdb_success, gtdb_msg = download_gtdb_taxonomy(gtdb_loc,
                                                        gtdb_source,
                                                        gtdb_overwrite,
                                                        gtdb_create_dirs)
        if gtdb_success is False:
            fail_encountered = True
            alert_msg+= [html.P(f"\nGTDB: {gtdb_msg}")]
    if kegg_check is True:
        print("fetching kegg...")
        kegg_success, kegg_msg = download_kegg_ko_map(kegg_overwrite,
                                                      kegg_create_dirs)
        if kegg_success is False:
            fail_encountered = True
            alert_msg+= [html.P(f"\nKEGG: {kegg_msg}")]
    
    print("validate download success...")
    
    # return alert based on success
    if fail_encountered is True:
        return alert_msg, True, None, False
    else:
        return None, False, None, True

@app.callback(
    Output("fetch_database_modal", "is_open"),
    Input("fetch_db_modal_open", "n_clicks"), 
)
def toggle_fetch_db_modal(n1):
    if n1:
        return True

@app.callback(
    Output("ncbi_taxonomy_db_source_url", "disabled"),
    Output("ncbi_taxonomy_overwrite_old_checkbox", "disabled"),
    Output("ncbi_taxonomy_create_parent_dirs_checkbox", "disabled"),
    Input("fetch_ncbi_taxonomy_checkbox", "value")
)
def toggle_ncbi_tax_fetch_disabled(checkbox):
    return tuple([not checkbox] * 3) 

@app.callback(
    Output("gtdb_taxonomy_db_source_url", "disabled"),
    Output("gtdb_taxonomy_overwrite_old_checkbox", "disabled"),
    Output("gtdb_taxonomy_create_parent_dirs_checkbox", "disabled"),
    Input("fetch_gtdb_taxonomy_checkbox", "value")
)
def toggle_gtdb_tax_fetch_disabled(checkbox):
    return tuple([not checkbox] * 3) 

@app.callback(
    Output("kegg_map_overwrite_old_checkbox", "disabled"),
    Output("kegg_map_create_parent_dirs_checkbox", "disabled"),
    Input("fetch_kegg_map_checkbox", "value")
)
def toggle_kegg_ko_fetch_disabled(checkbox):
    return tuple([not checkbox] * 2)


@app.callback(
    Output("start_database_download", "disabled"),
    Input("fetch_kegg_map_checkbox", "value"),
    Input("fetch_gtdb_taxonomy_checkbox", "value"),
    Input("fetch_ncbi_taxonomy_checkbox", "value")
)
def toggle_fetch_data_button(checkbox_1, checkbox_2, checkbox_3):
    if checkbox_1 or checkbox_2 or checkbox_3:
        return False
    else:
        return True
