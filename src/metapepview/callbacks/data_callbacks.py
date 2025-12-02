from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from copy import deepcopy

from metapepview.server import app

from metapepview.backend import *
from metapepview.backend.utils import truncate_end
from metapepview.constants import GlobalConstants as gc


@app.callback(
    Output('experiment_sample_table', 'data'),
    Output('experiment_sample_table', 'columns'),
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
    table_cols = deepcopy(gc.experiment_sample_table_cols)

    # remove unrelevant columns from table depending on dashboard function level
    if gc.display_db_search is False:
        table_cols.remove("DB Search Imported")
        table_cols.remove("De Novo Imported")
        table_cols.remove("Taxonomy DB Name")
        table_cols.remove("Functional Annotation DB Name")
    if gc.display_de_novo is False:
        table_cols.remove("DB Search Imported")
        table_cols.remove("De Novo Imported")

    # display message to import data if no peptides dataset is present
    if peptides_json is None:
        return (None,
                [{'id': c, 'name': c} for c in table_cols],
                "-",
                "-",
                "-",
                "-")

    peptides_obj = MetaPepTable.read_json(peptides_json)
    peptides_df = peptides_obj.data

    db_search_format = peptides_obj.db_search_format
    de_novo_format = peptides_obj.de_novo_format
    tax_db_format = peptides_obj.taxonomy_db_format
    func_db_format = peptides_obj.functional_db_format

    # fetch sample names and annotation db name + formats
    sample_df = peptides_df.drop_duplicates(subset=['Sample Name'], keep="first")
    sample_df = sample_df[table_cols]

    # ensure that subset selection is dataframe, should not ever be triggered.
    if not isinstance(sample_df, pd.DataFrame):
        raise TypeError("DataFrame subset selection did not result in a DataFrame.")

    # limit length of name columns in dataset
    for col in ['Sample Name', 'Taxonomy DB Name', 'Functional Annotation DB Name']:
        def text_processing(cell) -> str:
            trunc = truncate_end(cell, 40)
            if not isinstance(trunc, str):
                trunc = "None"
            return trunc
        sample_df[col] = sample_df[col].apply(text_processing)

    return (sample_df.to_dict('records'),
            [{'id': c, 'name': c} for c in table_cols],
            db_search_format,
            de_novo_format,
            tax_db_format,
            func_db_format)


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
    Output('experiment_name', 'data', allow_duplicate=True),
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
    Output("experiment_name", "data", allow_duplicate=True),
    Output('peptides', 'data', allow_duplicate=True),
    Output("mzml_data", "data", allow_duplicate=True),
    Output("mzml_peaks_data", "data", allow_duplicate=True),
    Output("mzml_metadata", "data", allow_duplicate=True),
    Output("features_data", "data", allow_duplicate=True),
    Output("features_metadata", "data", allow_duplicate=True),
    Output("db_search_qa_data", "data", allow_duplicate=True),
    Output("de_novo_qa_data", "data", allow_duplicate=True),
    Output('clear_peptides_data', 'n_clicks'),
    Input('clear_peptides_data', 'n_clicks'),
    prevent_initial_call=True
)
def clear_peptide_data(n_clicks):
    if n_clicks > 0:
        return (None,)*10 + (0,)
    raise PreventUpdate


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
