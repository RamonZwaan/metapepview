from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from metapepview.server import app
# import layout elements
from metapepview.layout.annotation_page import *
from metapepview.layout.sidebar import *

from metapepview.backend.annotation import *
from metapepview.backend.exceptions import AnnotationError
from metapepview.backend.io import *
from metapepview.backend.types import *
from metapepview.backend.type_operations import *
from metapepview.backend.utils import *

from metapepview.constants import GlobalConstants as gc



@app.callback(
    Output('db_search_psm_qa_valid', 'data'),
    Output('db_search_psm_qa_name', 'children'),
    Output('db_search_psm_qa_upload', 'contents'),
    Output('db_search_qa_format_alert', 'children', allow_duplicate=True),
    Output('db_search_qa_format_alert', 'is_open', allow_duplicate=True),
    Output('db_search_psm_qa_import_box', 'style'),
    Input('db_search_psm_qa_upload', 'contents'),
    Input('mzml_metadata', 'data'),
    State('db_search_psm_qa_format', 'value'),
    State('db_search_psm_qa_upload', 'filename'),
    State('db_search_psm_qa_upload', 'last_modified'),
    prevent_initial_call=True
)
def show_db_psm_search_qa_name(contents, 
                               mzml_metadata,
                               file_format,
                               name, 
                               date):
    """Display filename of annotated peptide dataset import.
    """
    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)

        if mzml_metadata is None:
            return validate_db_search(cont_buf, file_format)
        else:
            try:
                db_search_obj = load_metapep_db_search(cont_buf,
                                                    name,
                                                    file_format)
            except Exception as err:
                return False, err

            source_name = mzml_metadata["raw file name"] 
            if source_name not in db_search_obj.source_files:
                msg = "mzml and DB search not from same experiment."
                return False, msg
        return True, None
    
    (valid_data, 
     name, 
     content, 
     err_msg, 
     success, 
     import_box_style) = validate_single_file(contents, 
                                              name, 
                                              date, 
                                              valid_func, 
                                              drag_and_drop=True)
    
    open_alert = not success

    return (valid_data, name, content, err_msg, open_alert, import_box_style)



@app.callback(
    Output('denovo_qa_valid', 'data'),
    Output('denovo_qa_name', 'children'),
    Output('denovo_qa_upload', 'contents'),
    Output('de_novo_qa_format_alert', 'children'),
    Output('de_novo_qa_format_alert', 'is_open'),
    Output('denovo_qa_import_box', 'style'),
    Input('denovo_qa_upload', 'contents'),
    Input('denovo_qa_format', 'value'),
    Input("mzml_metadata", "data"),
    State('denovo_qa_upload', 'filename'),
    State('denovo_qa_upload', 'last_modified'))
def show_denovo_search_qa_name(contents, 
                               file_format, 
                               mzml_metadata,
                               name, 
                               date):
    """Display filename of annotated peptide dataset import.
    """
    if contents is None:
        raise PreventUpdate

    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)

        if mzml_metadata is None:
            return validate_de_novo(cont_buf, file_format)
        else:
            try:
                de_novo_obj = load_metapep_de_novo(cont_buf,
                                                   name,
                                                   file_format)
            except Exception as err:
                return False, err

            source_name = mzml_metadata["raw file name"] 
            if source_name not in de_novo_obj.source_files:
                msg = "mzml and de novo data not from same experiment."
                return False, msg
        return True, None
    

    (valid_data, 
     name, 
     content, 
     err_msg, 
     success, 
     import_box_style) = validate_single_file(contents, 
                                              name, 
                                              date, 
                                              valid_func, 
                                              drag_and_drop=True)
    
    open_alert = not success
    
    return (valid_data, name, content, err_msg, open_alert, import_box_style)


@app.callback(
    Output("mzml_name", "children"),
    Input("mzml_upload", "filename")
)
def show_mzml(filename):
    return filename


@app.callback(
    Output("features_name", "children"),
    Input("features_upload", "filename")
)
def show_features(filename):
    return filename


@app.callback(
    Output('start_spectral_import_button', 'disabled'),
    Output('spectral_import_hint', 'children'),
    Input('mzml_upload', 'contents'),
    Input("db_search_psm_qa_valid", "data"),
    Input("denovo_qa_valid", "data"),
)
def inactivate_spectral_import_button(mzml_content,
                                      db_search_valid,
                                      de_novo_valid):
    tooltip_target = "start_spectral_import_button_wrapper"

    if mzml_content is None:
        tooltip = dbc.Tooltip("Import spectral dataset (mzML)",
                              target=tooltip_target,
                              placement="bottom",
                              className="mt-1")
        return (True, [tooltip])
    elif db_search_valid is False:
        tooltip = dbc.Tooltip("Invalid db search file supplied. Check format or if it is from same experiment as mzml",
                              target=tooltip_target,
                              placement="bottom",
                              className="mt-1")
        return (True, [tooltip])
    elif de_novo_valid is False:
        tooltip = dbc.Tooltip("Invalid de novo file supplied. Check format or if it is from same experiment as mzml",
                              target=tooltip_target,
                              placement="bottom",
                              className="mt-1")
        return (True, [tooltip])

    else:
        tooltip = dbc.Tooltip(f"start import",
                                target="start_spectral_import_button_wrapper",
                                placement="bottom",
                                className="mt-1")
        return (False, [tooltip])


@app.callback(
    Output('spectral_import_loading_spot', 'children'),
    Output("mzml_data", "data", allow_duplicate=True),
    Output("mzml_peaks_data", "data", allow_duplicate=True),
    Output("mzml_metadata", "data", allow_duplicate=True),
    Output("features_data", "data", allow_duplicate=True),
    Output("features_metadata", "data", allow_duplicate=True),
    Output("db_search_qa_data", "data", allow_duplicate=True),
    Output("de_novo_qa_data", "data", allow_duplicate=True),
    Output("spectral_data_import_container", "children"),
    Output("qa_data_import_alert", "children", allow_duplicate=True),
    Output("qa_data_import_alert", "is_open", allow_duplicate=True),
    Input("start_spectral_import_button", "n_clicks"),
    State("mzml_upload", "contents"),
    State("mzml_upload", "filename"),
    State("features_upload", "contents"),
    State("features_upload", "filename"),
    State("db_search_psm_qa_upload", "contents"),
    State("db_search_psm_qa_format", "value"),
    State("db_search_psm_qa_valid", "data"),
    State("denovo_qa_upload", "contents"),
    State("denovo_qa_format", "value"),
    State("denovo_qa_valid", "data"),
    prevent_initial_call=True
)
def store_spectral_dataset(btn, 
                           mzml_content, 
                           mzml_filename,
                           features_content, 
                           features_filename,
                           db_search_content,
                           db_search_format,
                           db_search_valid, 
                           de_novo_content,
                           de_novo_format,
                           de_novo_valid):
    loading_status = None

    # only update if valid data uploaded
    if mzml_content is None:
        raise PreventUpdate    

    empty_data = (None,)*8

    # Import mzml data and return as compressed string data
    (mzml_data, 
     mzml_peaks, 
     mzml_metadata,
     mzml_valid) = import_mzml(mzml_content, mzml_filename)
    
    # require mzml data for any information to be stored
    if mzml_valid is False or any(x is False for x in [db_search_valid, de_novo_valid]):
        if mzml_valid is False:
            msg = "Failed to load mzml dataset"
        elif db_search_valid is False:
            msg = "Invalid db search file... Potentially due to invalid format or different experiment source from mzml"
        elif de_novo_valid is False:
            msg = "Invalid de novo file... Potentially due to invalid format or different experiment source from mzml"
        return empty_data + (spectral_data_import_container, msg, True)

    # import feature data and return as compressed string data
    if features_content is not None:
        (feature_data, 
        feature_metadata, 
        feature_valid) = import_features(features_content, 
                                         features_filename, 
                                         mzml_metadata)
    else:
        # no feature data is considered valid to prevent raising alert message
        feature_data = None
        feature_metadata = None
        feature_valid = True

    # compress db search and de novo data
    if db_search_content is not None:
        db_search_data = load_metapep_db_search(db_search_content, 
                                                "sample", 
                                                db_search_format)
        
        if mzml_metadata["raw file name"] not in db_search_data.source_files:
            msg = "mzml and DB search not from same experiment."
            return empty_data + (spectral_data_import_container, msg, True)

        db_search_json = db_search_data.to_json()
        db_search_store = compress_string(db_search_json)
    else:
        db_search_store = None
    if de_novo_content is not None:
        de_novo_data = load_metapep_de_novo(de_novo_content,
                                            "sample",
                                            de_novo_format)

        if mzml_metadata["raw file name"] not in de_novo_data.source_files:
            msg = "mzml and de novo not from same experiment."
            return empty_data + (spectral_data_import_container, msg, True)

        de_novo_json = de_novo_data.to_json()
        de_novo_store = compress_string(de_novo_json)
    else:
        de_novo_store = None

    if feature_valid is False:
        alert_open = True
        alert_msg = "Failed to load Feature dataset, spectral data imported without features."
    else:
        alert_open = False
        alert_msg = None

    return (loading_status,
            mzml_data,
            mzml_peaks,
            mzml_metadata,
            feature_data,
            feature_metadata,
            db_search_store,
            de_novo_store,
            spectral_data_import_container,
            alert_msg,
            alert_open)


@app.callback(
    Output("mzml_data", "data", allow_duplicate=True),
    Output("mzml_peaks_data", "data", allow_duplicate=True),
    Output("mzml_metadata", "data", allow_duplicate=True),
    Output("features_data", "data", allow_duplicate=True),
    Output("features_metadata", "data", allow_duplicate=True),
    Output("db_search_qa_data", "data", allow_duplicate=True),
    Output("de_novo_qa_data", "data", allow_duplicate=True),
    # Output("clear_spectral_dataset", "n_clicks"),
    Input("clear_spectral_dataset", "n_clicks"),
    prevent_initial_call=True,
)
def clear_spectral_data(n_clicks):
    if n_clicks > 0:
        return (None,)*7# + (0,)
    raise PreventUpdate


@app.callback(
    Output("mzml_store_name", "children"),
    Output("feature_store_valid", "className"),
    Output("db_search_qa_store_valid", "className"),
    Output("db_search_qa_store_format", "children"),
    Output("de_novo_qa_store_valid", "className"),
    Output("de_novo_qa_store_name", "children"),
    Input("mzml_metadata", "data"),
    Input("features_data", "data"),
    Input("db_search_qa_data", "data"),
    Input("de_novo_qa_data", "data")
)
def show_imported_spectra(mzml_metadata,
                          features,
                          db_search,
                          de_novo):
    # icon classname format
    icon_classname = "bi me-3 ms-3 fs-5 "
    failed_icon = icon_classname + "bi-x-circle-fill me-3 ms-3 fs-5 text-danger"
    success_icon = icon_classname + "bi-check-circle-fill me-3 ms-3 fs-5 text-success"
    
    if mzml_metadata is not None:
        spectra_name = mzml_metadata["raw file name"]
    else:
        spectra_name = "-"

    if db_search is not None:
        db_search_obj = MetaPepDbSearch.read_json(
            decompress_string(db_search)
        )
        db_search_format = db_search_obj.data_source
    else:
        db_search_format = "-"

    if de_novo is not None:
        de_novo_obj = MetaPepDeNovo.read_json(
            decompress_string(de_novo)
        )
        de_novo_format = de_novo_obj.data_source
    else:
        de_novo_format = "-"


    return (spectra_name,
            failed_icon if features is None else success_icon,
            failed_icon if db_search is None else success_icon,
            db_search_format,
            failed_icon if de_novo is None else success_icon,
            de_novo_format)
