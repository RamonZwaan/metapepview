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


@app.callback(
    Output('db_search_psm_qa_valid', 'data'),
    Output('db_search_psm_qa_name', 'children'),
    Output('db_search_psm_qa_upload', 'file_infos'),
    Output('db_search_qa_format_alert', 'children', allow_duplicate=True),
    Output('db_search_qa_format_alert', 'is_open', allow_duplicate=True),
    Output('db_search_psm_qa_import_box', 'style'),
    Input('db_search_psm_qa_upload', 'file_infos'),
    Input('mzml_metadata', 'data'),
    State('db_search_psm_qa_format', 'value'),
    prevent_initial_call=True
)
def show_db_psm_search_qa_name(contents, 
                               mzml_metadata,
                               file_format):
    """Display filename of annotated peptide dataset import.
    """
    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = upload_to_stringio(cont, archv)

        return validate_db_search(cont_buf, file_format)

        # mzml_metadata does not relate to currently loaded mzml for import
        # if mzml_metadata is None:
        #     return validate_db_search(cont_buf, file_format)
        # else:
        #     try:
        #         db_search_obj = load_metapep_db_search(cont_buf,
        #                                             name,
        #                                             file_format)
        #     except Exception as err:
        #         return False, err

        #     source_name = mzml_metadata["raw file name"] 
        #     if source_name not in db_search_obj.source_files:
        #         msg = "mzml and DB search not from same experiment."
        #         return False, msg
        # return True, None

    if contents == []:
        raise PreventUpdate

    path = Path(contents[0]["path"])
    filename = contents[0]["filename"]

    (valid_data, 
     name, 
     content, 
     err_msg, 
     success, 
     import_box_style) = validate_single_file(path, 
                                              filename, 
                                              valid_func, 
                                              drag_and_drop=True)
    
    # if validation failed, remove all data and set file_infos to []
    if content is None:
        for file in contents:
            Path(file["path"]).unlink(True)
        contents = []
    
    open_alert = not success

    return (valid_data, name, contents, err_msg, open_alert, import_box_style)



@app.callback(
    Output('denovo_qa_valid', 'data'),
    Output('denovo_qa_name', 'children'),
    Output('denovo_qa_upload', 'file_infos'),
    Output('de_novo_qa_format_alert', 'children'),
    Output('de_novo_qa_format_alert', 'is_open'),
    Output('denovo_qa_import_box', 'style'),
    Input('denovo_qa_upload', 'file_infos'),
    Input('denovo_qa_format', 'value'),
    Input("mzml_metadata", "data"))
def show_denovo_search_qa_name(contents, 
                               file_format, 
                               mzml_metadata):
    """Display filename of annotated peptide dataset import.
    """
    if contents == []:
        raise PreventUpdate

    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = upload_to_stringio(cont, archv)

        return validate_de_novo(cont_buf, file_format)
        
        # mzml_metadata does not relate to currently loaded mzml for import
        # if mzml_metadata is None:
        #     return validate_de_novo(cont_buf, file_format)
        # else:
        #     try:
        #         de_novo_obj = load_metapep_de_novo(cont_buf,
        #                                            name,
        #                                            file_format)
        #     except Exception as err:
        #         return False, err

        #     source_name = mzml_metadata["raw file name"] 
        #     if source_name not in de_novo_obj.source_files:
        #         msg = "mzml and de novo data not from same experiment."
        #         return False, msg
        # return True, None
    
    de_novo_path = Path(contents[0]["path"])
    filename = contents[0]["filename"]

    (valid_data, 
     name, 
     content, 
     err_msg, 
     success, 
     import_box_style) = validate_single_file(de_novo_path, 
                                              filename, 
                                              valid_func, 
                                              drag_and_drop=True)

    # if validation failed, remove all data and set file_infos to []
    if content is None:
        for file in contents:
            Path(file["path"]).unlink(True)
        contents = []

    open_alert = not success
    
    return (valid_data, name, contents, err_msg, open_alert, import_box_style)


@app.callback(
    Output("mzml_name", "children"),
    Input("mzml_upload", "file_infos")
)
def show_mzml(file_path):
    if (file_path != []) and (Path(file_path[0]["path"]).is_file()):
        file_name = file_path[0]["filename"]
        return truncate_end(file_name, 30)
    else:
        return "No file..."


@app.callback(
    Output("features_name", "children"),
    Input("features_upload", "file_infos")
)
def show_features(filename):
    if filename != []:
        return filename[0]["filename"]
    else:
        return "No file..."


@app.callback(
    Output('start_spectral_import_button', 'disabled'),
    Output('spectral_import_hint', 'children'),
    Input('mzml_upload', 'file_infos'),
    Input("db_search_psm_qa_valid", "data"),
    Input("denovo_qa_valid", "data"),
)
def inactivate_spectral_import_button(mzml_file_path,
                                      db_search_valid,
                                      de_novo_valid):
    tooltip_target = "start_spectral_import_button_wrapper"

    if len(mzml_file_path) == 0:
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
    State("mzml_upload", "file_infos"),
    State("features_upload", "file_infos"),
    State("db_search_psm_qa_upload", "file_infos"),
    State("db_search_psm_qa_format", "value"),
    State("db_search_psm_qa_valid", "data"),
    State("denovo_qa_upload", "file_infos"),
    State("denovo_qa_format", "value"),
    State("denovo_qa_valid", "data"),
    prevent_initial_call=True
)
def store_spectral_dataset(btn, 
                           mzml_file_path,
                           features_content, 
                           db_search_content,
                           db_search_format,
                           db_search_valid, 
                           de_novo_content,
                           de_novo_format,
                           de_novo_valid):
    loading_status = None

    # only update if valid data uploaded
    if (len(mzml_file_path) == 0) or (not Path(mzml_file_path[0]["path"]).is_file()):
        raise PreventUpdate    

    empty_data = (None,)*8

    # Import mzml data and return as compressed string data
    (mzml_data, 
     mzml_peaks, 
     mzml_metadata,
     mzml_valid) = import_mzml(Path(mzml_file_path[0]["path"]), mzml_file_path[0]["filename"])
    
    #TODO: Add check that DB search and de novo is related to raw mzml file

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
    if features_content != []:
        (feature_data, 
        feature_metadata, 
        feature_valid) = import_features(Path(features_content[0]["path"]), 
                                         features_content[0]["filename"], 
                                         mzml_metadata)
    else:
        # no feature data is considered valid to prevent raising alert message
        feature_data = None
        feature_metadata = None
        feature_valid = True

    # compress db search and de novo data
    if db_search_content != []:
        db_search_path = Path(db_search_content[0]["path"])
        db_search_archive_format = determine_archive_format(db_search_content[0]["filename"])
        if db_search_archive_format is not None:
            db_search_files, db_search_names = archive_to_file_list(db_search_path,
                                                                    db_search_archive_format)
            db_search_path = db_search_files[0]

        db_search_data = load_metapep_db_search(db_search_path, 
                                                "sample", 
                                                db_search_format)
        
        if mzml_metadata["raw file name"] not in db_search_data.source_files:
            msg = "mzml and DB search not from same experiment."
            return empty_data + (spectral_data_import_container, msg, True)

        db_search_json = db_search_data.to_json()
        db_search_store = compress_string(db_search_json)
    else:
        db_search_store = None
    if de_novo_content != []:
        de_novo_path = Path(de_novo_content[0]["path"])
        de_novo_archive_format = determine_archive_format(de_novo_content[0]["filename"])
        if de_novo_archive_format is not None:
            de_novo_files, de_novo_names = archive_to_file_list(de_novo_path,
                                                                de_novo_archive_format)
            de_novo_path = de_novo_files[0]

        de_novo_data = load_metapep_de_novo(de_novo_path,
                                            de_novo_content[0]["filename"],
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

    # store mzml data on server, store key to data in store
    mzml_peaks_stored = add_dataset_to_server_store(app, 
                                                    "mzml_peaks_data", 
                                                    mzml_peaks)
    mzml_data_stored = add_dataset_to_server_store(app, 
                                                   "mzml_data", 
                                                   mzml_data)

    return (loading_status,
            mzml_data_stored,
            mzml_peaks_stored,
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
        remove_dataset_from_server_store(app, "mzml_peaks_data")
        remove_dataset_from_server_store(app, "mzml_data")
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


    return (truncate_end(spectra_name, 30),
            failed_icon if features is None else success_icon,
            failed_icon if db_search is None else success_icon,
            db_search_format,
            failed_icon if de_novo is None else success_icon,
            de_novo_format)
