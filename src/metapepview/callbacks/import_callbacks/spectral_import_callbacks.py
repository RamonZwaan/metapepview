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



# @app.callback(
#     Output('db_search_psm_qa_valid', 'data'),
#     Output('db_search_psm_qa_name', 'children'),
#     Output('db_search_psm_qa_upload', 'contents'),
#     # Output('db_search_psm_qa_format_alert', 'children'),
#     # Output('db_search_psm_qa_format_alert', 'is_open'),
#     Output('db_search_psm_qa_import_box', 'style'),
#     Input('db_search_psm_qa_upload', 'contents'),
#     Input('db_search_psm_qa_format', 'value'),
#     State('db_search_psm_qa_upload', 'filename'),
#     State('db_search_psm_qa_upload', 'last_modified'))
# def show_db_psm_search_qa_name(contents, file_format, name, date):
#     """Display filename of annotated peptide dataset import.
#     """
#     # set validation function
#     def valid_func(cont, archv) -> Tuple[bool, str | None]:
#         cont_buf = memory_to_stringio(cont, archv)
#         return validate_db_search(cont_buf, file_format)
    
#     valid_data, name, content, err_msg, success, import_box_style = validate_single_file(contents, name, date, valid_func, drag_and_drop=True)
    
#     # update import box style
#     qa_box_style = deepcopy(StyleConstants.qa_import_box_style)
#     if "background-color" in import_box_style.keys():
#         qa_box_style["background-color"] = import_box_style["background-color"]
    
#     return (valid_data, name, content, qa_box_style)


# @app.callback(
#     Output('denovo_qa_valid', 'data'),
#     Output('denovo_qa_name', 'children'),
#     Output('denovo_qa_upload', 'contents'),
#     # Output('denovo_qa_format_alert', 'children'),
#     # Output('denovo_qa_format_alert', 'is_open'),
#     Output('denovo_qa_import_box', 'style'),
#     Input('denovo_qa_upload', 'contents'),
#     Input('denovo_qa_format', 'value'),
#     State('denovo_qa_upload', 'filename'),
#     State('denovo_qa_upload', 'last_modified'))
# def show_denovo_search_qa_name(contents, file_format, name, date):
#     """Display filename of annotated peptide dataset import.
#     """
#     # set validation function
#     def valid_func(cont, archv) -> Tuple[bool, str | None]:
#         cont_buf = memory_to_stringio(cont, archv)
#         return validate_de_novo(cont_buf, file_format)
#     valid_data, name, content, err_msg, success, import_box_style = validate_single_file(contents, name, date, valid_func, drag_and_drop=True)

#     # update import box style
#     qa_box_style = deepcopy(StyleConstants.qa_import_box_style)
#     if "background-color" in import_box_style.keys():
#         qa_box_style["background-color"] = import_box_style["background-color"]
    
#     return (valid_data, name, content, qa_box_style)


# @app.callback(
#     Output('start_annotation_button', 'disabled'),
#     Output('annotation_hint', 'children'),
#     Input('mzml_valid', 'data'),

# )
# def inactivate_spectral_import_button(mzml_valid):
#     tooltip_target = "start_spectral_import_button_wrapper"

#     if mzml_valid is not True:
#         tooltip = dbc.Tooltip("Import spectral dataset (mzML)",
#                               target=tooltip_target,
#                               placement="bottom",
#                               className="mt-1")
#         return (True, [tooltip])

#     else:
#         tooltip = dbc.Tooltip(f"start annotation",
#                                 target="start_annotation_button_wrapper",
#                                 placement="bottom",
#                                 className="mt-1")
#     return (False, [tooltip])


@app.callback(
    Output("mzml_data", "data"),
    Output("mzml_peaks_data", "data"),
    Output("mzml_metadata", "data"),
    Output("mzml_upload", "contents"),
    Output("mzml_name", "children"),
    Output('mzml_valid', 'data'),
    Output('mzml_import_box', 'style'),
    Output("qa_data_import_alert", "children", allow_duplicate=True),
    Output("qa_data_import_alert", "is_open", allow_duplicate=True),
    Input("start_spectral_import_button", "n_clicks"),
    State("mzml_upload", "contents"),
    State("mzml_upload", "filename"),
    prevent_initial_call=True
)
def store_mzml_dataset(content, filename):
    # only update if valid data uploaded
    if content is None:
        raise PreventUpdate
    print("Start mzML wrangling...")
    
    fields = [
        "scan number",
        "MS level",
        "peaks count",
        "retention time",
        "total ion current",
        "precursor intensity",
        "precursor scan number",
        "ion injection time",
        "precursor m/z",
        "m/z array",
        "intensity array"
    ]
    
    archive_format = determine_archive_format(filename)
    
    qa_box_style = deepcopy(StyleConstants.qa_import_box_style)
    
    try:
        data, metadata = mzml_to_df(memory_to_file_like(content, archive_format),
                                    fields)
        
        numeric_fields = [
            'scan number',
            'MS level',
            'peaks count',
            'retention time',
            'total ion current',
            'precursor intensity',
            'precursor m/z',
            'precursor scan number',
            'ion injection time'
        ]
        
        data[numeric_fields] = data[numeric_fields].astype(float)
        
        # store peaks inside separate dataset
        peaks_data = data[["m/z array", "intensity array"]]
        peaks_data = peaks_data.to_json(orient="index")
        data = data.drop(labels=["m/z array", "intensity array"], axis=1)
        
        metadata["total retention time"] = data.iloc[-1]['retention time']    
    except Exception as err:
        print(err)
        
        qa_box_style["background-color"] = StyleConstants.import_failed_color
        return (
            None, 
            None, 
            None, 
            None, 
            import_single_file(None, None, drag_and_drop=True),
            False,
            qa_box_style,
            f"Failed to import mzML data: {err}",
            True
        )
    print("Finished wrangling...")
    # After dataset is imported, remove upload data to save memory
    # However, do keep the filename for display in dashboard
    
    qa_box_style["background-color"] = StyleConstants.import_success_color
    return (
        compress_string(data.to_json()), 
        compress_string(peaks_data), 
        metadata, 
        None, 
        import_single_file(filename, None, drag_and_drop=True),
        True,
        qa_box_style,
        None,
        False
    )


@app.callback(
    Output("features_data", "data"),
    Output("features_metadata", "data"),
    Output("features_upload", "contents"),
    Output("features_name", "children"),
    Output('features_valid', 'data'),
    Output('features_import_box', 'style'),
    Output("qa_data_import_alert", "children", allow_duplicate=True),
    Output("qa_data_import_alert", "is_open", allow_duplicate=True),
    Input("features_upload", "contents"),
    Input("features_upload", "filename"),
    State("mzml_metadata", "data"),
    prevent_initial_call=True
)
def store_features_dataset(content, filename, mzml_metadata):
    # only update if valid data uploaded
    if content is None:
        raise PreventUpdate
    qa_box_style = deepcopy(StyleConstants.qa_import_box_style)
    
    # only process features file after mzml file is imported
    if mzml_metadata is None:
        qa_box_style["background-color"] = StyleConstants.import_failed_color
        return (
            None, 
            None, 
            None, 
            import_single_file(None, None, drag_and_drop=True),
            False,
            qa_box_style,
            "Need to import mzml data first...",
            True
        )
    
    archive_format = determine_archive_format(filename)
    
    print("Process feature data...")
    try:
        data, metadata = featurexml_to_df(memory_to_file_like(content, archive_format), 
                                          None)
    except Exception as err:
        print(err)
        
        qa_box_style["background-color"] = StyleConstants.import_failed_color
        return (
            None, 
            None, 
            None, 
            import_single_file(None, None, drag_and_drop=True),
            False,
            qa_box_style,
            f"Failed to import features: {err}",
            True
        )
    print("Finished feature processing...")
    # After dataset is imported, remove upload data to save memory
    # However, do keep the filename for display in dashboard
    
    qa_box_style["background-color"] = StyleConstants.import_success_color
    return (
        compress_string(data.to_json()), 
        metadata, 
        None, 
        import_single_file(filename, None, drag_and_drop=True),
        True,
        qa_box_style,
        None,
        False
    )

