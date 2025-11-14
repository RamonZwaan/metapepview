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
    Output("db_search_modal", "is_open"),
    Input("db_search_modal_open", "n_clicks"),
)
def toggle_db_search_filters(n1):
    if n1:
        return True

@app.callback(
    Output("de_novo_modal", "is_open"),
    Input("de_novo_modal_open", "n_clicks"),
)
def toggle_de_novo_filters(n1):
    if n1:
        return True

@app.callback(
    Output("taxonomy_map_modal", "is_open"),
    Input("taxonomy_map_modal_open", "n_clicks"),
)
def toggle_tax_map_filters(n1):
    if n1:
        return True

@app.callback(
    Output("function_map_modal", "is_open"),
    Input("function_map_modal_open", "n_clicks"),
)
def toggle_func_map_filters(n1):
    if n1:
        return True


@app.callback(
    Output('peptides_name', 'children'),
    Input('peptides_upload', 'contents'),
    State('peptides_upload', 'filename'),
    State('peptides_upload', 'last_modified'))
def show_peptides_names(contents, names, dates):
    """Display filename of annotated peptide dataset import.
    """
    return import_single_file(names, dates, max_name_len=30, drag_and_drop=False)



@app.callback(
    Output('db_search_psm_valid', 'data'),
    Output('db_search_psm_name', 'children', allow_duplicate=True),
    Output("sample_name_import", "value"),
    Output('db_search_psm_upload', 'contents', allow_duplicate=True),
    Output('db_search_format_alert', 'children', allow_duplicate=True),
    Output('db_search_format_alert', 'is_open', allow_duplicate=True),
    Output('db_search_import_box', 'style'),
    Input('db_search_psm_upload', 'contents'),
    State('db_search_psm_upload', 'filename'),
    State('db_search_psm_format', 'value'),
    State('db_search_psm_upload', 'last_modified'),
    State("sample_name_import", "value"),
    prevent_initial_call=True)
def show_db_search_psm_names(contents, names, format, dates, current_sample_name):
    """Display filename of db search psm import
    """
        # set up validator, including checking of compression type
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)
        return validate_db_search(cont_buf, format)
    valid_data, name_list, contents, msg, success_status, box_style = validate_multiple_files(
        contents,
        names,
        dates,
        valid_func)

    # if no sample name given, give it the first item from PSM files
    if current_sample_name is None and valid_data is True:
        current_sample_name = names[0]
    return (valid_data, name_list, current_sample_name, contents, msg, False,
            box_style)


@app.callback(
    Output('db_search_psm_name', 'children', allow_duplicate=True),
    Output('db_search_psm_upload', 'contents', allow_duplicate=True),
    Input('db_search_psm_format', 'value'),
    prevent_initial_call=True
)
def reset_db_search_import_box(new_format):
    return None, None


@app.callback(
    Output('denovo_name', 'children', allow_duplicate=True),
    Output('denovo_upload', 'contents', allow_duplicate=True),
    Input('de_novo_format', 'value'),
    prevent_initial_call=True
)
def reset_de_novo_import_box(new_format):
    return None, None


@app.callback(
    Output('de_novo_valid', 'data'),
    Output('denovo_name', 'children'),
    Output('denovo_upload', 'contents'),
    Output('de_novo_format_alert', 'children', allow_duplicate=True),
    Output('de_novo_format_alert', 'is_open', allow_duplicate=True),
    Output('de_novo_import_box', 'style'),
    Input('denovo_upload', 'contents'),
    State('denovo_upload', 'filename'),
    State('de_novo_format', 'value'),
    State('denovo_upload', 'last_modified'),
    prevent_initial_call=True)
def show_denovo_names(contents, names, de_novo_format, dates):
    """Display filename of denovo dataset import.
    """
    # set up validator, including checking of compression type
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)
        return validate_de_novo(cont_buf, de_novo_format)

    return validate_multiple_files(contents,
                                   names,
                                   dates,
                                   valid_func)



@app.callback(
    Output('func_annot_db_valid', 'data'),
    Output('func_annot_name', 'children'),
    Output('func_annot_db_upload', 'contents', allow_duplicate=True),
    Output('functional_db_format_alert', 'children'),
    Output('functional_db_format_alert', 'is_open'),
    Output('functional_db_import_box', 'style'),
    Input('func_annot_db_upload', 'contents'),
    State('func_annot_db_upload', 'filename'),
    State('func_annot_db_format', 'value'),
    State('func_annot_db_upload', 'last_modified'),
    prevent_initial_call=True)
def show_functional_db_name(contents, name, func_db_format, dates):
    """Display filename of protein database import.
    """
    # set up validator, including checking of compression type
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)
        return validate_func_map(cont_buf, func_db_format)

    return validate_single_file(
        contents,
        name,
        dates,
        valid_func,
    )


@app.callback(
    Output('taxonomy_db_valid', 'data'),
    Output('taxonomy_db_name', 'children'),
    Output('taxonomy_db_upload', 'contents'),
    Output('taxonomy_db_format_alert', 'children'),
    Output('taxonomy_db_format_alert', 'is_open'),
    Output('taxonomy_db_import_box', 'style'),
    Input('taxonomy_db_upload', 'contents'),
    Input('acc_tax_map_delim', 'value'),
    Input('acc_tax_map_acc_idx', 'value'),
    Input('acc_tax_tax_idx', 'value'),
    Input('taxonomy_db_format', 'value'),
    Input('taxonomy_id_format_checkbox', 'value'),
    State('taxonomy_db_upload', 'filename'),
    State('taxonomy_db_upload', 'last_modified'),
    prevent_initial_call=True)
def show_taxonomy_db_name(contents,
                          delim,
                          acc_idx,
                          tax_idx,
                          tax_format,
                          tax_element_format,
                          name,
                          date):
    """Display filename of functional annotation import
    """
    valid_func = lambda cont, archv: validate_acc_tax_map(cont,
                                                          acc_idx,
                                                          tax_idx,
                                                          delim,
                                                          tax_format,
                                                          tax_element_format,
                                                          archv)
    return validate_single_file(contents,
                                name,
                                date,
                                valid_func)


@app.callback(
    Output('peptides', 'data', allow_duplicate=True),
    Output('experiment_name_field', 'value'),
    Output('peptides_import_format_alert', 'children'),
    Output('peptides_import_format_alert', 'is_open'),
    Input('peptides_obj_upload', 'contents'),
    Input('peptides_obj_upload', 'filename'),
    prevent_initial_call=True
)
def process_peptides_data(peptide_data, file_name):
    if peptide_data is None:
        raise PreventUpdate
    archive_format = determine_archive_format(file_name)
    data_str = memory_to_str(peptide_data, archive_format)

    valid, err_msg = MetaPepTable.validate_json(data_str)
    if not valid:
        return None, None, err_msg, True

    # experiment name taken from filename except last suffix
    file_name = "".join(file_name.split('.')[:-1])

    return data_str, file_name, None, False


@app.callback(
    Output('sample_name_import', 'disabled'),
    Input('merge_psm_switch', 'value'),
)
def inactivate_sample_name(merge_psm_switch):
    if merge_psm_switch is True:
        return False
    else:
        return True


@app.callback(
    Output('merge_psm_switch', 'disabled'),
    Output('merge_psm_switch', 'value'),
    Input('db_search_psm_valid', 'data'),
    Input('merge_psm_switch', 'value'),
)
def inactivate_psm_merge(db_search_valid, current_val):
    if db_search_valid is True:
        return (False, current_val)
    else:
        return (True, True)


@app.callback(
    Output('start_annotation_button', 'disabled'),
    Output('annotation_hint', 'children'),
    Input('db_search_psm_valid', 'data'),
    Input('de_novo_valid', 'data'),
    Input('taxonomy_db_valid', 'data'),
    Input('func_annot_db_valid', 'data'),
    Input('database_present_status', 'data'),
    Input('global_taxonomy_annotation_checkbox', 'value'),
    Input('sample_name_import', 'value'),
    State('merge_psm_switch', 'value')
)
def inactivate_annotation_button(psm_valid,
                                 denovo_valid,
                                 prot_db_valid,
                                 func_annot_db_valid,
                                 db_presence,
                                 global_tax_annot,
                                 sample_name,
                                 merge_psm):
    # ensure a sample name is given
    if sample_name is None and merge_psm is True:
        tooltip = dbc.Tooltip(f"Input sample name",
                                target="start_annotation_button_wrapper",
                                placement="bottom",
                                className="mt-1")
        return (True, [tooltip])


    # ensure metaproteomics data is provided
    if gc.display_db_search is True and gc.display_de_novo is True:
        tooltip_txt = "Import db search or de novo data"
    elif gc.display_db_search is True:
        tooltip_txt = "Import db search data"
    elif gc.display_de_novo is True:
        tooltip_txt = "Import de novo data"
    else:
        raise ValueError("invalid dashboard function mode, none of db search or de novo provided...")

    if all(i is not True for i in [psm_valid, denovo_valid]):
        tooltip = dbc.Tooltip(tooltip_txt,
                              target="start_annotation_button_wrapper",
                              placement="bottom",
                              className="mt-1")
        return (True, [tooltip])

    # if no db search data, global taxonomy annotation should be performed
    elif psm_valid is not True and global_tax_annot is False:
        tooltip = dbc.Tooltip(f"Add db search data, or check 'global annotation of peptides'",
                                target="start_annotation_button_wrapper",
                                placement="bottom",
                                className="mt-1")
        return (True, [tooltip])
    # finally, if db search data is present, at least one classification step has to be performed
    elif prot_db_valid is not True and\
        func_annot_db_valid is not True and\
        global_tax_annot is False:
        tooltip = dbc.Tooltip(f"Add taxonomy or function classification data, or check 'global annotation of peptides'",
                                target="start_annotation_button_wrapper",
                                placement="bottom",
                                className="mt-1")
        return (True, [tooltip])

    # valid ncbi taxonomy database is required for mapping of taxa to db search
    elif prot_db_valid is True and db_presence["ncbi_taxonomy"] is False:
        tooltip = dbc.Tooltip(f"Import ncbi taxonomy database or perform global taxonomy classification only",
            target="start_annotation_button_wrapper",
            placement="bottom",
            className="mt-1")
        return (True, [tooltip])

    else:
        tooltip = dbc.Tooltip(f"start annotation",
                                target="start_annotation_button_wrapper",
                                placement="bottom",
                                className="mt-1")
    return (False, [tooltip])


@app.callback(
    Output('tax_acc_map_delimiter_container', 'className'),
    Output('tax_acc_map_acc_type_container', 'className'),
    Output('tax_acc_map_acc_idx_container', 'className'),
    Output('tax_acc_map_tax_idx_container', 'className'),
    Output('tax_acc_map_tax_type_container', 'className'),
    Output('gtdb_genome_to_ncbi_container', 'className'),
    Output('global_taxonomy_annotation_checkbox', 'disabled'),
    Input('taxonomy_db_format', 'value'),
    Input('gtdb_genome_to_ncbi_checkbox', 'value')
)
def disable_taxonomy_annotations_options(tax_db_format, gtdb_to_ncbi):
    """Manage display and interactivity of taxonomy map settings depending on
    input format and dashboard version. 
    Accession pattern element is managed separately.
    """
    # Some elements will be disabled in "simple version"
    advanced = gc.show_advanced_settings

    # set style settings if items would be visible
    delim_style = "d-flex justify-content-start mb-3"
    acc_type_style = "d-flex justify-content-start align-items-center mb-3"
    acc_idx_style = "d-flex justify-content-start mb-3"
    tax_idx_style = "d-flex justify-content-start mb-3"
    tax_type_style = "d-flex justify-content-start align-items-center"
    gtbd_to_ncbi_style = "d-flex mt-4 justify-content-start align-items-center"
    unipept_checkbox = False

    # disable elements depending on formats and dashboard functionality
    # if GTDB format, only disable Unipept annotation switch if no conversion to NCBI
    if tax_db_format == "GTDB" and gtdb_to_ncbi is False:
        unipept_checkbox = True
    # For NCBI format, disable gtdb-ncbi conversion switch
    elif tax_db_format == "NCBI":
        gtbd_to_ncbi_style = "d-none"
    # For GhostKOALA format: hide all options except accession pattern
    elif tax_db_format == "GhostKOALA":
        delim_style = "d-none"
        acc_type_style = "d-none"
        acc_idx_style = "d-none"
        tax_idx_style = "d-none"
        tax_type_style = "d-none"
        gtbd_to_ncbi_style = "d-none"
    elif tax_db_format != "GTDB":
        raise ValueError("Invalid db format given...")

    return (
        delim_style,
        acc_type_style,
        acc_idx_style,
        tax_idx_style,
        tax_type_style,
        gtbd_to_ncbi_style,
        unipept_checkbox
    )

@app.callback(
    Output('acc_tax_map_acc_pattern_container', 'className'),
    Output('acc_tax_map_acc_pattern', 'value'),
    Input('accession_parser_checkbox', 'value'),
    State('acc_tax_map_acc_pattern', 'value'),
)
def taxonomy_annotations_accession_parser(parser_option, current_value):
    if parser_option == "Custom regex":
        return ("d-flex justify-content-start mb-3", current_value)
    elif parser_option == "Full string":
        return ("d-none", None)
    elif parser_option == "Up to first white-space":
        return ("d-none", r"\S+")


@app.callback(
    Output('func_annot_acc_pattern_container', 'className'),
    Output('func_annot_acc_pattern', 'value'),
    Input('func_accession_parser_checkbox', 'value'),
    State('func_annot_acc_pattern', 'value'),
)
def functional_annotations_accession_parser(parser_option, current_value):
    if parser_option == "Custom regex":
        return ("d-flex justify-content-start mb-3", current_value)
    elif parser_option == "Full string":
        return ("d-none", None)
    elif parser_option == "Up to first white-space":
        return ("d-none", r"\S+")


@app.callback(
    Output('db_search_accession_pattern_container', 'className'),
    Output('db_search_accession_pattern', 'value'),
    Input('db_search_accession_parser_items', 'value'),
    State('db_search_accession_pattern', 'value'),
)
def db_search_accession_parser(parser_option, current_value):
    if parser_option == "Custom regex":
        return ("d-flex justify-content-start mb-3", current_value)
    elif parser_option == "Full string":
        return ("d-none", None)
    elif parser_option == "Up to first white-space":
        return ("d-none", r"\S+")


#TODO: Wrap filter settings into metadata, store these in a dictionary dcc.Store object
@app.callback(
    Output('peptides', 'data', allow_duplicate=True),
    Output('peptides_metadata', 'data', allow_duplicate=True),
    Output('annotation_loading_spot', 'children'),
    Output('data_import_container', 'children'),
    Output('annotation_error_alert', 'children'),
    Output('annotation_error_alert', 'is_open'),
    Input('start_annotation_button', 'n_clicks'),
    State('peptides', 'data'),
    State('peptides_metadata', 'data'),
    State('sample_name_import', 'value'),
    State('merge_psm_switch', 'value'),
    State('db_search_psm_upload', 'contents'),
    State('db_search_psm_upload', 'filename'),
    State('db_search_psm_format', 'value'),
    State('db_search_accession_pattern', 'value'),
    State('db_search_psm_score_threshold', 'value'),
    State('db_search_filter_crap', 'value'),
    State('denovo_upload', 'contents'),
    State('denovo_upload', 'filename'),
    State('de_novo_format', 'value'),
    State('de_novo_score_threshold', 'value'),
    State('de_novo_filter_crap', 'value'),
    State('taxonomy_db_upload', 'contents'),
    State('taxonomy_db_format', 'value'),
    State('taxonomy_id_format_checkbox', 'value'),
    State('taxonomy_db_upload', 'filename'),
    State('acc_tax_map_delim', 'value'),
    State('tax_accession_format_radio', 'value'),
    State('acc_tax_map_acc_idx', 'value'),
    State('acc_tax_map_acc_pattern', 'value'),
    State('acc_tax_tax_idx', 'value'),
    State('gtdb_genome_to_ncbi_checkbox', 'value'),
    State('global_taxonomy_annotation_checkbox', 'value'),
    State('func_annot_db_upload', 'contents'),
    State('func_annot_db_format', 'value'),
    State('func_annot_db_upload', 'filename'),
    State('func_annot_combine', 'value'),
    State('func_annot_acc_pattern', 'value'),
    State('current_taxonomy_db_loc', 'data'),
    State('ncbi_taxonomy_db_loc', 'value'),
    prevent_initial_call=True
)
def process_manual_annotation(n_clicks,
                              current_peptides,
                              current_metadata,
                              sample_name,
                              merge_psms,
                              psm_list,
                              psm_names,
                              psm_format,
                              psm_acc_pattern,
                              psm_score_threshold,
                              psm_filter_crap,
                              denovo_data,
                              denovo_names,
                              denovo_format,
                              denovo_score_threshold,
                              denovo_filter_crap,
                              acc_tax_map,
                              acc_tax_map_format,
                              acc_tax_map_elem_format,
                              acc_tax_map_name,
                              acc_tax_map_delimiter,
                              acc_tax_map_acc_format,
                              acc_tax_map_acc_idx,
                              acc_tax_map_acc_pattern,
                              acc_tax_map_tax_idx,
                              gtdb_to_ncbi,
                              global_tax_annot,
                              func_annot_db,
                              func_annot_db_format,
                              func_annot_db_name,
                              func_annot_combine,
                              func_annot_acc_pattern,
                              tax_db_loc,
                              ncbi_tax_db_loc):
    """Annotate and combine imported datasets into one peptide dataset
    """
    # default state of output
    alert_content = None
    alert_open = False
    loading_status = None
    
    no_update = (current_peptides, 
                 current_metadata, 
                 loading_status, 
                 data_import_container, 
                 alert_content, 
                 alert_open)

    # keep multi file import consistent with empty list as opposed to nonetype
    psm_list = [] if psm_list is None else psm_list
    psm_names = [] if psm_names is None else psm_names
    denovo_data = [] if denovo_data is None else denovo_data
    denovo_names = [] if denovo_names is None else denovo_names

    # at least one db search or de novo dataset has to be imported to start
    if all(i is None or len(i) == 0 for i in [psm_list, denovo_data]):
        return no_update
    # if no db search data, global taxonomy annotation should be performed
    elif (psm_list is None or len(psm_list) == 0) and global_tax_annot is False:
        return no_update
    # finally, if db search data is present, at least one classification step has to be performed
    elif acc_tax_map is None and\
        func_annot_db is None and\
        global_tax_annot is False:
        return no_update

    # if archive given for psm files and de novo files, extract:
    if len(psm_list) != 0:
        try:
            archive_format = determine_archive_format(psm_names[0])
            if len(psm_list) == 1 and archive_format is not None:
                psm_list, psm_names = archive_to_file_list(psm_list[0],
                                                        archive_format)
            elif archive_format is not None:
                return no_update
        except:
            alert_msg = "Failed to extract DB search data. Ensure the archive only contains valid DB search files in consistent format."
            alert_open = True
            return (current_peptides, 
                    current_metadata, 
                    loading_status, 
                    data_import_container, 
                    alert_msg, 
                    alert_open)

    # extract denovo if archive given
    if len(denovo_data) != 0:
        try:
            archive_format = determine_archive_format(denovo_names[0])
            if len(denovo_data) == 1 and archive_format is not None:
                denovo_data, denovo_names = archive_to_file_list(denovo_data[0],
                                                                archive_format)
            elif archive_format is not None:
                return no_update
        except:
            alert_msg = "Failed to extract de novo data. Ensure the archive only contains valid de novo files in consistent format."
            alert_open = True
            return (current_peptides, 
                    current_metadata, 
                    loading_status, 
                    data_import_container, 
                    alert_msg, 
                    alert_open)

    # Import current peptides dataset, if data present, add new data
    if current_peptides is None:
        current_peptides_df = None
        current_sample_names = []
    else:
        try:
            current_peptides_df = MetaPepTable.read_json(current_peptides)
            current_sample_names = current_peptides_df.data["Sample Name"].unique().tolist()
        except:
            alert_msg = "Failed to parse project data..."
            alert_open = True
            return (current_peptides, 
                    current_metadata, 
                    loading_status, 
                    data_import_container, 
                    alert_msg, 
                    alert_open)
           
            
    # deduplicate (or filter) input file (names) and deduplicate sample name
    try:
        if merge_psms is False:
            sample_name, psm_list, psm_names = deduplicate_input_lists(current_sample_names,
                                                                       sample_name,
                                                                       psm_list,
                                                                       psm_names)
        else:
            sample_name = deduplicate_strings(sample_name, current_sample_names)
    except:
        alert_msg = "Failed to parse db search filenames..."
        alert_open = True
        return (current_peptides, 
                current_metadata, 
                loading_status, 
                data_import_container, 
                alert_msg, 
                alert_open)

    # set data formats to None if no data supplied
    if len(psm_list) == 0: psm_format = None
    if len(denovo_data) == 0: denovo_format = None
    if acc_tax_map is None and global_tax_annot is False: acc_tax_map_format = None
    if func_annot_db is None: func_annot_db_format = None

    # set variable for taxonomy database used (NCBI, GTDB)
    tax_db_format = "GTDB" if acc_tax_map_format == "GTDB" else "NCBI"

    # construct options object containing all filter settings
    options = AnnotationOptions(
        psm_format,
        psm_acc_pattern,
        psm_score_threshold,
        psm_filter_crap,
        denovo_format,
        denovo_score_threshold,
        denovo_filter_crap,
        GlobalConstants.min_pept_len,
        acc_tax_map_delimiter,
        acc_tax_map_name,
        acc_tax_map_acc_format,
        acc_tax_map_format,
        tax_db_format,
        acc_tax_map_elem_format,
        gtdb_to_ncbi,
        ncbi_tax_db_loc,
        global_tax_annot,                  # de novo tax search
        func_annot_db_name,
        func_annot_db_format,
        acc_tax_map_acc_pattern,
        acc_tax_map_acc_idx,
        acc_tax_map_tax_idx,
        func_annot_combine,
        func_annot_acc_pattern,
        merge_psms,
    )

    # perform taxonomy and functional annotation to psm data
    try:
        new_peptides = annotate_peptides(sample_name,
                                         psm_list,
                                         psm_names,
                                         denovo_data, # type: ignore
                                         acc_tax_map,
                                         func_annot_db,
                                         tax_db_loc,
                                         options)
    except AnnotationError as err:
        alert_open = True
        return (current_peptides, 
                current_metadata, 
                loading_status, 
                data_import_container, 
                err.msg, 
                alert_open)
    except Exception as err:
        alert_msg = "Failed annotation due to unexpected problem: {err}"
        alert_open = True
        return (current_peptides, 
                current_metadata, 
                loading_status, 
                data_import_container, 
                alert_msg, 
                alert_open)

    # merge new data to current data if present and store as json
    if current_peptides_df is not None:
        try:
            new_peptides = MetaPepTable.concat_tables([current_peptides_df, new_peptides])
        except ValueError as e:
            alert_msg = f"Failed to add sample to project table: {e}"
            alert_open = True
            return (current_peptides, 
                    current_metadata, 
                    loading_status, 
                    data_import_container, 
                    alert_msg, 
                    alert_open)

    try:
        new_peptides_dump = new_peptides.to_json()
    except:
        alert_msg = "failed to store project table..."
        alert_open = True
        return (current_peptides, 
                current_metadata, 
                loading_status, 
                data_import_container, 
                alert_msg, 
                alert_open)
        
        
    return (new_peptides_dump,
            current_metadata,
            None,
            data_import_container,
            None,
            False)


def global_taxonomy_annotation_only():
    ...



@app.callback(
    Output("download_peptides_csv", "data"),
    Input("export_peptides_csv", "n_clicks"),
    State("peptides", "data"),
    State("experiment_name_field", "value"),
    prevent_initial_call=True
)
def download_annotated_dataset_csv(n_clicks,
                                   peptide_json,
                                   experiment_name):
    if peptide_json is None:
        return

    if experiment_name is None:
        experiment_name = "peptides_df"

    # import json into metapep object
    peptides_obj = MetaPepTable.read_json(peptide_json)

    # Download dataframe from object
    return dcc.send_data_frame(peptides_obj.data.to_csv,
                               f"{experiment_name}.csv")


@app.callback(
    Output("download_peptides_json", "data"),
    Input("export_peptides_json", "n_clicks"),
    State("peptides", "data"),
    State("experiment_name_field", "value"),
    prevent_initial_call=True
)
def download_annotated_dataset_json(n_clicks,
                                    peptide_json,
                                    experiment_name):
    if peptide_json is None:
        return

    if experiment_name is None:
        experiment_name = "peptides"

    # import json into metapep object
    peptides_obj = MetaPepTable.read_json(peptide_json)
    file_name = f'{experiment_name}.json'

    # Download dataframe from object
    return dict(content=peptides_obj.to_json(), filename=file_name)


@app.callback(
    Output("current_taxonomy_db_loc", "data"),
    Input("ncbi_taxonomy_db_loc", "value"),
    Input("gtdb_taxonomy_db_loc", "value"),
    Input('taxonomy_db_format', 'value'),
)
def set_taxonomy_db_loc(ncbi_loc,
                        gtdb_loc,
                        taxonomy_format):
    match taxonomy_format:
        case "GTDB":
            return gtdb_loc
        case _:
            return ncbi_loc
