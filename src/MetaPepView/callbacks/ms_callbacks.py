from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import json
import numpy as np
from io import StringIO

from server import app

# import layout elements
from layout.style_constants import *
from layout.validation_page import *
from layout.func_annot_page import *

from backend import *
from backend.type_operations import *
from backend.plots import tic_over_rt_plot, \
    ms2_from_signal_arrays, \
    mz_over_rt_plot, \
    scan_tic_dist_plot, \
    ms1_int_over_ms2_int, \
    confidence_dist_plot, \
    ref_score_dist_plot, \
    ref_score_threshold_plot, \
    ref_intensity_dist_plot, \
    ref_score_threshold_barplot, \
    ref_transmission_scatter_plot


def show_spectra_name(name):
    """Display filename of annotated peptide dataset import.
    """
    # provide own implementation that does not look for nonetype at content
    # but at name, as content will not be stored in upload object
    if name is None:
        return html.P("No file...")
    # update name in sidebar
    else:
        if len(name) > 30:
            name = name[:30-3] + '...'
        return html.P(name, className="ms-1") 


@app.callback(
    Output('db_search_psm_qa_valid', 'data'),
    Output('db_search_psm_qa_name', 'children'),
    Output('db_search_psm_qa_upload', 'contents'),
    # Output('db_search_psm_qa_format_alert', 'children'),
    # Output('db_search_psm_qa_format_alert', 'is_open'),
    Output('db_search_psm_qa_import_box', 'style'),
    Input('db_search_psm_qa_upload', 'contents'),
    Input('db_search_psm_qa_format', 'value'),
    State('db_search_psm_qa_upload', 'filename'),
    State('db_search_psm_qa_upload', 'last_modified'))
def show_db_psm_search_qa_name(contents, file_format, name, date):
    """Display filename of annotated peptide dataset import.
    """
    print("parse db search")
    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)
        return validate_db_search(cont_buf, file_format)
    
    valid_data, name, content, err_msg, success, import_box_style = validate_single_file(contents, name, date, valid_func, drag_and_drop=True)
    
    # update import box style
    qa_box_style = qa_import_box_style
    if "background-color" in import_box_style.keys():
        qa_box_style["background-color"] = import_box_style["background-color"]
    
    return (valid_data, name, content, qa_box_style)


@app.callback(
    Output('denovo_qa_valid', 'data'),
    Output('denovo_qa_name', 'children'),
    Output('denovo_qa_upload', 'contents'),
    # Output('denovo_qa_format_alert', 'children'),
    # Output('denovo_qa_format_alert', 'is_open'),
    Output('denovo_qa_import_box', 'style'),
    Input('denovo_qa_upload', 'contents'),
    Input('denovo_qa_format', 'value'),
    State('denovo_qa_upload', 'filename'),
    State('denovo_qa_upload', 'last_modified'))
def show_denovo_search_qa_name(contents, file_format, name, date):
    """Display filename of annotated peptide dataset import.
    """
    print("parse de novo")
    # set validation function
    def valid_func(cont, archv) -> Tuple[bool, str | None]:
        cont_buf = memory_to_stringio(cont, archv)
        return validate_de_novo(cont_buf, file_format)
    valid_data, name, content, err_msg, success, import_box_style = validate_single_file(contents, name, date, valid_func, drag_and_drop=True)

    # update import box style
    qa_box_style = qa_import_box_style
    if "background-color" in import_box_style.keys():
        qa_box_style["background-color"] = import_box_style["background-color"]
    
    print(f"de novo valid? {valid_data}")
    return (valid_data, name, content, qa_box_style)


@app.callback(
    Output("mzml_data", "data"),
    Output("mzml_peaks_data", "data"),
    Output("mzml_metadata", "data"),
    Output("mzml_upload", "contents"),
    Output("mzml_name", "children"),
    Output('mzml_valid', 'data'),
    Output('mzml_import_box', 'style'),
    Input("mzml_upload", "contents"),
    Input("mzml_upload", "filename"),
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
    
    qa_box_style = qa_import_box_style
    
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
        
        qa_box_style["background-color"] = IMPORT_FAILED_COLOR
        return (
            None, 
            None, 
            None, 
            None, 
            import_single_file(None, None, drag_and_drop=True),
            False,
            qa_box_style
        )
    print("Done wrangling...")
    # After dataset is imported, remove upload data to save memory
    # However, do keep the filename for display in dashboard
    
    qa_box_style["background-color"] = IMPORT_SUCCESS_COLOR
    return (
        compress_string(data.to_json()), 
        compress_string(peaks_data), 
        metadata, 
        None, 
        import_single_file(filename, None, drag_and_drop=True),
        True,
        qa_box_style
    )

@app.callback(
    Output("tic_sec_param_int_cutoff_container", "hidden"),
    Output("tic_sec_param_conf_cutoff_container", "hidden"),
    Input("tic_ms_secondary_y", "value"),
)
def hide_tic_over_rt_options(secondary_param):
    if secondary_param == "Peak Count":
        return (False, True)
    elif secondary_param == "DB Search Confidence":
        return (True, False)
    elif secondary_param == "De Novo Confidence":
        return (True, False)
    else:
        return (True, True)

@app.callback(
    Output("tic_over_rt_div", "children"),
    Output("tic_over_rt_div", "style"),
    Input("mzml_data", "data"),
    Input("mzml_peaks_data", "data"),
    Input("mzml_metadata", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("tic_ms_level", "value"),
    Input("tic_ms_sma_range", "value"),
    Input("tic_ms_red_fact", "value"),
    Input("tic_ms_secondary_y", "value"),
    Input("tic_sec_param_int_cutoff", "value"),
    Input("tic_sec_param_conf_cutoff", "value"),
)
def show_tic_over_rt(dataset,
                     peaks,
                     mzml_metadata,
                     db_search_psm,
                     db_search_psm_format,
                     de_novo,
                     de_novo_format,
                     ms_level,
                     sma_range,
                     reduction_factor,
                     secondary_param,
                     peak_int_cutoff,
                     metapep_confidence_cutoff):
    if peak_int_cutoff is None:
        peak_int_cutoff = 0
    if metapep_confidence_cutoff is None:
        metapep_confidence_cutoff = 0
    
    # only update once mzml is uploaded
    if dataset is None:
        raise PreventUpdate
    dataset = decompress_string(dataset)
    dataset = pd.read_json(StringIO(dataset))
    
    # only load peaks data if required
    prot_data = None
    if secondary_param == "Peak Count" and peak_int_cutoff > 0:
        peaks = decompress_string(peaks)
    # only load metapep data if required
    elif secondary_param == "DB Search Confidence" and db_search_psm is not None:
        prot_data = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
        secondary_param = "Confidence"
    elif secondary_param == "De Novo Confidence" and de_novo is not None:
        prot_data = load_metapep_de_novo(de_novo, "filename", de_novo_format)
        secondary_param = "Confidence"

    # Obtain TIC + RT for all MS1 spectra
    fig = tic_over_rt_plot(dataset,
                           peaks,
                           mzml_metadata['compression type'],
                           mzml_metadata['binary type'],
                           int(ms_level),
                           secondary_param,
                           sma_range,
                           reduction_factor,
                           prot_data,
                           metapep_confidence_cutoff,
                           peak_int_cutoff)
    
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="tic_over_rt_fig", style={'height': '100%'})
    
    print("finish")
    return (graph, {"display": 'block', 'height': '20rem'})


@app.callback(
    Output("ms2_spec_div", "children"),
    Output("ms2_spec_div", "style"),
    Input("tic_over_rt_fig", "clickData"),
    State("mzml_data", "data"),
    State("mzml_peaks_data", "data")
)
def show_ms2_spectrum(tic_data, content, peaks):
    print("tic_data:", tic_data)
    if content is None or tic_data is None:
        raise PreventUpdate
    content = decompress_string(content)
    peaks = decompress_string(peaks)
    
    print(tic_data)
    
    print("start ms2 processing")
    content = pd.read_json(StringIO(content))
    click_rt = tic_data["points"][0]["x"]
    
    xvals = []
    yvals = []
    
    print("parse content")
    # parse mzxml for closest ms2 scan
    for idx, scan in content.iterrows():
            
        # skip scan if mslevel not of interest
        if scan['MS level'] != 2:
            continue
        
        rt = scan['retention time']
        
        # if scan closest to selected rt, retrieve spectra data
        if float(rt) > click_rt:
            peaks = json.loads(peaks).get(str(idx))
            if peaks is not None:
                mz_array = peaks["m/z array"]
                int_array = peaks["intensity array"]
                peak_number = scan['peaks count']
                xvals = decode_mzml_peaks(mz_array, peak_number)
                yvals = decode_mzml_peaks(int_array, peak_number)
            break
    
    print("show peaks")
    fig = ms2_from_signal_arrays(xvals, yvals)
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="ms2_spec_fig", style={'height': '100%'})
    
    return (graph, {"display": 'block', "height": "19rem"})


@app.callback(
    Output("mz_over_rt_div", "children"),
    Output("mz_over_rt_div", "style"),
    Input("mzml_data", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("mz_over_rt_int_cutoff", "value"),
    Input("mz_over_rt_ident_frac", "value")
)
def show_mz_over_rt(mzml_content,
                    db_search_psm,
                    db_search_psm_format,
                    de_novo,
                    de_novo_format,
                    int_cutoff,
                    ident_val):
    # only update once mzxml is uploaded
    if mzml_content is None:
        raise PreventUpdate
    
    mzml_content = decompress_string(mzml_content)
    
    mzml_content = pd.read_json(StringIO(mzml_content))
    
    # add identification data based on selected option and delivered datasets
    if ident_val == "DB search" and db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    if ident_val == "De novo" and de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
    
    fig = mz_over_rt_plot(mzml_content, db_search_psm, de_novo, int_cutoff)
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="tic_over_rt_fig", style={'height': '100%'})
    return (graph, {"display": 'block', 'height': '19rem'})


@app.callback(
    Output("mz_over_rt_ident_frac", "disabled"),
    Output("mz_over_rt_ident_frac", "options"),
    Input("db_search_psm_qa_valid", "data"),
    Input("denovo_qa_valid", "data"),
)
def update_ident_dropdown_mz_over_rt(db_search_psm_valid,
                                     de_novo_valid):
    disable = True
    options = ["None"]
    
    print(f"db search psm valid in mz over rt: {db_search_psm_valid}")
    print(f"de novo valid in mz over rt: {de_novo_valid}")
    if db_search_psm_valid is True:
        options.append("DB search")
        disable = False
    if de_novo_valid is True:
        options.append("De novo")
        disable = False
        
    return disable, options


@app.callback(
    Output("scan_int_dist_alc_cutoff", "disabled"),
    Input("denovo_qa_valid", "data"),
)
def update_alc_cutoff_scan_int(de_novo_valid):
    print(f"de novo valid in alc cutoff: {de_novo_valid}")
    if de_novo_valid is True:
        return False
    else:
        return True


@app.callback(
    Output("ms1_over_ms2_ident_frac", "disabled"),
    Output("ms1_over_ms2_ident_frac", "options"),
    Input("db_search_psm_qa_valid", "data"),
    Input("denovo_qa_valid", "data"),
)
def update_ident_dropdown_ms1_ms2(db_search_psm_valid,
                                  de_novo_valid):
    disable = True
    options = ["None"]

    if db_search_psm_valid is True:
        options.append("DB search")
        disable = False
    if de_novo_valid is True:
        options.append("De novo")
        disable = False
        
    return disable, options


@app.callback(
    Output("scan_int_dist_div", "children"),
    Output("scan_int_dist_div", "style"),
    Input("mzml_data", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("scan_int_dist_ms_level", "value"),
    Input("scan_int_dist_alc_cutoff", "value")
)
def show_tic_dist(mzml_content,
                  db_search_psm,
                  db_search_psm_format,
                  de_novo,
                  de_novo_format,
                  ms_level,
                  alc_cutoff):
    # only update once mzxml is uploaded
    if mzml_content is None:
        raise PreventUpdate
    
    mzml_content = decompress_string(mzml_content)
    
    mzml_content = pd.read_json(StringIO(mzml_content))
    
    # add identification data based on selected option and delivered datasets
    if db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    if de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
    
    
    fig = scan_tic_dist_plot(mzml_content, int(ms_level), db_search_psm, de_novo, alc_cutoff)
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="scan_int_dist_fig", style={'height': '100%'})
    return (graph, {"display": 'block', 'height': '19rem'})


@app.callback(
    Output("ms1_over_ms2_int_div", "children"),
    Output("ms1_over_ms2_int_div", "style"),
    Input("mzml_data", "data"),
    Input("mzml_peaks_data", "data"),
    Input("mzml_metadata", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("ms1_over_ms2_mz_cutoff", "value"),
    Input("ms1_over_ms2_ident_frac", "value")
)
def show_frag_eff(mzml_content,
                  peaks,
                  metadata,
                  db_search_psm,
                  db_search_psm_format,
                  de_novo,
                  de_novo_format,
                  mz_cutoff,
                  ident_frac):
    # only update once mzxml is uploaded
    if mzml_content is None:
        raise PreventUpdate
    
    mzml_content = decompress_string(mzml_content)
    peaks = decompress_string(peaks)
    
    mzml_content = pd.read_json(StringIO(mzml_content))
   
    # add identification data based on selected option and delivered datasets
    if ident_frac == "DB search" and db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    if ident_frac == "De novo" and de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
   
    
    fig = ms1_int_over_ms2_int(mzml_content,
                               peaks,
                               db_search_psm,
                               de_novo,
                               min_mz=mz_cutoff,
                               peaks_compression=metadata['compression type'],
                               peaks_precision=metadata['binary type'])
    
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="ms1_over_ms2_int_fig", style={'height': '100%'})
    return (graph, {"display": 'block', 'height': '19rem'})


@app.callback(
    Output("pept_confidence_div", "children"),
    Output("pept_confidence_div", "style"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("pept_conf_metric", "value"),
    Input("pept_conf_lgp_cutoff", "value"),
    Input("pept_conf_alc_cutoff", "value")
)
def show_confidence_dist(db_search_psm,
                         db_search_psm_format,
                         de_novo,
                         de_novo_format,
                         metric,
                         lgp_cutoff,
                         alc_cutoff):
    # wait for valid value of cutoff items
    if lgp_cutoff is None or alc_cutoff is None:
        raise PreventUpdate
    
    if metric in ["DB search", "DB search / De novo", "All"] and db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    
    if metric in ["De novo", "DB search / De novo", "All"] and de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
    
    if metric == "All":
        plot_all = True
    else:
        plot_all = False
        
    fig = confidence_dist_plot(db_search_psm, de_novo, lgp_cutoff, alc_cutoff, plot_all)
    if fig is None or lgp_cutoff is None or alc_cutoff is None:
        raise PreventUpdate
    
    fig.update_layout(autosize=True)
    graph = dcc.Graph(figure=fig, id="pept_confidence_fig", style={'height': '100%'})
    return (graph, {"display": 'block', 'height': '24rem'})
   
   

################################################################################ 
# Reference Benchmark callbacks
################################################################################ 

@app.callback(
    Output("ref_statistics_dropdown", "options"),
    Input("ref_statistics", "data")
)
def add_ref_files_dropdown(ref_data):
    ref_dict = json.loads(ref_data)
    
    samples = list(ref_dict.keys())
    return samples


@app.callback(
    Output("current_ref_statistics_store", "data"),
    Output("ref_statistics_db_search_format", "children"),
    Output("ref_statistics_de_novo_format", "children"),
    Input("ref_statistics", "data"),
    Input("ref_statistics_dropdown", "value"),
    Input("custom_ref_dataset", "contents")
)
def load_ref_data(total_ref_stat, ref_dropdown_option, custom_ref):
    if ref_dropdown_option is None and custom_ref is None:
        return None, "...", "..."
    
    # directly extract sample based on option key
    if custom_ref is not None:
        ref_dict = json.loads(memory_to_str(custom_ref))
    else:
        ref_dict = json.loads(total_ref_stat)[ref_dropdown_option]
    
    try: 
        metadata = ref_dict["metadata"]
        db_search_format = metadata["db search format"]
        de_novo_format = metadata["de novo format"]
    except KeyError:
        raise ValueError("Metadata not present in experimental benchmark dataset")
    
        
    return json.dumps(ref_dict), db_search_format, de_novo_format


@app.callback(
    Output("custom_ref_statistics_name", "children"),
    Input("custom_ref_dataset", "filename"),
    State("custom_ref_dataset", "last_modified"),
)
def show_custom_ref_name(names, dates):
    return import_single_file(names, dates, max_name_len=40, drag_and_drop=False)


@app.callback(
    Output("reference_confidence_div", "children"),
    Output("reference_confidence_div", "style"),
    Input("current_ref_statistics_store", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("mzml_metadata", "data"),
    Input("ref_score_conf_metric", "value"),
    Input("ref_score_conf_x_scaling", "value"),
)
def show_reference_score_dist(ref_data,
                              db_search_psm,
                              db_search_psm_format,
                              de_novo,
                              de_novo_format,
                              spectral_metadata,
                              score_type,
                              x_normalization):
    if ref_data is None:
        raise PreventUpdate
    
    if db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    
    if de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
        
    if spectral_metadata is not None:
        sample_ms2_count = spectral_metadata.get('MS2 spectrum count')
    else:
        sample_ms2_count = None
        
    ref_dict = json.loads(ref_data)
    
    # set normalization params based on selected option
    match_norm, ms2_norm = False, False
    if x_normalization == "normalize_matches": match_norm = True
    elif x_normalization == "normalize_ms2": ms2_norm = True
    
    fig = ref_score_dist_plot(stat_dict=ref_dict,
                              sample_db_search=db_search_psm,
                              sample_de_novo=de_novo,
                              format=score_type,
                              normalize_matches=match_norm,
                              normalize_scans=ms2_norm,
                              sample_ms2_count=sample_ms2_count)
    
    graph = dcc.Graph(figure=fig,
                      id="reference_confidence_fig",
                      style={'height': '100%'})
    
    return (graph, {"display": 'block', 'height': '24rem'})


@app.callback(
    Output("reference_score_thres_div", "children"),
    Output("reference_score_thres_div", "style"),
    Input("current_ref_statistics_store", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("mzml_metadata", "data"),
    Input("ref_conf_dist_norm", "value")
)
def show_score_threshold_dist(ref_data,
                              db_search_psm,
                              db_search_psm_format,
                              de_novo,
                              de_novo_format,
                              spectral_metadata,
                              normalization_option):
    if ref_data is None:
        raise PreventUpdate
    
    # configure score formats to process
    formats = ['db search', 'de novo', 'de novo only']
    
    if db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    
    if de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
    
    # configure normalization
    norm_scan, norm_rt = False, False
    if normalization_option == "scans":
        norm_scan = True
    elif normalization_option == "RT":
        norm_rt = True
    
    # directly extract sample based on option key
    ref_dict = json.loads(ref_data)
    
    fig = ref_score_threshold_plot(stat_dict=ref_dict,
                                   formats=formats,
                                   normalize_psm=norm_scan,
                                   normalize_rt=norm_rt,
                                   sample_db_search=db_search_psm,
                                   sample_de_novo=de_novo,
                                   spectral_metadata=spectral_metadata)
    
    graph = dcc.Graph(figure=fig,
                      id="reference_score_thres_fig",
                      style={'height': '100%'})
    
    return (graph, {"display": 'block', 'height': '24rem'})


@app.callback(
    Output("reference_intensity_div", "children"),
    Output("reference_intensity_div", "style"),
    Input("current_ref_statistics_store", "data"),
    Input("mzml_data", "data"),
)
def show_intensity_dist(ref_data,
                        mzml_content):
    if ref_data is None:
        raise PreventUpdate
    
    # only update once mzml is uploaded
    if mzml_content is not None:
        mzml_content = decompress_string(mzml_content)
        mzml_content = pd.read_json(StringIO(mzml_content))

    # directly extract sample based on option key
    ref_dict = json.loads(ref_data)
    
    fig = ref_intensity_dist_plot(stat_dict=ref_dict, spectral_data=mzml_content)
    
    graph = dcc.Graph(figure=fig,
                      id="reference_intensity_fig",
                      style={'height': '100%'})
    
    return (graph, {"display": 'block', 'height': '24rem'})


@app.callback(
    Output("reference_transmission_div", "children"),
    Output("reference_transmission_div", "style"),
    Input("current_ref_statistics_store", "data"),
    Input("mzml_data", "data"),
    Input("scale_ion_injection_time", "value")
)
def show_transmission_dist(ref_data,
                           mzml_content,
                           scale_ion_inj):
    if ref_data is None:
        raise PreventUpdate
    
    # only update once mzml is uploaded
    if mzml_content is not None:
        mzml_content = decompress_string(mzml_content)
        mzml_content = pd.read_json(StringIO(mzml_content))
    
    # directly extract sample based on option key
    ref_dict = json.loads(ref_data)
    
    fig = ref_transmission_scatter_plot(ref_dict, mzml_content, scale_ion_inj)
    
    graph = dcc.Graph(figure=fig,
                      id="reference_transmission_fig",
                      style={'height': '100%'})
    
    return (graph, {"display": 'block', 'height': '24rem'})


@app.callback(
    Output("reference_score_barplot_div", "children"),
    Output("reference_score_barplot_div", "style"),
    Input("current_ref_statistics_store", "data"),
    Input("db_search_psm_qa_upload", "contents"),
    Input("db_search_psm_qa_format", "value"),
    Input("denovo_qa_upload", "contents"),
    Input("denovo_qa_format", "value"),
    Input("mzml_metadata", "data"),
    Input("threshold_barplot_scaling", "value"),
    Input("threshold_barplot_fill_bars", "value")
)
def show_score_barplot_dist(ref_data,
                            db_search_psm,
                            db_search_psm_format,
                            de_novo,
                            de_novo_format,
                            spectral_metadata,
                            normalization,
                            fill_bars):
    if ref_data is None:
        raise PreventUpdate
    
    # configure score formats to process
    formats = ['db search', 'de novo', 'de novo only']
    
    if db_search_psm is not None:
        db_search_psm = load_metapep_db_search(db_search_psm, "filename", db_search_psm_format)
    else:
        db_search_psm = None
    
    if de_novo is not None:
        de_novo = load_metapep_de_novo(de_novo, "filename", de_novo_format)
    else:
        de_novo = None
    
    # set normalization options   
    normalize_scans = False
    normalize_rt = False
    if normalization == "normalize_psm":
        normalize_scans = True
    elif normalization == "normalize_rt":
        normalize_rt = True
    
    # directly extract sample based on option key
    ref_dict = json.loads(ref_data)
    
    fig = ref_score_threshold_barplot(stat_dict=ref_dict,
                                      formats=formats,
                                      normalize_psm=normalize_scans,
                                      normalize_rt=normalize_rt,
                                      normalize_fill=fill_bars,
                                      sample_db_search=db_search_psm,
                                      sample_de_novo=de_novo,
                                      spectral_metadata=spectral_metadata)
    
    graph = dcc.Graph(figure=fig,
                      id="reference_score_barplot_fig",
                      style={'height': '100%'})
    
    return (graph, {"display": 'block', 'height': '24rem'})
