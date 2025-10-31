from dash import html, dcc
import dash_bootstrap_components as dbc

from metapepview.html_templates import qa_importer_block
from metapepview.constants import GlobalConstants, StyleConstants



mzml_import = qa_importer_block(
    "Spectral File (mzML)",
    "mzml_upload",
    "mzml_name",
    "mzml_valid",
)

features_import = qa_importer_block(
    "Features (featureXML)",
    "features_upload",
    "features_name",
    "features_valid",
)


psm_file_validation_import = qa_importer_block(
    "DB search Import",
    "db_search_psm_qa_upload",
    "db_search_psm_qa_name",
    "db_search_psm_qa_valid",
    format_options=GlobalConstants.db_search_dropdown_options,
    format_id="db_search_psm_qa_format"
)


denovo_file_validation_import = qa_importer_block(
    "De Novo Import",
    "denovo_qa_upload",
    "denovo_qa_name",
    "denovo_qa_valid",
    format_options=GlobalConstants.de_novo_dropdown_options,
    format_id="denovo_qa_format"
)


ref_set_import = qa_importer_block(
    "Reference Data Import",
    "ref_prot_upload",
    "ref_prot_name",
    "ref_prot_valid",
)



################################################################################
# Figure blocks
################################################################################


tic_over_rt = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Total ion current", style={"margin-right": "3rem"}),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("MS level:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Select(options=[{"label": "MS1", "value": 1}, {"label": "MS2", "value": 2}],
                                                    value=1, id="tic_ms_level", style={"width": "6rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("SMA window:", 
                                               id="sma_window_text",
                                               style={"text-align": "center", 
                                                      "margin-right": "1rem",
                                                      "margin-left": "1rem",
                                                      "text-decoration-line": "underline", 
                                                      "text-decoration-style": "dotted"}),
                                        dbc.Input(type="number",
                                                  min=1,
                                                  max=1000,
                                                  step=1,
                                                  value=20,
                                                  id="tic_ms_sma_range",
                                                  debounce=True,
                                                  style={'width': '6rem'}),
                                        dbc.Popover("""
                                            Set window size for single moving average smoothing
                                            """,
                                            id="sma_window_popover",
                                            target="sma_window_text",
                                            trigger="hover",
                                            placement='bottom',
                                            className="p-2"
                                        )       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Data reduction factor:", 
                                               id="data_red_fact_text",
                                               style={"text-align": "center", 
                                                      "margin-right": "1rem", 
                                                      "margin-left": "1rem",
                                                      "text-decoration-line": "underline", 
                                                      "text-decoration-style": "dotted"}),
                                        dbc.Input(type="number",
                                                  min=1,
                                                  max=10,
                                                  step=1,
                                                  value=3,
                                                  id="tic_ms_red_fact",
                                                  debounce=True,
                                                  style={'width': '4rem'}),
                                        dbc.Popover("""
                                            For computational speedup, keep only every n'th data point, discarding the rest.
                                            """,
                                            id="data_red_factor_popover",
                                            target="data_red_fact_text",
                                            trigger="hover",
                                            placement='bottom',
                                            className="p-2"
                                        )        
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Right axis param:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dcc.Dropdown(
                                            ["None",
                                             "DB Search Counts",
                                             "De Novo Counts",
                                             "Peak Count",
                                             "Peak Width (FWHM)",
                                             "Feature Quality",
                                             "Ion injection time",
                                             "topN MS2"],
                                            value="None",
                                            id="tic_ms_secondary_y",
                                            style={"width": "12rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("int. cutoff:", 
                                               id="tic_over_rt_int_text",
                                               style={"text-align": "center", 
                                                      "margin-right": "1rem", 
                                                      "margin-left": "1rem",
                                                      "text-decoration-line": "underline", 
                                                      "text-decoration-style": "dotted"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  step=1,
                                                  value=0,
                                                  id="tic_sec_param_int_cutoff",
                                                  debounce=True,
                                                  style={'width': '8rem'}),
                                        dbc.Popover("""
                                            Ignore peaks below intensity cutoff
                                            """,
                                            id="tic_over_rt_int_popover",
                                            target="tic_over_rt_int_text",
                                            trigger="hover",
                                            placement='bottom',
                                            className="p-2"
                                        )             
                                    ],
                                    id="tic_sec_param_int_cutoff_container",
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"},
                                    hidden=True
                                ),
                                html.Div(
                                    [
                                        html.B("Confidence cutoff:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  step=0.01,
                                                  value=0,
                                                  id="tic_sec_param_conf_cutoff",
                                                  debounce=True,
                                                  style={'width': '6rem'}),       
                                    ],
                                    id="tic_sec_param_conf_cutoff_container",
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"},
                                    hidden=True
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        )
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="tic_over_rt_fig"), id="tic_over_rt_div", style={'display': 'none'}),
            ], style={"display": "block", "flex-wrap": "wrap"}
        )
    ],
    color="light",
    style={"margin": "0.5rem 0.5rem", "height": "30rem"}
)


feature_map = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("m/z over RT"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("int. cutoff:", 
                                               id="mz_over_rt_int_cutoff_text",
                                               style={"margin-right": "1rem", 
                                                      "text-align": "center",
                                                      "text-decoration-line": "underline", 
                                                      "text-decoration-style": "dotted"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  max=1e8, 
                                                  step=1, 
                                                  value=0, 
                                                  id="mz_over_rt_int_cutoff", 
                                                  debounce=True,
                                                  style={'width': '8rem'}),
                                        dbc.Popover("""
                                            Ignore peaks below intensity cutoff
                                            """,
                                            id="mz_over_rt_int_popover",
                                            target="mz_over_rt_int_cutoff_text",
                                            trigger="hover",
                                            placement='bottom',
                                            className="p-2"
                                        )             
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Show identified:", style={"margin-right": "1rem", "margin-left": "1rem", "text-align": "center"}),
                                        dbc.Select(options=["None", "DB search", "De novo"],
                                                    value="None", disabled=True, id="mz_over_rt_ident_frac", style={"width": "8rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        )
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="mz_over_rt_fig"), id="mz_over_rt_div", style={'display': 'none'}),
            ]
        )
    ],
    color="light",
    style={"height": "25rem"}
)


intensity_hist = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Identification distribution"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("MS level:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Select(options=[{"label": "MS1", "value": 1}, {"label": "MS2", "value": 2}],
                                                    value=2, id="scan_int_dist_ms_level", style={"width": "6rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Confidence cutoff (de novo):", style={"margin-right": "1rem", "margin-left": "1rem", "text-align": "center"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  max=100,
                                                  step=0.01,
                                                  value=80,
                                                  disabled=True,
                                                  id="scan_int_dist_alc_cutoff",
                                                  debounce=True,
                                                  style={'width': '5rem'}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("normalize bars:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Checkbox(
                                            id="int_dist_fig_norm_bars",
                                            disabled=False,
                                            value=False,
                                            #style={"width": "15rem"}
                                        ),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        )
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="scan_int_dist_fig"), id="scan_int_dist_div", style={'display': 'none'}),
            ]
        )
    ],
    color="light",
    style={"height": "25rem"}
)


charge_dist = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Charge distribution (features)"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("Confidence cutoff (de novo):", style={"margin-right": "1rem", "margin-left": "1rem", "text-align": "center"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  max=100,
                                                  step=1,
                                                  value=80,
                                                  disabled=True,
                                                  id="feature_charge_dist_alc_cutoff",
                                                  debounce=True,
                                                  style={'width': '5rem'}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        )
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="feature_charge_dist_fig"), id="feature_charge_dist_div", style={'display': 'none'}),
            ]
        )
    ],
    color="light",
    style={"height": "25rem"}
)


transmission_scatter = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("MS2 ion transmission"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("min m/z:", 
                                               id="ion_transm_min_mz_text",
                                               style={"margin-right": "1rem", 
                                                      "text-align": "center",
                                                      "text-decoration-line": "underline", 
                                                      "text-decoration-style": "dotted"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  max=1400,
                                                  step=1,
                                                  value=0,
                                                  id="ms1_over_ms2_mz_cutoff",
                                                  debounce=True,
                                                  style={'width': '6rem'}),
                                        dbc.Popover("""
                                            Ignore signals below m/z threshold in MS2 intensity calculation
                                            """,
                                            id="ion transm_popover",
                                            target="ion_transm_min_mz_text",
                                            trigger="hover",
                                            placement='bottom',
                                            className="p-2"
                                        )    ,
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Show identified:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Select(options=["None", "DB search", "De novo"],
                                                    value="None", disabled=True, id="ms1_over_ms2_ident_frac", style={"width": "8rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        )
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="ms1_over_ms2_int_fig"), id="ms1_over_ms2_int_div", style={'display': 'none'}),
            ]
        )
    ],
    color="light",
    style={"height": "25rem"}
)


ident_score_rank_dist = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Confidence Distribution", style={"margin-right": "3rem"}),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("Metric:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Select(options=["DB search",
                                                            "De novo",
                                                            "DB search / De novo",
                                                            "All"],
                                                    value="DB search", id="pept_conf_metric", style={"width": "14rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Confidence cutoff DB search:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  step=1,
                                                  value=35,
                                                  id="pept_conf_lgp_cutoff",
                                                  debounce=True,
                                                  style={'width': '6rem'}),       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Confidence cutoff de novo:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number",
                                                  min=0,
                                                  max=200,
                                                  step=1,
                                                  value=35,
                                                  id="pept_conf_alc_cutoff",
                                                  debounce=True,
                                                  style={'width': '6rem'}),       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "justify-content": "space-evenly"}
                        ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="pept_confidence_fig"), id="pept_confidence_div", style={'display': 'none'}),
            ], style={"display": "block", "flex-wrap": "wrap"}
        )
    ],
    color="light",
    style={"margin": "0.5rem 0.5rem", "height": "35rem"}
)


ref_confidence_dist_line = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Confidence ranked distribution"),
                        html.Div(
                            [
                                html.B("Metric:", style={"margin-right": "1rem", "text-align": "center"}),
                                dbc.Select(options=[{"label": "DB search", "value": "db search"},
                                                    {"label": "De novo", "value": "de novo"}],
                                            value="db search",
                                            id="ref_score_conf_metric",
                                            style={"width": "14rem", 
                                                   "margin-right": "4rem"}),
                                html.B("x-axis normalization:", style={"margin-right": "1rem", "text-align": "center"}),
                                dbc.RadioItems(options=[{"label": "None", "value": None},
                                                        {"label": "Total matches", "value": "normalize_matches"},
                                                        {"label": "Total ms2", "value": "normalize_ms2"}],
                                                inline=True,
                                                value=None,
                                                id="ref_score_conf_x_scaling")
                            ],
                            style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                        ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_confidence_fig"), id="reference_confidence_div", style={'display': 'none'})
            ],
            style={"display": "block", "flex-wrap": "wrap"}
        )
    ],
    color="light",
    style={"margin": "1rem 0rem 0.5rem 1rem", "height": "30rem"}
)


ref_confidence_dist_scatter = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Peptide confidence scatterplot"),
                        html.Div(
                            [
                                html.B("Normalize:", style={"margin-right": "1rem", "text-align": "center"}),
                                dbc.Select(options=[{"label": "None", "value": None},
                                                    {"label": "MS2 scans", "value": "scans"},
                                                    {"label": "Retention time", "value": "RT"}],
                                            value="db search psm", id="ref_conf_dist_norm", style={"width": "14rem"}),
                            ],
                            style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                        ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_score_thres_fig"), id="reference_score_thres_div", style={'display': 'none'})
            ]
        )
    ],
    color="light",
    style={"margin": ".5rem 0rem 0rem 1rem", "height": "30rem"}
)


ref_intensity_dist_scatter = dbc.Card(
    [
        dbc.CardBody(
            [
                html.H4("Scan intensity distribution"),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_intensity_fig"), id="reference_intensity_div", style={'display': 'none'})
            ]
        )
    ],
    color="light",
    style={"margin": ".5rem 0rem 0rem 0rem", "height": "30rem"}
)


ref_transmission_scatter = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Ion transmission"),
                        html.Div(
                            [
                                html.B("Scale ion injection time:", style={"margin-right": "1rem", "text-align": "center"}),
                                dbc.Checkbox(
                                    id="scale_ion_injection_time",
                                    value=False,
                                    #style={"width": "15rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                        ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_transmission_fig"), id="reference_transmission_div", style={'display': 'none'})
            ]
        )
    ],
    color="light",
    style={"margin": "1rem 0rem 0rem 1rem", "height": "30rem"}
)


ref_miscleavage_dist = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Miscleavage distribution (DB search)"),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="miscleavage_dist_fig"), id="miscleavage_dist_div", style={'display': 'none'})
            ]
        )
    ],
    color="light",
    style={"margin": "1rem 0rem 0rem 0rem", "height": "30rem"}
)


ref_metrics_comp = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Experiment metrics overview"),
                        # html.Div(
                        #     [

                        #     ],
                        #     style={'display': 'flex', "align-items": "center"} 
                        # ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_metrics_scores_fig"), id="reference_metrics_scores_div", style={'display': 'none'})
            ],
            style={"display": "block", "flex-wrap": "wrap"}
        )
    ],
    color="light",
    style={"margin": "0rem 0rem 0.5rem 1rem", "height": "28rem"}
)


ref_confidence_dist_bar = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Confidence thresholds barplot"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("y-axis scaling:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.RadioItems(options=[{"label": "match count", "value": None},
                                                                {"label": "fraction total scans", "value": "normalize_psm"},
                                                                {"label": "matches per min", "value": "normalize_rt"}],
                                                        inline=True,
                                                        value=None, id="threshold_barplot_scaling"),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 3rem 0rem 0rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("db search %:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Checkbox(value=False, 
                                                    id="threshold_barplot_fill_bars"),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Filter de novo only:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Checkbox(value=False, 
                                                     id="threshold_barplot_filter_de_novo_only"),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                            ],
                            style={'display': 'flex', "align-items": "center"} 
                        ),
                    ],
                    style={'display': 'flex', "align-items": "center", "justify-content": "space-between"}
                ),
                html.Hr(),
                html.Div(dcc.Graph(id="reference_score_barplot_fig"), id="reference_score_barplot_div", style={'display': 'none'})
            ],
            style={"display": "block", "flex-wrap": "wrap"}
        )
    ],
    color="light",
    style={"margin": "1rem 0rem 0.5rem 1rem", "height": "30rem"}
)


################################################################################
# Page layout
################################################################################


ms_spectra_tab = [
    dbc.Row(
        id="qa_metric_values",
        className="d-flex justify-content-start",
        style={"margin": "0.5rem 0.5rem"}),
    tic_over_rt,
    dbc.Row([
        dbc.Col(
            feature_map,
            width=6,
            style={"padding": "0rem 0.5rem 0rem 0rem", "margin": "0rem"}    
        ),
        dbc.Col(
            transmission_scatter,
            width=6,
            style={"padding": "0rem 0rem 0rem 0.5rem", "margin": "0rem"}
        )
        # dbc.Col(
        #     ms2_spectrum
        # )
    ],
    style={"margin": "1rem 0.5rem 0.5rem 0.5rem", "padding": "0rem"}
    ),
    dbc.Row([
        dbc.Col(
            intensity_hist,
            width=6,
            style={"padding": "0rem 0.5rem 0rem 0rem", "margin": "0rem"}
        ),
        dbc.Col(
            charge_dist,
            width=6,
            style={"padding": "0rem 0.5rem 0rem 0.5rem", "margin": "0rem"}
        )
    ],
    style={"margin": "1rem 0.5rem 0.5rem 0.5rem", "padding": "0rem"}
    )
]


peptide_identification = [
    ident_score_rank_dist
]


ref_benchmark = [
    html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.P("Reference dataset:", className="me-4 mb-0 fs-6"),
                            dcc.Dropdown(
                                [],
                                id="ref_statistics_dropdown",
                                style={"width": "15rem"}
                            )
                        ],
                        className="d-flex align-items-center mx-2"
                    ),
                    html.B("or", className="mx-4 mb-0 fs-6"),
                    html.P("import custom dataset:", className="mx-2 mb-0 fs-6"),
                    dcc.Upload(
                        id="custom_ref_dataset",
                        children=dbc.Button("Import", outline=True, color="primary",
                                            className="px-4"),
                        className=" mx-1 py-1",
                    ),
                    html.I("No file...", id="custom_ref_statistics_name", className="ms-2"),
                    dcc.Store(id="current_ref_statistics_store", data=None),
                ],
                className="d-flex p-1 ms-4 align-items-center justify-content-start"
            ),
            html.Div(
                html.A(
                    html.I("Guide: how to create reference dataset", 
                           style={"font-size": "0.9vw"},
                    ),
                    href=GlobalConstants.docs_ref_data_prep_url, 
                    target="_blank"), 
                className="d-flex p-1 me-5 align-items-center justify-content-start"
            ),
        ],
        className="d-flex align-items-center justify-content-between"

    ),
    html.Div(
        [
            html.Div(
                [
                    html.P("DB search format: ", className="fst-italic me-3"),
                    html.P("...", id="ref_statistics_db_search_format", className="fw-bold"),
                ],
                className="d-flex ms-2 w-25"
            ),
            html.Div(
                [
                    html.P("De novo format: ", className="fst-italic me-3"),
                    html.P("...", id="ref_statistics_de_novo_format", className="fw-bold"),
                ],
                className="d-flex w-25"
            ),
        ],
        className="d-flex pt-1 ms-5 align-items-center justify-content-start"
    ),
    ref_metrics_comp,
    ref_confidence_dist_bar,
    ref_confidence_dist_line,
    dbc.Row(
        [
            dbc.Col(
                ref_confidence_dist_scatter,
                width=6
            ),
            dbc.Col(
                ref_intensity_dist_scatter,
                width=6
            )
        ],
    ),
    dbc.Row(
        [
            dbc.Col(
                ref_transmission_scatter,
                width="6",
                style={"margin": "0rem 0rem 0.5rem 0rem"}
            ),
            dbc.Col(
                ref_miscleavage_dist,
                width=6,
                style={"margin": "0rem 0rem 0.5rem 0rem"}
            )
        ],
    )
]


# Protein DB page
ms_performance = html.Div(
    [
        html.Div(
            [
                html.Div(
                    mzml_import,
                    id='mzml_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=StyleConstants.qa_import_box_style
                ),
                html.Div(
                    features_import,
                    id='features_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=StyleConstants.qa_import_box_style
                ),
                html.Div(
                    psm_file_validation_import,
                    id='db_search_psm_qa_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=StyleConstants.qa_import_box_style
                ),
                html.Div(
                    denovo_file_validation_import,
                    id='denovo_qa_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=StyleConstants.qa_import_box_style
                ),
            ],
            className="d-flex mb-3 justify-content-start"
        ),
        dbc.Alert(
            id="qa_data_import_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Tabs(
            [            
                dbc.Tab(
                    ref_benchmark,
                    label="Experiment comparison",
                    id="Reference Benchmark"
                ),
                # dbc.Tab(
                #     peptide_identification,
                #     label="Identification performance",
                #     id="Peptide Identification"
                # ),
                dbc.Tab(
                    ms_spectra_tab,
                    label="Experimental quality",
                    id="MS Spectra",
                    style={"margin": "0.5rem 0.5rem"}
                )
            ],
            id="ms_performance_tabs",
            # active_tab="Reference Benchmark"
        ),
    ],
    style={"margin-left": "0rem", "padding": "0rem 1rem"},
)






Protein_db = [
    dbc.Card(
        [
            dbc.CardHeader("Annotation Pattern"),
            dbc.CardBody(
                [
                    html.Div(html.P("Barplot figure"), id='denovo_fraction_barplot_graph'),
                    dbc.Row(
                        [
                            dbc.Col(html.P("y-axis scale: ", style={"margin-top": "0.5rem"}), width=1),
                            dbc.Col(
                                dbc.RadioItems(
                                    options=[
                                        {"label": "Absolute", "value": False},
                                        {"label": "Fraction", "value": True},
                                    ],
                                    value=False,
                                    id="de_novo_fraction_radio",
                                    inline=True,
                                    style={"width": "15rem", "margin-top": "0.5rem"}
                                ),
                                width=2
                            ),
                            dbc.Col(html.P("ALC cutoff: ", style={"margin-top": "0.5rem"}), width=1),
                            dbc.Col(
                                dbc.Input(value=85, 
                                          debounce=True,
                                          id="de_novo_fraction_alc_cutoff"),
                                width=1
                            )
                        ],
                    ),
                ]
            )  
        ],
        color="light",
        style={"margin": "5rem 0rem 1rem 1rem", "height": "35rem"}
    ),    
]
