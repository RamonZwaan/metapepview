from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from backend.html_templates import qa_importer_block

from .style_constants import *
from constants import GlobalConstants




mzml_import = qa_importer_block(
    "Spectral File (MzML)",
    "mzml_upload",
    "mzml_name",
    "mzml_valid",
)



psm_file_validation_import = qa_importer_block(
    "PSM Import",
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
                        html.H4("Total Ion Current", style={"margin-right": "3rem"}),
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
                                        html.B("SMA window:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=1, max=1000, step=1, value=20, id="tic_ms_sma_range", style={'width': '6rem'}),       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Data reduction factor:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=1, max=10, step=1, value=3, id="tic_ms_red_fact", style={'width': '4rem'}),       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("Secondary param:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dcc.Dropdown(
                                            ["None",
                                             "DB Search Confidence",
                                             "De Novo Confidence",
                                             "Peak Count"],
                                            value="None",
                                            id="tic_ms_secondary_y",
                                            style={"width": "12rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("int. cutoff:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=0, step=1, value=0, id="tic_sec_param_int_cutoff", style={'width': '8rem'}),       
                                    ],
                                    id="tic_sec_param_int_cutoff_container",
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"},
                                    hidden=True
                                ),
                                html.Div(
                                    [
                                        html.B("Confidence cutoff:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=1, step=1, value=0, id="tic_sec_param_conf_cutoff", style={'width': '6rem'}),       
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
    style={"margin": "1rem 0rem 0.5rem 1rem", "height": "30rem"}
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
                                        html.B("int. cutoff:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Input(type="number", min=0, max=1e8, step=1, value=0, id="mz_over_rt_int_cutoff", style={'width': '8rem'}),
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
    style={"margin": ".5rem 0rem 0rem 1rem", "height": "25rem"}
)


intensity_hist = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Scan Intensities"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("MS level:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Select(options=[{"label": "MS1", "value": 1}, {"label": "MS2", "value": 2}],
                                                    value=1, id="scan_int_dist_ms_level", style={"width": "6rem"}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("ALC cutoff (de novo):", style={"margin-right": "1rem", "margin-left": "1rem", "text-align": "center"}),
                                        dbc.Input(type="number", min=0, max=100, step=1, value=80, disabled=True, id="scan_int_dist_alc_cutoff", style={'width': '6rem'}),
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                )
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
    style={"margin": "1rem 0rem 1rem 1rem", "height": "25rem"}
)

transmission_scatter = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("MS2 Efficiency"),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.B("min m/z:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Input(type="number", min=0, max=1400, step=1, value=0, id="ms1_over_ms2_mz_cutoff", style={'width': '6rem'}),
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
    style={"margin": ".5rem 0rem .5rem 0rem", "height": "25rem"}
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
                                        html.B("-10lgP cutoff:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=0, max=200, step=1, value=35, id="pept_conf_lgp_cutoff", style={'width': '6rem'}),       
                                    ],
                                    style={'display': 'flex', "align-items": "center", "margin": "0rem 0.5rem"}
                                ),
                                html.Div(
                                    [
                                        html.B("ALC cutoff:", style={"text-align": "center", "margin-right": "1rem", "margin-left": "1rem"}),
                                        dbc.Input(type="number", min=0, max=200, step=1, value=35, id="pept_conf_alc_cutoff", style={'width': '6rem'}),       
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
    style={"margin": "1rem 0rem 0.5rem 1rem", "height": "35rem"}
)


ref_confidence_dist_line = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Score Confidence Distribution"),
                        html.Div(
                            [
                                html.B("Metric:", style={"margin-right": "1rem", "text-align": "center"}),
                                dbc.Select(options=[{"label": "DB search", "value": "db search"},
                                                    {"label": "De novo", "value": "de novo"}],
                                            value="db search", id="ref_score_conf_metric", style={"width": "14rem"}),
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
                        html.H4("Score Threshold Distribution"),
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
                html.H4("Scan Intensity Distribution"),
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
                        html.H4("Ion transmission loss"),
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



ref_confidence_dist_bar = dbc.Card(
    [
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.H4("Spectral Matches"),
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
                                        html.B("fill bars:", style={"margin-right": "1rem", "text-align": "center"}),
                                        dbc.Checkbox(value=False, 
                                                    id="threshold_barplot_fill_bars"),
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
    tic_over_rt,
    dbc.Row([
        dbc.Col(
            feature_map,
            width=6
        ),
        dbc.Col(
            transmission_scatter,
            width=6
        )
        # dbc.Col(
        #     ms2_spectrum
        # )
    ]
    ),
    dbc.Row([
        dbc.Col(
            intensity_hist,
            width=6
        )
    ]
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
                    html.P("Reference dataset:", className="mx-4 mb-0 fs-6"),
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
            dcc.Store(id="current_ref_statistics_store", data=None)
            #html.H5("Add own statistics file")
        ],
        className="d-flex p-1 align-items-center justify-content-start"
    ),
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
                width="6"
            )
        ],
    ),
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
                    style=qa_import_box_style
                ),
                html.Div(
                    psm_file_validation_import,
                    id='db_search_psm_qa_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=qa_import_box_style
                ),
                html.Div(
                    denovo_file_validation_import,
                    id='denovo_qa_import_box',
                    className="ms-3 mw-25 pt-2 border shadow-sm", 
                    style=qa_import_box_style
                ),
                # html.Div(
                #     ref_set_import,
                #     className="ms-3 mw-25 pt-2 border shadow-sm", 
                #     style={"padding": ".0rem .5rem", "background-color": PLOT_COLOR, "border-radius": ".5rem"}
                # ),
            ],
            className="d-flex mb-3 justify-content-start"
        ),
        dbc.Tabs(
            [            
                dbc.Tab(
                    ms_spectra_tab,
                    label="Identification performance",
                    id="MS Spectra"
                ),
                dbc.Tab(
                    peptide_identification,
                    label="Experimental quality",
                    id="Peptide Identification"
                ),
                dbc.Tab(
                    ref_benchmark,
                    label="Reference Benchmark",
                    id="Reference Benchmark"
                )
            ],
            id="ms_performance_tabs",
            active_tab="MS Spectra",
            style={}
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
                                dbc.Input(value=85, id="de_novo_fraction_alc_cutoff"),
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
