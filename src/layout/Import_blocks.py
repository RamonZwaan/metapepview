from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from .style_constants import *


peptide_dataset = [
    html.H5("Peptides Dataset"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                   "margin-bottom": ".5rem"}
                    )
                    ],
                    id="peptides_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.Div(
                    [
                    html.P("No file...", style={"padding": "0rem .5rem",
                                                "text-align": "left"})
                    ],
                    id="peptides_name"
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    )
]

mzxml_spectra = [
    html.H5("Spectra (mzXML)"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                   "margin-bottom": ".5rem"}
                    )
                    ],
                    id="mzml_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.Div(
                    [
                    html.P("No file...", style={"padding": "0rem .5rem",
                                                "text-align": "left"})
                    ],
                    id="mzml_name", 
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    )
]

db_search_psm = [
    html.H5("db search PSM"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                   "margin-bottom": ".5rem"}
                    )
                    ],
                    multiple=True,
                    id="db_search_psm_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.Div(
                    [
                    html.P("No file...", style={"padding": "0rem .5rem",
                                                "text-align": "left"})
                    ],
                    id="db_search_psm_name", 
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    )
]


de_novo_dataset = [
    html.H5("De Novo Dataset"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                   "margin-bottom": ".5rem"}
                    )
                    ],
                    multiple=True,
                    id="denovo_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.Div(
                    html.P("No file...", style={"padding": "0rem .5rem",
                                                "text-align": "left"}
                    ),
                    id="denovo_name"
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    )
]


taxonomy_db = [
    html.H5("Taxonomy DB"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                "margin-bottom": ".5rem"}
                    )
                    ],
                    id="taxonomy_db_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.P("No file...", id="taxonomy_db_name", style={"padding": "0rem .5rem",
                                                                    "text-align": "left"}
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    ),
    html.H6("File format:", style={"font-weight": 'bold'}),
    dbc.RadioItems(
        options=[
            {"label": "fasta", "value": "fasta"},
            {"label": "accession-taxonomy map", "value": "accession-taxonomy map"},
            {"label": "sequence-taxonomy map", "value": "sequence-taxonomy map"},
        ],
        value=1,
        id="taxonomy_db_format_radio",
        inline=False,
    )
]


functional_annotation = [
    html.H5("Functional Annotation"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                "margin-bottom": ".5rem"}
                    )
                    ],
                    id="func_annot_db_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.P("No file...", id="func_annot_name", style={"padding": "0rem .5rem",
                                                                "text-align": "left"}
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    ),
    html.H6("Source:", style={"font-weight": 'bold'}),
    dbc.RadioItems(
        options=[
            {"label": "KEGG", "value": 'KEGG'},
            {"label": "EggNOG", "value": 'EggNOG'},
        ],
        value=1,
        id="func_annot_db_format_radio",
        inline=False,
    )
]


ref_set_files = [
    html.H5("Reference set (DB Search + De Novo)"),
    dbc.Row(
        [
            dbc.Col(
                dcc.Upload(
                    [
                    dbc.Button("browse", color="secondary", style={"margin-top": ".5rem",
                                                                   "margin-bottom": ".5rem"}
                    )
                    ],
                    multiple=True,
                    id="reference_files_upload"
                ),
                width=5
            ),
            dbc.Col(
                html.Div(
                    [
                    html.P("No file...", style={"padding": "0rem .5rem",
                                                "text-align": "left"})
                    ],
                    id="reference_files_name", 
                )
            )
        ],
        style={'display': "flex", "align-items": "center"}
    )
]