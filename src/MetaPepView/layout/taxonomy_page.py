from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from .style_constants import *
from constants import *


# buttongroup to select the type of graph to show
graph_type_selection = dbc.Row(
    [
        dbc.ButtonGroup(
            [
                dbc.Button("Stacked barplot", id="taxonomy_stacked_barplot_button", active=True, color='secondary', outline=True, style={'width': '50%'}),
                dbc.Button("Heatmap", id="taxonomy_heatmap_button", active=False, color='secondary', outline=True, style={'width': '50%'})
            ],
            style={'display': 'flex'}
        )
    ],
    style={"display": 'flex', 'justify-content': 'center'}
)


sample_selector = dbc.Row(
    [
        dbc.Col(html.B("Select Sample")),
        dbc.Col(
            dcc.Dropdown(
                options=[],
                value=None,
                id="barplot_de_novo_sample_items",
                style={"width": "15rem"}
            )
        )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)



# dropdown menu to select taxonomy names to display
taxa_dropdown_selector = dbc.Row(
    [
        dbc.Col(html.B("Select taxa")),
        dbc.Col(
            dcc.Dropdown(
                options=[],
                value=[],
                multi=True,
                id="barplot_custom_taxa_items",
                style={"width": "15rem"}
            )
        )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)


taxa_group_display_selector = dbc.Row(
    [
        dbc.Col(html.B("Taxa display type")),
        dbc.Col(
            dbc.RadioItems(
                options=[
                    {"label": "Top 10", "value": 1},
                    {"label": "Top 25", "value": 3},
                    {"label": "Custom taxa", "value": 2},
                ],
                value=1,
                id="barplot_taxa_selector_radio",
                inline=False,
                style={"width": "15rem"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem", "display": "flex", "align-items": "center"}
)


taxonomy_rank_selector = dbc.Row(
    [
        dbc.Col(html.B("Display rank")),
        dbc.Col(
            dcc.Dropdown(
                GlobalConstants.standard_lineage_ranks,
                value='Phylum',
                id='barplot_taxa_rank_items',
                style={"width": "15rem"}
            )
        )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)


normalize_abundance_selector = dbc.Row(
    [
        dbc.Col(html.B("Normalize abundances")),
        dbc.Col(
            dbc.Checkbox(
                id="barplot_taxa_fraction_checkbox",
                value=False,
                style={"width": "15rem"}#, "display": "flex", "justify-content": "flex-end"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)


include_unannotated_selector = dbc.Row(
    [
        dbc.Col(html.B("Include unannotated")),
        dbc.Col(
            dbc.Checkbox(
                id="barplot_taxa_unannotated_checkbox",
                value=False,
                style={"width": "15rem"}#, "display": "flex", "justify-content": "flex-end"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)

global_fallback_selector = dbc.Row(
    [
        dbc.Col(html.B("Allow global annotation")),
        dbc.Col(
            dbc.Checkbox(
                id="barplot_taxa_allow_global_annot_checkbox",
                value=False,
                style={"width": "15rem"}#, "display": "flex", "justify-content": "flex-end"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)

global_annot_de_novo_only_selector = dbc.Row(
    [
        dbc.Col(html.B("De novo only global annotation")),
        dbc.Col(
            dbc.Checkbox(
                id="global_annot_de_novo_only_checkbox",
                value=False,
                style={"width": "15rem"}#, "display": "flex", "justify-content": "flex-end"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)


# modal body of barplot
tax_sample_comp_options = [
    html.H3("Plot options"),
    html.Hr(),
    graph_type_selection,
    taxa_dropdown_selector,
    taxa_group_display_selector,
    taxonomy_rank_selector,
    normalize_abundance_selector,
    include_unannotated_selector,
    global_fallback_selector
]


tax_sample_de_novo_options = [
    html.H3("Plot options"),
    html.Hr(),
    graph_type_selection,
    sample_selector,
    taxa_dropdown_selector,
    taxa_group_display_selector,
    taxonomy_rank_selector,
    global_annot_de_novo_only_selector,
    normalize_abundance_selector,
    include_unannotated_selector
]


clade_filter = [
    html.H5("Filter by Clade", style={"margin-top": "2.5rem"}),
    html.Hr(),
    dbc.Row(
        [
            dbc.Col(html.B("Clade rank"), style={"text-align": 'left'}),
            dbc.Col(dcc.Dropdown(
                ['Root',
                 'Superkingdom',
                 'Phylum',
                 'Class',
                 'Order',
                 'Family',
                 'Genus',
                 'Species'],
                id='tax_barplot_clade_selection_rank',
                style={"width": "15rem"}
                )
            )
        ],
        className="d-flex align-items-center",
        style={"margin": "1.5rem 0rem"}
    ),
    dbc.Row(
        [
            dbc.Col(html.B("Root taxonomy"), style={"text-align": 'left'}),
            dbc.Col(dcc.Dropdown(
                id='tax_barplot_clade_selection_taxa',
                disabled=True,
                style={"width": "15rem"}
                )
            ),
        ],
        className="d-flex align-items-center",
        style={"margin": "1.5rem 0rem"}
    ),

]


export_button = [
    dbc.Row(
        [
            dbc.Button("Export taxonomy",
                        id="export_taxonomy_button",
                        className="")
        ],
        style={"margin": "2.5rem 5rem"}
    )
]


# contents of box that displays taxonomy barplot
taxonomy_barplot = [
    html.H3("Figure", id="taxonomy_figure_title"),
    html.Hr(),
    html.Div(dcc.Graph(id="taxonomy_barplot_figure"), id='taxa_barplot_graph', style={'display': 'None'}),
]


taxonomy_de_novo_barplot = [
    html.H3("Figure", id="taxonomy_de_novo_figure_title"),
    html.Hr(),
    html.Div([
        dcc.Graph(id="taxonomy_barplot_de_novo_figure"),
        dcc.Graph(id="taxonomy_dif_barplot_de_novo_figure"),
    ], id='taxa_barplot_de_novo_graph', style={'display': 'None'}),
]


taxonomic_dropoff_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("taxonomic dropoff", id="taxonomic_dropoff_title")),
        dbc.ModalBody(
            [
                html.Div(
                    [
                        html.P("Cumulative annotation drop from:",
                                className="mb-0 me-3 fs-5 text"),
                        dcc.Dropdown(
                            ["Root"] + GlobalConstants.standard_lineage_ranks,
                            value="Root",
                            clearable=False,
                            id='taxonomy_cumulative_dropoff_rank',
                            style={"width": "20rem"}),
                    ],
                    className="d-flex align-items-center justify-content-between mb-3 mx-4",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.P("Lineage drop:",
                                       className="mb-0 fs-5 text"),
                                html.P(id="lineage_tax_dropoff_text",
                                       className="mb-0 ms-1 fs-5 text fw-bold")
                            ],
                            className="d-flex ms-3"
                        ),
                        html.Div(
                            [
                                html.P("Clade average drop:",
                                       className="mb-0 fs-5 text"),
                                html.P(id="global_av_dropoff_text",
                                       className="mb-0 ms-1 fs-5 text fw-bold"),
                            ],
                            className="d-flex ms-3"
                        )
                    ],
                    className="d-flex align-items-center justify-content-between mb-3 mx-4",
                ),
                html.Hr(className="p-0"),
                dbc.Checkbox(id="taxonomic_dropoff_normalize",
                             label="abundance distribution percentage",
                             value=False,
                             className="ms-4"),
                html.Div(
                    [
                        dcc.Graph(id="taxonomic_dropoff_figure"),
                    ],
                    id='taxonomic_dropoff_graph', style={'display': 'None'})
            ], 
            className="px-1 pb-1")
    ],
    id="taxonomic_dropoff_modal",
    size="lg",
    is_open=False
)


# TODO: Move function to html templates
def taxonomy_page_constructor(
    filter_settings,
    graph,
    description,
    hidden_components = []
):
    return [
    dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    [
                        dbc.CardBody(filter_settings)
                    ],
                    className="shadow-sm",
                    style={"margin": "0rem 0rem 1rem 1rem",
                           "height": "64rem",
                           "overflow-y": "scroll"}
                ),
                width={'size': 4}
            ),
            dbc.Col(
                dbc.Card( # stacked barplot card
                    [
                        dbc.CardBody(graph)    
                    ],
                    className="shadow-sm",
                    style={"margin": "0rem 0rem 1rem 0rem", "height": "64rem"}
                ),
                width={'size': 8, 'order': 'last'}
            ),
        ],
        style={"height": "100%"}
    ),
    ] + hidden_components


# Taxonomy annotation page
taxonomy_sample_analysis = taxonomy_page_constructor(
    tax_sample_comp_options + clade_filter + export_button,
    taxonomy_barplot,
    [
        html.P("...")
    ],
    [
        # download component for taxonomy export
        dcc.Download(id="download_taxonomy_composition_csv"),
        taxonomic_dropoff_modal
    ]
)

# Taxonomy annotation page
taxonomy_de_novo_analysis = taxonomy_page_constructor(
    tax_sample_de_novo_options + clade_filter,
    taxonomy_de_novo_barplot,
    [
        html.P("...")
    ],
    [
        # download component for taxonomy export
        dcc.Download(id="download_taxonomy_composition_csv")
    ]
)

