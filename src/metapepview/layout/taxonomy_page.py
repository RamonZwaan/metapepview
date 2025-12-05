from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from metapepview.constants import GlobalConstants as gc

from metapepview.html_templates import abundance_counting_selector,\
    normalize_abundance_selector,\
    include_unannotated_selector,\
    global_fallback_selector


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
        dbc.Col(html.B("Display taxa")),
        dbc.Col(
            dbc.RadioItems(
                options=[
                    {"label": "Top 10 abundant", "value": 1},
                    {"label": "Top 20 abundant", "value": 3},
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
                gc.standard_lineage_ranks,
                value='Phylum',
                id='barplot_taxa_rank_items',
                style={"width": "15rem"}
            )
        )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)


global_annot_de_novo_only_selector = dbc.Row(
    [
        dbc.Col(html.B("Unipept composition de novo only",
                       id="unipept_de_novo_only_text",
                       style={"text-decoration-line": "underline", 
                              "text-decoration-style": "dotted"})),
        dbc.Col(
            dbc.Checkbox(
                id="global_annot_de_novo_only_checkbox",
                value=False,
                style={"width": "15rem"}#, "display": "flex", "justify-content": "flex-end"}
            )
        ),
        dbc.Popover(
            """
                For the community composition from Unipept taxonomy, only use de novo
                peptides. If unchecked, all DB search and de novo peptides with Unipept
                annotation are included in the Unipept based community composition.
                Setting this option eliminates quantification biases that arise due to 
                over- or under-representation of species in the local database, which
                are reflected in DB search results. However, filtering out DB search
                peptides strongly reduces the number of confident peptide annotations.
            """,
            id="glob_annot_de_novo_popover",
            target="unipept_de_novo_only_text",
            trigger="hover",
            placement='right',
            className="p-2"
        )
    ],
    style={"margin": "1.5rem 0rem"}
)


facet_plot_switch = dbc.Row(
    [
        dbc.Col(html.B("Enable side-by-side plot")),
        dbc.Col(
            dbc.Switch(id="activate_taxonomy_facet", 
                       value=False,
                       className="align-self-center",
                       style={"width": "15rem"})
        )
    ],
    style={"margin": "1.5rem 0rem"}
)


tax_sample_facet_collapse = [
    facet_plot_switch,
    dbc.Collapse(
        dbc.Card(
            dbc.CardBody(
                [
                    html.H5("Facet plot options"),
                    abundance_counting_selector(True),
                    normalize_abundance_selector(True),
                    include_unannotated_selector(True),
                    global_fallback_selector(True)
                ]
            )
        ),
        id="tax_comp_facet_options",
        is_open=False
    )
]

clade_filter = [
    # html.H5("Filter by clade", style={"margin-top": "2.5rem"}),
    # html.Hr(),
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


# modal body of barplot
tax_sample_comp_options = [
    html.H3("Plot options"),
    html.Hr(),
    graph_type_selection,
    taxa_dropdown_selector,
    taxa_group_display_selector,
    taxonomy_rank_selector,
    abundance_counting_selector(),
    normalize_abundance_selector(),
    include_unannotated_selector(),
    global_fallback_selector()
]


tax_sample_comp_options_acc = [
    html.H3("Plot options"),
    html.Hr(),
    graph_type_selection,
    dbc.Accordion(
        [
            dbc.AccordionItem(
                html.Div(
                    [
                        taxa_dropdown_selector,
                        taxa_group_display_selector,
                        taxonomy_rank_selector,
                    ]
                ),
                title="Select taxa"
            ),
            dbc.AccordionItem(
                html.Div(
                    [
                        abundance_counting_selector(),
                        normalize_abundance_selector(),
                        include_unannotated_selector(),
                        global_fallback_selector()
                    ]
                ),
                title="Configure quantification"
            ),
            dbc.AccordionItem(
                html.Div(
                    [
                        facet_plot_switch,
                        dbc.Collapse(
                            dbc.Card(
                                dbc.CardBody(
                                    [
                                        html.H5("Facet plot options"),
                                        abundance_counting_selector(True),
                                        normalize_abundance_selector(True),
                                        include_unannotated_selector(True),
                                        global_fallback_selector(True)
                                    ]
                                )
                            ),
                            id="tax_comp_facet_options",
                            is_open=False
                        )
                    ]
                ),
                title="Configure secondary plot"
            ),
            dbc.AccordionItem(
                html.Div(
                    clade_filter,
                ),
                title="Configure taxonomy clade filter"
            )
        ],
        always_open=True,
        start_collapsed=True,
        className="my-3"
    )
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
    abundance_counting_selector(),
    normalize_abundance_selector(),
    include_unannotated_selector()
]


tax_sample_de_novo_options_acc = [
    html.H3("Plot options"),
    html.Hr(),
    graph_type_selection,
    sample_selector,

    dbc.Accordion(
        [
            dbc.AccordionItem(
                html.Div(
                    [
                        taxa_dropdown_selector,
                        taxa_group_display_selector,
                        taxonomy_rank_selector,
                    ]
                ),
                title="Select taxa"
            ),
            dbc.AccordionItem(
                html.Div(
                    [
                        global_annot_de_novo_only_selector,
                        abundance_counting_selector(),
                        normalize_abundance_selector(),
                        include_unannotated_selector()
                    ]
                ),
                title="Configure quantification"
            ),
        ],
        always_open=True,
        start_collapsed=True,
        className="my-3"
    )
]


export_button = [
    html.H5("Export data"),
    html.Hr(),
    dbc.Row(
        [
            dbc.Button("Export complete compositions",
                        id="export_taxonomy_button",
                        className="")
        ],
        style={"margin": "2rem 5rem"}
    )
]

comp_eval_export_button = [
    html.H5("Export data"),
    html.Hr(),
    dbc.Row(
        [
            dbc.Button("Export complete compositions",
                        id="export_taxonomy_eval_button",
                        className="")
        ],
        style={"margin": "2rem 5rem"}
    )
]

# contents of box that displays taxonomy barplot
taxonomy_barplot = [
    html.Div(
        [
            html.H5("Figure", id="taxonomy_figure_title"),
            dbc.Button("Export figure data", 
                       id="export_taxonomy_figure_data",
                       className="me-3",
                       disabled=True)
        ],
        className="d-flex align-items-center justify-content-between"
    ),
    html.Hr(className="my-2"),
    html.Div(dcc.Graph(id="taxonomy_barplot_figure"), id='taxa_barplot_graph', style={'display': 'None'}),
    dcc.Store(id="taxonomy_barplot_figure_data", data=None)
]


taxonomy_de_novo_barplot = [
    html.Div(
        [
            html.H5("Figure", id="taxonomy_de_novo_figure_title"),
            dbc.Button("Export figure data", 
                       id="export_taxonomy_de_novo_figure_data",
                       className="me-3",
                       disabled=True)
        ],
        className="d-flex align-items-center justify-content-between"
    ),

    html.Hr(className="my-2"),
    html.Div([
        dcc.Graph(id="taxonomy_barplot_de_novo_figure"),
        dcc.Graph(id="taxonomy_dif_barplot_de_novo_figure"),
    ], id='taxa_barplot_de_novo_graph', style={'display': 'None'}),
    dcc.Store(id="taxonomy_de_novo_figure_data")
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
                            ["Root"] + gc.standard_lineage_ranks,
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
    # tax_sample_comp_options + tax_sample_facet_collapse + clade_filter + export_button,
    tax_sample_comp_options_acc + export_button,
    taxonomy_barplot,
    [
        html.P("...")
    ],
    [
        # download component for taxonomy export
        dcc.Download(id="download_taxonomy_composition_csv"),
        dcc.Download(id="download_taxonomy_figure_data_csv"),
        taxonomic_dropoff_modal
    ]
)

# Taxonomy annotation page
taxonomy_de_novo_analysis = taxonomy_page_constructor(
    # tax_sample_de_novo_options + comp_eval_export_button, # + clade_filter,
    tax_sample_de_novo_options_acc + comp_eval_export_button,
    taxonomy_de_novo_barplot,
    [
        html.P("...")
    ],
    [
        # download component for taxonomy export
        dcc.Download(id="download_taxonomy_composition_de_novo_csv")
    ]
)
