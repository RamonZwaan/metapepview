from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from metapepview.layout.taxonomy_page import taxonomy_barplot, clade_filter



id_display_format = dbc.Row(
    [
    dbc.Col(html.B("KEGG ID display format:")),
    dbc.Col(
        dbc.RadioItems(
            options=[
                "KO",
                "EC",
                "Module",
                "Protein/Gene name"
            ],
            value="Protein/Gene name",
            id="kegg_display_format_radio",
            inline=False,
            style={"width": "15rem"}
        ),
    )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)

kegg_group_option_selector = dbc.Row(
    [
    dbc.Col(html.B("Select groups:")),
    dbc.Col(
        dbc.RadioItems(
            options=[
                "Pathway",
                "BRITE",
                "Manual"
            ],
            value="Pathway",
            id="kegg_group_type",
            style={"width": "15rem"}
        ),
    )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)


brite_selector = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.B("Select BRITE group:")),
            ],
            style={"margin": "0rem 0rem 1rem 0rem"}
        ),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="brite_group_dropdown",
                        placeholder="Select BRITE group",
                        disabled=False,
                        className="w-100",
                    ),
                )
            ],
            style={"margin": "0rem 0rem 1.5rem 0rem"}
        )
    ],
    id="brite_group_dropdown_container",
    hidden=True
)


pathway_selector = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.B("Select pathway:")),
            ],
            style={"margin": "0rem 0rem 1rem 0rem"}
        ),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        id="pathway_dropdown",
                        placeholder="Select pathway",
                        disabled=False,
                        className="w-100",
                    ),
                )
            ],
            style={"margin": "0rem 0rem 1.5rem 0rem"}
        )
    ],
    id="pathway_dropdown_container",
    hidden=True
)


module_selector = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.B("Select module:")),
            ],
            style={"margin": "0rem 0rem 1rem 0rem"}
        ),
        dbc.Row(
            [
                dbc.Col(
                    dcc.Dropdown(
                        [],
                        id="kegg_module_dropdown",
                        placeholder="Select module",
                        disabled=False,
                        className="w-100",
                    ),
                )
            ],
            style={"margin": "0rem 0rem 1.5rem 0rem"}
        )
    ],
    id="kegg_module_container",
    hidden=True
)


custom_protein_selector = html.Div(
    [
        dbc.Row(
            [
                dbc.Col(html.B("Custom proteins:")),
                dbc.Col(
                    dcc.Dropdown(
                        options=[],
                        multi=True,
                        id="custom_pathway_items",
                        disabled=False,
                        style={"width": "15rem"}
                    ),
                )
            ],
            style={"margin": "1.5rem 0rem"}
        )
    ],    
    id="custom_pathway_container",
    hidden=True
)

abundance_counting_selector = dbc.Row(
    [
        dbc.Col(html.B("Quantification")),
        dbc.Col(
            dbc.RadioItems(
                options=[
                    {"label": "PSM count", "value": "Match Count"}, 
                    {"label": "Combined signal intensity", "value":"Area"}
                ],
                value="Match Count",
                id="barplot_func_quantification_column",
                inline=False,
                style={"width": "15rem"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem", "display": "flex", "align-items": "center"}
)

normalize_samples = dbc.Row(
    [
        dbc.Col(html.B("Normalize sample:")),
        dbc.Col(
            dcc.Dropdown(
                id="pathway_normalize_sample",
                placeholder="Select...",
                disabled=False,
                style={"width": "15rem"}
            )
        )
    ],
    className="d-flex align-items-center",
    style={"margin": "1.5rem 0rem"}
)

combine_annot_duplicates = dbc.Row(
    [
        dbc.Col(html.B("Combine multiple annotations:")),
        dbc.Col(
            dbc.Checkbox(
                id="func_annot_combine_duplicates",
                value=False,
                style={"width": "5.5rem"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)

include_taxonomies = dbc.Row(
    [
        dbc.Col(html.B("Include taxonomies",
                       id="func_annot_incl_tax_text",
                       style={"text-decoration-line": "underline", 
                              "text-decoration-style": "dotted"})),
        dbc.Col(
            dbc.Checkbox(
                id="barplot_pathway_include_taxa_checkbox",
                value=True,
                style={"width": "5.5rem"}
            )
        ),
        dbc.Popover(
            """
                Show for each functional group which taxonomies express them.
                Bars in the graph are partitioned by Taxonomy name. The taxonomy 
                is shown by hovering over a bar element. 
            """,
            id="func_annot_incl_tax_popover",
            target="func_annot_incl_tax_text",
            trigger="hover",
            placement='right',
            className="p-2"
        )
    ],
    style={"margin": "1.5rem 0rem"}
)

fractional_abundances = dbc.Row(
    [
        dbc.Col(html.B("Fractional abundances")),
        dbc.Col(
            dbc.Checkbox(
                id="barplot_pathway_fraction_checkbox",
                value=False,
                style={"width": "5.5rem"}
            )
        )
    ],
    style={"margin": "1.5rem 0rem"}
)


kegg_export_block = [
    html.H5("Visualize pathway map"),#, style={"margin-top": "2.5rem"}),
    html.Hr(),
    dbc.Row(
        [
            dbc.Col(html.B("Select Samples (Max 4)")),
            dbc.Col(
                dcc.Dropdown(
                    options=[],
                    value=[],
                    id="kegg_export_samples",
                    style={"width": "15rem"},
                    multi=True
                )
            )
        ],
        className="d-flex align-items-center",
        style={"margin": "1.5rem 0rem"}
    ),
    html.Div(
        id="sample_color_table_func_export"
    ),
    dbc.Row(
        [
            dbc.Button("Show pathway map",
                        id="kegg_pathway_map_link",
                        target="_blank",
                        external_link=True)
        ],
        style={"margin": "1rem 5rem"}
    )
]

data_export_button = [
    html.H5("Export data"),
    html.Hr(),
    dbc.Row(
        [
            dbc.Button("Export complete functions",
                        id="export_functions_button",
                        className="")
        ],
        style={"margin": "2rem 5rem"}
    )
]


pathway_filters = [
    html.H3("Plot options"),
    html.Hr(),
    id_display_format,
    html.Hr(),
    kegg_group_option_selector,
    brite_selector,
    pathway_selector,
    module_selector,
    custom_protein_selector,
    html.Hr(),
    abundance_counting_selector,
    normalize_samples,
    combine_annot_duplicates,
    include_taxonomies,
    fractional_abundances,
] + clade_filter + data_export_button + kegg_export_block


pathway_filters_acc = [
    html.H3("Plot options"),
    html.Hr(),
    dbc.Accordion(
        [
            dbc.AccordionItem(
                html.Div(
                    [
                        id_display_format,
                    ]
                ),
                title="Select display format"
            ),
            dbc.AccordionItem(
                html.Div(
                    [
                        kegg_group_option_selector,
                        brite_selector,
                        pathway_selector,
                        module_selector,
                        custom_protein_selector,
                    ]
                ),
                title="Select function group"
            ),
            dbc.AccordionItem(
                html.Div(
                    [
                        abundance_counting_selector,
                        normalize_samples,
                        combine_annot_duplicates,
                        include_taxonomies,
                        fractional_abundances,
                    ]
                ),
                title="Configure quantification"
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
] + data_export_button + kegg_export_block


functional_annotation_barplot = [
    html.Div(
        [
            html.H5("Figure", id="func_annot_figure_title"),
            dbc.Button("Export figure data", 
                       id="export_func_annot_figure_data",
                       className="me-3",
                       disabled=True)
        ],
        className="d-flex align-items-center justify-content-between"
    ),
    html.Hr(className="my-2"),
    html.Div(dcc.Graph(id="pathway_barplot_figure"),
             id='pathway_barplot_graph', 
             style={"display":"None"}),
]


# Functional annotation page
functional_annotation_page = [
    dbc.Row(
        [
            dbc.Col(
                dbc.Card(
                    [
                        # dbc.CardHeader("Filter settings"),
                        dbc.CardBody(pathway_filters_acc)
                    ],
                    className="shadow-sm",
                    style={"margin": "0rem 0rem 1rem 1rem",
                           "height": "64rem",
                           "overflow-y": "scroll"}
                ),
                width={'size': 4}
            ),
            dbc.Col(
                dbc.Card(
                    [
                        dbc.CardBody(
                            functional_annotation_barplot
                        )
                    ],
                    className="shadow-sm",
                    style={"margin": "0rem 0rem 1rem 0rem", "height": "64rem"}
                ),
                width={'size': 8, 'order': 'last'}
            ),
            dcc.Download(id="download_community_functions_csv"),
        ],
        style={"height": "100%"}
    ),
]
