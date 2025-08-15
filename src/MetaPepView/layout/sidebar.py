from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from constants import GlobalConstants as gc

"""
Block elements for sidebar to display
"""


fetch_db_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Public DB Settings")),
        dbc.ModalBody(
            [
                html.H4("Database Location", className="mb-3"),
                dbc.Label("NCBI Taxonomy path", className="fst-italic"),
                dbc.Input(id='ncbi_taxonomy_db_loc',
                    value=gc.ncbi_taxonomy_dir,
                    size="sm mb-4"),

                # gtdb module will only be provided in full functionality mode
                html.Div(
                    [
                        dbc.Label("GTDB Taxonomy path", className="fst-italic"),
                        dbc.Input(id='gtdb_taxonomy_db_loc',
                            value=gc.gtdb_taxonomy_dir,
                            size="sm mb-4"),
                    ],
                    className="" if gc.show_advanced_settings is True else "d-none"
                ),
                # KEGG module will be hidden in "de novo only" mode
                html.Div(
                    [
                        dbc.Label("KEGG KO map path", className="fst-italic"),
                        dbc.Input(id='kegg_map_loc', value=gc.kegg_map_dir, size="sm mb-4"),
                    ],
                    className="" if gc.display_db_search is True else "d-none"
                ),
                html.Hr(),
                html.H4("Fetch public database", className="mb-3"),
                dbc.Alert(
                    id="db_download_status_alert",
                    dismissable=True,
                    is_open=False,
                    color="danger"
                ),
                dbc.Checkbox(
                    id="fetch_ncbi_taxonomy_checkbox",
                    label="NCBI Taxonomy",
                    value=False
                ),
                html.Div(
                    [
                        dbc.Checkbox(
                            id="ncbi_taxonomy_overwrite_old_checkbox",
                            label="Overwrite existing data",
                            value=False
                        ),
                        dbc.Checkbox(
                            id="ncbi_taxonomy_create_parent_dirs_checkbox",
                            label="Create parent dir",
                            value=False,
                            className="ms-5"
                        ),
                    ],
                    className="d-flex"
                ),
                dbc.Label("Source URL", className="fst-italic"),
                dbc.Input(id='ncbi_taxonomy_db_source_url', value=gc.ncbi_taxonomy_url, size="sm"),
                dbc.FormText("Database stored in path location specified above.", className="fst-italic"),

                # gtdb download module only shown in full functional interface
                html.Div(
                    [
                        dbc.Checkbox(
                            id="fetch_gtdb_taxonomy_checkbox",
                            label="GTDB Taxonomy",
                            value=False,
                            className="mt-4"
                        ),
                        html.Div(
                            [
                                dbc.Checkbox(
                                    id="gtdb_taxonomy_overwrite_old_checkbox",
                                    label="Overwrite existing data",
                                    value=False
                                ),
                                dbc.Checkbox(
                                    id="gtdb_taxonomy_create_parent_dirs_checkbox",
                                    label="Create parent dir",
                                    value=False,
                                    className="ms-5"
                                ),
                            ],
                            className="d-flex"
                        ),
                        dbc.Label("Source URL (Parent Dir)", className="fst-italic"),
                        dbc.Input(id='gtdb_taxonomy_db_source_url', value=gc.gtdb_taxonomy_url, size="sm"),
                        dbc.FormText("Database stored in path location specified above.", className="fst-italic"),
                    ],
                    className="" if gc.show_advanced_settings is True else "d-none"
                ),

                # KEGG download module hidden in "de novo only" mode
                html.Div(
                    [
                        dbc.Checkbox(
                            id="fetch_kegg_map_checkbox",
                            label="KEGG KO mapping",
                            value=False,
                            className="mt-4"
                        ),
                        html.Div(
                            [
                                dbc.Checkbox(
                                    id="kegg_map_overwrite_old_checkbox",
                                    label="Overwrite existing data",
                                    value=False
                                ),
                                dbc.Checkbox(
                                    id="kegg_map_create_parent_dirs_checkbox",
                                    label="Create parent dir",
                                    value=False,
                                    className="ms-5"
                                ),
                            ],
                            className="d-flex"
                        ),
                        dbc.FormText("Database stored in path location specified above.", className="fst-italic mb-xs"),
                    ],
                    className="" if gc.display_db_search is True else "d-none"
                ),

                html.Div(
                    [
                        dbc.Spinner(html.Div(id="loader_fetch_db", style={"width": "3rem"}), size="sm"),
                        dbc.Button("Fetch data", id="start_database_download", className="", n_clicks=0),
                        dcc.Store(id="data_imported_trigger", data=False),
                    ],
                    className="ms-auto d-flex"
                )

            ],
            className="vstack"
        )
    ],
    id="fetch_database_modal",
    scrollable=True,
    is_open=False,
)


project_modules = [
    dbc.NavLink(
        [
            html.I(className="bi bi-gear me-2"),
            "Create project"
        ],
        id="sidebar_data_button", href="/", active=True, className="mb-2"),
]

community_modules = [
    dbc.NavLink(
        [
            html.I(className="bi bi-graph-up me-2"),
            "Community composition"
        ],
        id="sidebar_taxonomy_button", href="/", className="mb-2"),
    html.Div(
        dbc.NavLink(
            [
                html.I(className="bi bi-graph-up me-2"),
                "Community functions"
            ],
            id="sidebar_functional_button",
            href="/",
            className="mb-2"),
        hidden = not gc.display_db_search),
     ]


evaluation_modules = [
    html.Div(
        dbc.NavLink(
            [
                html.I(className="bi bi-speedometer2 me-2"),
                "Experimental performance"
            ],
            id="sidebar_validation_button", href="/", className="mb-2"),
        hidden= not gc.display_qa_page
    ),
    html.Div(
        dbc.NavLink(
            [
                html.I(className="bi bi-clipboard-data me-2"),
                "Evaluate community composition"
            ],
            id="sidebar_taxonomy_de_novo_button", href="/", className="mb-2"),
        hidden= not gc.display_db_search or not gc.display_de_novo
    )
    ]

# filter out presented modules based on


new_sidebar = [
    html.Div(
        [
            html.H4("Menu", className="display-4"),
            html.Hr(),
            html.H3("Project management"),
            dbc.Nav(
                project_modules,
                vertical=True,
                pills=True,
                className="mb-4"
            ),
            html.H3("Community analysis"),
            dbc.Nav(
                community_modules,
                vertical=True,
                pills=True,
                className="mb-4"
            ),
            html.H3("Experiment validation"),
            dbc.Nav(
                evaluation_modules,
                vertical=True,
                pills=True,
            ),
        ]
    ),
    html.Div(
        [
            html.H2("Databases", className="m-2 mb-4"),
            html.Hr(),
            dcc.Store(
                id="database_present_status",
                data={"ncbi_taxonomy": False,
                      "kegg_map": False}),
            html.Div(
                [
                    # html.I(className="bi bi-check-circle-fill me-3 mt-1 fs-6 text-success"),
                    html.I(className="bi bi-x-circle-fill me-3 ms-3 fs-5 text-danger"),
                    html.P("NCBI Taxonomy", className="fs-5")
                ],
                className="d-inline-flex align-items-start w-100",
                id="ncbi_db_presence_check"
            ),
            html.Div(
                [
                    # html.I(className="bi bi-check-circle-fill me-3 mt-1 fs-6 text-success"),
                    html.I(className="bi bi-x-circle-fill me-3 ms-3 fs-5 text-danger"),
                    html.P("GTDB Taxonomy", className="fs-5")
                ],
                className="d-inline-flex align-items-start w-100"\
                    if gc.show_advanced_settings is True else "d-none",
                id="gtdb_db_presence_check"
            ),
            html.Div(
                [
                    html.I(className="bi bi-check-circle-fill me-3 ms-3 fs-5 text-success"),
                    html.P("KEGG KO Map", className="fs-5")
                ],
                className="d-inline-flex align-items-start w-100"\
                    if gc.display_db_search is True else "d-none",
                id="kegg_map_presence_check"
            ),
            html.Div(
                [
                    html.B("Fetch source:", className="mb-0 fs-5"),
                    dbc.Button("Open", id="fetch_db_modal_open"),
                ],
                className="d-flex align-items-center justify-content-between m-2"
            ),
            dcc.Interval(
                id='start_db_check',
                interval=6*10000, # in milliseconds
                n_intervals=0
            ),
            fetch_db_modal
        ],
        # className="border rounded-3 shadow-sm",
        # style={"background-color": "#f9f9f9"}
    )
]
