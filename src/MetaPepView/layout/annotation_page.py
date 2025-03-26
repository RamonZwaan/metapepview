from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from MetaPepView.html_templates import *



sample_options_block = [
    html.Div(
        [
            html.H3("New sample"),
            html.Div(
                [
                    html.H4("Name: ", className="align-self-center me-4"),
                    dbc.Input(id="sample_name_import", type="text"),
                ],
                className="d-flex justify-content-between align-items-center",
                style={"width": "30rem"}
            ),
            html.Div(
                [
                    html.H4("merge db search files", className="align-self-center me-4"),
                    dbc.Switch(id="merge_psm_switch", className="align-self-center", value=True),
                ],
                className="d-flex justify-content-start"
            ),
            html.Div(
                [
                    dbc.Spinner(html.Div(id="annotation_loading_spot", style={"width": "5rem"})),
                    html.Div(
                        dbc.Button("Add Data", id="start_annotation_button",
                                color="primary",
                                className="px-3"),
                        id="start_annotation_button_wrapper",
                        className="me-1"
                    ),
                    html.Div(id="annotation_hint")
                ],
                className="d-flex align-items-center",
            ),
        ],
        className="d-flex justify-content-between"
    ),
    # proteomics_data,
    # de_novo_data,
    # taxonomic_annotation,
    # functional_annotation,
    
    dcc.Store(id="current_taxonomy_db_loc"),
    # dcc.Store(id="valid_import_data"),
]

peptide_data_block = html.Div(
    [
        html.Div(
            [
                html.H3("Import/Export experiment"),
                html.Div(
                    [
                        html.H4("Experiment name: ", className="align-self-center"),
                        dbc.Input(id="experiment_name_field",
                                  type="text",
                                  style={'width': '25rem'},
                                  className="fw-bold"),
                    ],
                    className="d-flex justify-content-between align-items-center",
                    style={"width": "40rem"}
                ),
                html.Div(
                    [
                        dbc.Button('Clear data', id='clear_peptides_data', className="px-3 me-2"),
                        dbc.Button('Export csv', id='export_peptides_csv', className="px-3 me-2"),
                        dbc.Button('Export json', id='export_peptides_json', className="px-3 me-2"),
                        dcc.Upload(
                            [
                                dbc.Button('Import json', id='import_peptides',
                                           outline=True,
                                           color="primary",
                                           className="px-3 me-1 my-0",
                                           style={"height": "100%"}),
                            ],
                            className="p-0 my-0",
                            multiple=False,
                            id="peptides_obj_upload",
                            style={"height": "100%"}
                        ),
                        dcc.Download(id="download_peptides_csv"),
                        dcc.Download(id="download_peptides_json"),
                    ],
                    className="d-flex"
                ),
            ],
            style={'float': 'bottom'},
            className="d-flex justify-content-between my-3 mx-3"
        ),
        html.Hr(className="my-0 py-0"),
        html.Div(
            [
                html.Div(
                    [
                        html.H4("DB search format", className="text-secondary"),
                        html.H6("-", id="peptides_db_search_format", className="ps-4 mb-1")
                    ],
                    className="w-25 py-2 px-4"
                ),
                html.Div(className="vr my-0"),
                html.Div(
                    [
                        html.H4("De novo format", className="text-secondary"),
                        html.H6("-", id="peptides_de_novo_format", className="ps-4 mb-1")
                    ],
                    className="w-25 py-2 px-4"
                ),
                html.Div(className="vr"),
                html.Div(
                    [
                        html.H4("Taxonomy db format", className="text-secondary"),
                        html.H6("-", id="peptides_taxonomy_db_format", className="ps-4 mb-1")
                    ],
                    className="w-25 py-2 px-4"
                ),
                html.Div(className="vr"),
                html.Div(
                    [
                        html.H4("Function db format", className="text-secondary"),
                        html.H6("-", id="peptides_function_db_format", className="ps-4 mb-1")
                    ],
                    className="w-25 py-2 px-4"
                ),
            ],
            style={'float': 'bottom'},
            className="d-flex justify-content-between my-0 mx-0"
        ),
    ],
    className="border border-1 bg-light rounded-3 shadow-sm",
)



db_search_options_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("DB Search filter settings")),
        dbc.ModalBody(
            [
                html.Div(
                    [   
                        html.B("Confidence threshold:", className="me-3 align-self-center"),
                        dbc.Input(value=30, min=0, id="db_search_psm_score_threshold", size="sm", type="number", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-between mb-3"
                ),
                html.Div(
                    [   
                        dbc.Checkbox(label= html.B("filter cRAP",
                                                   id="db_search_filter_crap_text",
                                                   className="ms-2 me-3 align-top text-decoration-underline"),
                                     id="db_search_filter_crap",
                                     value=False),
                        dbc.Popover("""
                                    Ignore peptides that occur within the cRAP dataset, which
                                    is a dataset that contains sequences common encountered
                                    in the environment during sample preparation.
                                    """,
                            id="db_search_filter_crap_info",
                            target="db_search_filter_crap",
                            trigger="hover",
                            placement='top',
                            className="px-3 py-1"
                        )
                    ],
                    className="d-flex mt-4 justify-content-start align-items-center"
                ),
                
            ],
            className="vstack p-5"
        )
    ],
    id="db_search_modal",
    scrollable=True,
    is_open=False,
)


de_novo_options_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("De novo filter settings")),
        dbc.ModalBody(
            [
                html.Div(
                    [   
                        html.B("Confidence threshold:", className="me-3 align-self-center"),
                        dbc.Input(value=30, min=0, id="de_novo_score_threshold", size="sm", type="number", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-between mb-3"
                ),
                html.Div(
                    [   
                        dbc.Checkbox(label= html.B("filter cRAP",
                                                   id="De_novo_filter_crap_text",
                                                   className="ms-2 me-3 align-top text-decoration-underline"),
                                     id="de_novo_filter_crap",
                                     value=False),
                        dbc.Popover("""
                                    Ignore peptides that occur within the cRAP dataset, which
                                    is a dataset that contains sequences common encountered
                                    in the environment during sample preparation.
                                    """,
                            id="de_novo_filter_crap_info",
                            target="de_novo_filter_crap",
                            trigger="hover",
                            placement='top',
                            className="px-3 py-1"
                        )
                    ],
                    className="d-flex mt-4 justify-content-start align-items-center"
                ),
                
            ],
            className="vstack p-5"
        )
    ],
    id="de_novo_modal",
    scrollable=True,
    is_open=False,
)


taxonomy_map_options_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Taxonomy map filter settings")),
        dbc.ModalBody(
            [
                html.Div(
                    [   
                        html.B("Delimiter:", className="me-3 align-self-center"),
                        dbc.Input(value=r'\t', id="acc_tax_map_delim", size="sm", type="text", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-between mb-3",
                    style={"width": "20rem"}
                ),
                html.Div(
                    [   
                        html.Div(
                            [
                                html.B("Accession column index:", className="me-3 align-self-center"),
                                dbc.Input(value=0, min=0, id="acc_tax_map_acc_idx", size="sm", type="number", style={"width": "4rem"}),
                            ],
                            className="d-flex justify-content-between",
                            style={"width": "20rem"}
                        ),
                        html.Div(
                            [
                                html.B("Pattern match (regex): ", className="ms-5 me-3 align-self-center"),
                                dbc.Input(id="acc_tax_map_acc_pattern", size="sm", type="text", style={"width": "10rem"}),
                            ],
                            className="d-flex justify-content-between",
                            style={"width": "25rem"}
                        )
                    ],
                    className="d-flex justify-content-between mb-3",
                    style={"width": "45rem"}
                ),
                html.Div(
                    [   
                        html.B("Taxonomy column index:", className="me-3 align-self-center"),
                        dbc.Input(value=1, min=0, id="acc_tax_tax_idx", size="sm", type="number", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-between",
                    style={"width": "20rem"}
                ),
                html.Div(
                    [   
                        # html.B("Combine multiple annotations",
                        #        id="func_annot_combine_text",
                        #        className="me-3 align-top"),
                        dbc.Checkbox(label= html.B("global annotation of peptides",
                                        id="global_taxonomy_annotation_text",
                                        className="ms-2 me-3 align-top text-decoration-underline"
                                        ),
                                        id="global_taxonomy_annotation_checkbox",
                                        value=False),
                        dbc.Popover("""
                                    As an additional step to annotation of protein id's, match
                                    peptide sequences directly to the UniprotKB protein
                                    database through the Unipept api.
                            """,
                            id="global_taxonomy_annotation_info",
                            target="global_taxonomy_annotation_text",
                            trigger="hover",
                            placement='top',
                            className="px-3 py-1"
                        )
                    ],
                    className="d-flex mt-4 justify-content-start align-items-center"
                ),
            ],
            className="vstack p-5"
        )
    ],
    id="taxonomy_map_modal",
    size="lg",
    scrollable=True,
    is_open=False,
)


function_map_options_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Function map filter settings")),
        dbc.ModalBody(
            [
                html.Div(
                    [   
                        # html.B("Combine multiple annotations",
                        #        id="func_annot_combine_text",
                        #        className="me-3 align-top"),
                        dbc.Checkbox(label= html.B("Combine multiple annotations",
                                        id="func_annot_combine_text",
                                        className="ms-2 me-3 align-top text-decoration-underline"
                                        ),
                                        id="func_annot_combine",
                                        value=False),
                        dbc.Popover("""
                            When a peptide matches to multiple accessions, they might receive multiple functional annotations.
                            this option will combine values of all accessions within single table cells.
                            This way, no functional information is discarded in the data. However,
                            these combined values may distort functional annotation visualization within the dashboard.
                            """,
                            id="func_annot_combine_info",
                            target="func_annot_combine_text",
                            trigger="hover",
                            placement='top',
                            className="px-3 py-1"
                        )
                    ],
                    className="d-flex justify-content-start align-items-center"
                )
            ],
            className="vstack p-5"
        )
    ],
    id="function_map_modal",
    scrollable=True,
    is_open=False,
)


db_search_import_block = [
    html.Div(
        [
            html.Div(
                [
                    html.H4("DB search"),
                    dbc.Button("Options", id="db_search_modal_open", size="sm", color="secondary", outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "db_search_psm_upload", "db_search_psm_import_txt", "db_search_psm_valid",
                format_options=GlobalConstants.db_search_dropdown_options,
                format_id="db_search_psm_format",
                allow_multiple=True
            ),
        ],
        className="mx-4 mb-3"
    ),
    html.Hr(className="w-100 my-0"),
    html.Div(
        ["No file..."],
        id="db_search_psm_name",
        className="px-4 pt-3 overflow-auto",
        style={"height": "10rem"}
    ),
]


de_novo_import_block = [
    html.Div(
        [
            html.Div(
                [
                    html.H4("De novo"),
                    dbc.Button("Options", id="de_novo_modal_open", size="sm", color="secondary", outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "denovo_upload", "de_novo_import_txt", "de_novo_valid",
                format_options=GlobalConstants.de_novo_dropdown_options,
                format_id="de_novo_format",
                allow_multiple=True
            )
        ],
        className="mx-4 mb-3"
    ),
    html.Hr(className="w-100 my-0"),
    html.Div(
        ["No file..."],
        id="denovo_name",
        className="px-4 pt-3 overflow-auto",
        style={"height": "10rem"}
    ),
]



taxonomy_map_import_block = [
    html.Div(
        [
            html.Div(
                [
                    html.H4("Taxonomy annotations"),
                    dbc.Button("Options", id="taxonomy_map_modal_open", size="sm", color="secondary", outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "taxonomy_db_upload", "tax_map_import_txt", "taxonomy_db_valid",
                format_options=[{'label': 'NCBI', 'value': 'NCBI'},
                                {'label': 'GTDB', 'value': 'GTDB'}],
                format_id="taxonomy_db_format"
            )
        ],
        className="mx-4 mb-3"
    ),
    html.Hr(className="w-100 my-0"),
    html.Div(
        ["No file..."],
        id="taxonomy_db_name",
        className="px-4 pt-3 overflow-auto",
        style={"height": "10rem"}
    ),
]


function_map_import_block = [
    html.Div(
        [
            html.Div(
                [
                    html.H4("Functional annotations"),
                    dbc.Button("Options", id="function_map_modal_open", size="sm", color="secondary", outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "func_annot_db_upload", "func_map_import_txt", "func_annot_db_valid",
                format_options=[{'label': 'EggNOG', 'value': 'EggNOG'},
                                {'label': 'gKOALA', 'value': 'gKOALA'}],
                format_id="func_annot_db_format"
            )
        ],
        className="mx-4 mb-3"
    ),
    html.Hr(className="w-100 my-0"),
    html.Div(
        ["No file..."],
        id="func_annot_name",
        className="px-4 pt-3 overflow-auto",
        style={"height": "10rem"}
    ),
]


data_import_container = [
    html.Div(
        db_search_import_block,
        id="db_search_import_box",
        className="py-3 w-25 overflow-visible",
    ),
    html.Div(className="vr"),
    html.Div(
        de_novo_import_block,
        id="de_novo_import_box",
        className="py-3 w-25 overflow-visible",
    ),
    html.Div(className="vr"),
    html.Div(
        taxonomy_map_import_block,
        id="taxonomy_db_import_box",
        className="py-3 w-25 overflow-visible",
    ),
    html.Div(className="vr"),
    html.Div(
        function_map_import_block,
        id="functional_db_import_box",
        className="py-3 w-25 overflow-visible",
    )
]


import_block = html.Div(
    [
        peptide_data_block,

        dbc.Alert(
            id="taxonomy_db_format_alert",
            dismissable=True,
            is_open=False,
            color="danger"
        ),
        dbc.Alert(
            id="functional_db_format_alert",
            dismissable=True,
            is_open=False,
            color="danger"
        ),
        dbc.Alert(
            id="de_novo_format_alert",
            dismissable=True,
            is_open=False,
            color="danger"
        ),
        dbc.Alert(
            id="db_search_format_alert",
            dismissable=True,
            is_open=False,
            color="danger"
        ),
        db_search_options_modal,
        de_novo_options_modal,
        taxonomy_map_options_modal,
        function_map_options_modal,
        html.Div(
            [
                html.Div(
                    sample_options_block,
                    id="peptide_import_container",
                    className="p-3",
                    #style={"height": "31rem"}
                ),
                html.Hr(className="m-0"),                
                html.Div(
                    data_import_container,
                    id="data_import_container",
                    className="d-flex",
                ),
                # html.Hr(className="m-0"),

            ],
            className="border border-2 my-3 bg-light rounded-3 shadow-sm",
        ),
        html.Div(
            [
                html.H3("Sample table", className="ps-3 pt-3"),
                html.Div(
                    [
                        dash_table.DataTable(
                            columns=[{'id': c, 'name': c} for c in GlobalConstants.experiment_sample_table_cols],
                            id="experiment_sample_table",
                            style_data={'table-layout': 'fixed', 'backgroundColor': 'rgba(0, 0, 0, 0)'},
                            style_header={'backgroundColor': 'rgb(216, 216, 235)',
                                          'color': 'black',
                                          'fontWeight': 'bold'},
                            style_cell={'textAlign': 'left', 'font-family': 'Arial'},
                            style_cell_conditional=[
                                {'if': {'column_id': 'DB Search Imported'},
                                'width': '12%'},
                                {'if': {'column_id': 'De Novo Imported'},
                                'width': '12%'}],
                            style_as_list_view=True,
                            row_deletable=True,
                            page_size=10)
                    ],
                    id="sample_table",
                    className="",
                    style={"margin": "0rem 0rem", "padding": "0.5rem"}
                ),
            ],
            className="border border-2 my-3 bg-light rounded-3 shadow-sm",
        ),
    ],
    id="import_block",
    style={"margin": "1rem"}
)
