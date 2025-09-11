from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from metapepview.html_templates import *
from metapepview.constants import GlobalConstants as gc



sample_options_block = [
    html.Div(
        [
            html.H3("Add sample"),
            html.Div(
                [
                    html.H5("Sample name:", className="align-self-center me-4"),
                    dbc.Input(id="sample_name_import", 
                              type="text",
                              style={"width": "20rem"}),
                ],
                className="d-flex justify-content-between align-items-center",
                style={"width": "30rem"}
            ),
            html.Div(
                [
                    html.H5("Merge DB search files",
                        className="align-self-center me-4"),
                    dbc.Switch(id="merge_psm_switch", className="align-self-center", value=True),
                ],
                className="d-flex justify-content-start" if gc.show_advanced_settings is True else "d-none"
            ),
            html.Div(
                [
                    html.H5("Annotate to Unipept",
                            id="global_taxonomy_annotation_text",
                            className="ms-2 me-3 align-top",
                            style={'text-decoration-line': "underline", 
                                   'text-decoration-style': "dotted"}
                    ),
                    dbc.Switch(id="global_taxonomy_annotation_checkbox", 
                               className="align-self-center", 
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
                        className="p-2"
                    )
                ],
                id="global_taxonomy_annotation_container",
                className="d-flex justify-content-start align-items-center"\
                    if gc.show_advanced_settings is True else "d-none",
            ),
            html.Div(
                [
                    dbc.Spinner(html.Div(id="annotation_loading_spot", style={"width": "5rem"})),
                    html.Div(
                        dbc.Button("Add sample", id="start_annotation_button",
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
                html.H3("Project"),
                html.Div(
                    [
                        html.H4("Name: ", className="align-self-center me-3"),
                        dbc.Input(id="experiment_name_field",
                                  type="text",
                                  style={'width': '25rem'},
                                  className="fw-bold"),
                    ],
                    className="d-flex justify-content-center align-items-center",
                    # style={"width": "40rem"}
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
            configure_metadata_format_container(),
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
                html.Div(accession_pattern_options(
                    "db_search_acc_selection_text",
                    "db_search_accession_parser_items",
                    "db_search_accession_popover",
                    "db_search_accession_pattern",
                    "db_search_accession_pattern_container",
                    """
                    Specify what part of the accession (protein id) should be used
                    as identifier to match accession id to DB taxonomy/function data.
                    """
                    )
                ),

                html.Div(
                    [
                        html.B("Confidence threshold:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.Input(value=30, min=0, id="db_search_psm_score_threshold", size="sm", type="number", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-start mt-3 mb-3"
                ),
                html.Div(
                    [
                        dbc.Checkbox(label= html.B("filter cRAP",
                                                   id="db_search_filter_crap_text",
                                                   className="ms-2 me-3 align-top",
                                                   style={"text-decoration-line": "underline", 
                                                          "text-decoration-style": "dotted"}),
                                     id="db_search_filter_crap",
                                     value=True),
                        dbc.Popover("""
                                    Ignore peptides that occur within the cRAP dataset, which
                                    is a dataset that contains sequences common encountered
                                    in the environment during sample preparation.
                                    """,
                            id="db_search_filter_crap_info",
                            target="db_search_filter_crap",
                            trigger="hover",
                            placement='top',
                            className="p-2"
                        )
                    ],
                    className="d-flex mt-4 justify-content-start align-items-center"
                ),

            ],
            className="vstack p-5"
        )
    ],
    size="lg",
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
                        html.B("Confidence threshold:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.Input(value=30, min=0, id="de_novo_score_threshold", size="sm", type="number", style={"width": "4rem"}),
                    ],
                    className="d-flex justify-content-start mb-3"
                ),
                html.Div(
                    [
                        dbc.Checkbox(label= html.B("filter cRAP",
                                                   id="De_novo_filter_crap_text",
                                                   className="ms-2 me-3 align-top",
                                                   style={"text-decoration-line": "underline", 
                                                          "text-decoration-style": "dotted"}),
                                     id="de_novo_filter_crap",
                                     value=True),
                        dbc.Popover("""
                                    Ignore peptides that occur within the cRAP dataset, which
                                    is a dataset that contains sequences common encountered
                                    in the environment during sample preparation.
                                    """,
                            id="de_novo_filter_crap_info",
                            target="de_novo_filter_crap",
                            trigger="hover",
                            placement='top',
                            className="p-2"
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
                        html.B("Delimiter:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.Input(value=r'\t', id="acc_tax_map_delim", size="sm", type="text", style={"width": "4rem"}),
                    ],
                    id="tax_acc_map_delimiter_container",
                    # className CONFIGURED IN CALLBACK
                    #style={"width": "20rem"}
                ),
                html.Div(
                    [
                        html.B("Accession type:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.RadioItems(options=[
                            {"label": "Protein id", "value": "Accession"}, 
                            {"label": "Peptide sequence", "value":"Sequence"}
                            ],
                            value="Accession",
                            id="tax_accession_format_radio",
                            inline=True,
                        )
                    ],
                    id="tax_acc_map_acc_type_container",
                    # className CONFIGURED IN CALLBACK
                ),
                html.Div(accession_pattern_options(
                            "accession_parser_text",
                            "accession_parser_checkbox",
                            "accession_parser_info",
                            "acc_tax_map_acc_pattern", 
                            "acc_tax_map_acc_pattern_container",
                            """
                                Specify what part of the accession should be used
                                as identifier to match accession id to DB search data.
                            """
                )),
                html.Div(
                    [
                        html.B("Accession column index:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.Input(value=0, 
                                  min=0, 
                                  id="acc_tax_map_acc_idx", 
                                  size="sm", 
                                  type="number", 
                                  style={"width": "4rem"}),
                        dbc.FormText("index starts at 0", className="ms-3 fst-italic align-self-center")
                    ],
                    id="tax_acc_map_acc_idx_container",
                    # className CONFIGURED IN CALLBACK
                    #style={"width": "20rem"}
                ),
                html.Div(
                    [
                        html.B("Taxonomy column index:", 
                               className="align-self-center",
                               style={"width": "18rem"}),
                        dbc.Input(value=1, 
                                  min=0, 
                                  id="acc_tax_tax_idx", 
                                  size="sm", 
                                  type="number", 
                                  style={"width": "4rem"}),
                        dbc.FormText("index starts at 0", className="ms-3 fst-italic align-self-center")
                    ],
                    id="tax_acc_map_tax_idx_container",
                    # className CONFIGURED IN CALLBACK
                   # style={"width": "20rem"}
                ),
                html.Div(
                    [
                        html.B("Taxonomy element format:",
                               id="taxonomy_id_format_text",
                               className="align-self-center align-top",
                               style={"text-decoration-line": "underline", 
                                      "text-decoration-style": "dotted",
                                      "width": "18rem"}),
                        dbc.RadioItems(options=[
                                       "taxonomy id",
                                       "taxonomy name"
                                       ],
                                       id="taxonomy_id_format_checkbox",
                                       inline=True,
                                       value="taxonomy id"),
                        dbc.Popover("""
                                    Specify if taxonomy column provides NCBI or GTDB
                                    taxonomy id's or taxonomy names.
                            """,
                                    # In NCBI, id's
                                    # are integer values that represent a taxonomy element.
                                    # For example, the Genus Escherichia has the id: '562'.
                                    # If the taxonomy column contains names, then these will be matched against
                                    # scientific names stored in the NCBI taxonomy database.
                                    # These are guaranteed unique for each taxonomy element.\n
                                    # For GTDB, taxonomy id's are stored as taxonomy name with a rank suffix.
                                    # For example, the genus Escherichia has id: 'g__Escherichia'.
                                    # If the taxonomy column contains names, the values will be matched against gtdb id's without the rank suffix: 'Escherichia'.
                                    # At the strain level, GTDB only stores NCBI genome names. These are valid taxonomy identifiers for both names and id's.
                            # """,
                            id="taxonomy_id_format_info",
                            target="taxonomy_id_format_text",
                            trigger="hover",
                            placement='top',
                            className="p-2" ,
                            style={"width": "100rem"}
                        )
                    ],
                    id="tax_acc_map_tax_type_container",
                    hidden=not gc.show_advanced_settings,
                    # className CONFIGURED IN CALLBACK
                ),
                html.Div(
                    [
                        dbc.Checkbox(label= html.B("To NCBI taxonomy id (genome id's only)",
                                        id="gtdb_genome_to_ncbi_text",
                                        className="ms-2 me-3 align-top",
                                        style={"text-decoration-line": "underline", 
                                               "text-decoration-style": "dotted"}
                                        ),
                                        id="gtdb_genome_to_ncbi_checkbox",
                                        value=False),
                        dbc.Popover("""
                            Convert genome id from GTDB database to NCBI taxonomy id.
                            After conversion, proteins will be mapped to the NCBI taxonomy DB.
                            NOTE: Only Genome id's present in the GTB database can be mapped to NCBI taxonomy format.
                            Any protein mapped to a GTDB taxonomy id will be discarded!""",
                            id="gtdb_genome_to_ncbi_info",
                            target="gtdb_genome_to_ncbi_text",
                            trigger="hover",
                            placement='top',
                            className="p-2"
                        )
                    ],
                    id="gtdb_genome_to_ncbi_container",
                    hidden=True,
                    # className CONFIGURED IN CALLBACK
                )
            ],
            className="vstack p-5"
        )
    ],
    id="taxonomy_map_modal",
    size="lg" if gc.show_advanced_settings is True else False,
    scrollable=True,
    is_open=False,
)


function_map_options_modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Function map filter settings")),
        dbc.ModalBody(
            [
                html.Div(accession_pattern_options(
                    "func_accession_parser_text",
                    "func_accession_parser_checkbox",
                    "func_accession_parser_info",
                    "func_annot_acc_pattern", 
                    "func_annot_acc_pattern_container",
                    """
                        Specify what part of the protein id string should be used
                        as identifier to match id to DB search data.
                    """
                )),
                html.Div(
                    [
                        # html.B("Combine multiple annotations",
                        #        id="func_annot_combine_text",
                        #        className="me-3 align-top"),
                        dbc.Checkbox(label= html.B("Combine multiple annotations",
                                        id="func_annot_combine_text",
                                        className="ms-2 me-3 align-top",
                                        style={"text-decoration-line": "underline", 
                                               "text-decoration-style": "dotted"}
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
                            className="p-2"
                        )
                    ],
                    className="d-flex justify-content-start align-items-center mt-3"
                )
            ],
            className="vstack p-5"
        )
    ],
    id="function_map_modal",
    scrollable=True,
    is_open=False,
    size="lg"
)


db_search_import_block = [
    html.Div(
        [
            html.Div(
                [
                    html.H4("DB search"),
                    dbc.Button("Options",
                               id="db_search_modal_open",
                               size="sm",
                               className="" if gc.show_advanced_settings is True else "d-none",
                               color="secondary",
                               outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "db_search_psm_upload", "db_search_psm_import_txt", "db_search_psm_valid",
                format_options=gc.db_search_dropdown_options,
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
                    dbc.Button("Options",
                               id="de_novo_modal_open",
                               size="sm",
                               className="" if gc.show_advanced_settings is True else "d-none",
                               color="secondary",
                               outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "denovo_upload", "de_novo_import_txt", "de_novo_valid",
                format_options=gc.de_novo_dropdown_options,
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
                    html.H4("Taxonomy annotation"),
                    dbc.Button("Options", id="taxonomy_map_modal_open", size="sm", color="secondary", outline=True),
                ],
                className="d-flex justify-content-between mb-4 align-items-center"
            ),
            annotation_mini_importer_block(
                "taxonomy_db_upload", "tax_map_import_txt", "taxonomy_db_valid",
                format_options=[{'label': 'gKOALA', 'value': 'gKOALA'},
                                {'label': 'NCBI', 'value': 'NCBI'}] +
                    ([{'label': 'GTDB', 'value': 'GTDB'}] if gc.show_advanced_settings is True\
                    else []),

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
                    html.H4("Functional annotation"),
                    dbc.Button("Options",
                        id="function_map_modal_open",
                        className="" if gc.show_advanced_settings is True else "d-none",
                        size="sm",
                        color="secondary",
                        outline=True),
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


data_import_container = configure_import_container(
    db_search_block=db_search_import_block,
    de_novo_block=de_novo_import_block,
    taxonomy_block=taxonomy_map_import_block,
    function_block=function_map_import_block
)


import_block = html.Div(
    [
        peptide_data_block,
        dbc.Alert(
            id="peptides_import_format_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Alert(
            id="taxonomy_db_format_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Alert(
            id="functional_db_format_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Alert(
            id="de_novo_format_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Alert(
            id="db_search_format_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
        ),
        dbc.Alert(
            id="annotation_error_alert",
            dismissable=True,
            is_open=False,
            color="danger",
            className="mt-3"
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
                html.H3("Samples in project", className="ps-3 pt-3"),
                html.Div(
                    [
                        dash_table.DataTable(
                            id="experiment_sample_table",
                            style_data={'table-layout': 'fixed', 'backgroundColor': 'rgba(0, 0, 0, 0)'},
                            style_header={'backgroundColor': 'rgb(216, 216, 235)',
                                          'color': 'black',
                                          'fontWeight': 'bold'},
                            style_cell={'textAlign': 'left', 'font-family': 'Arial'},
                            style_cell_conditional=[
                                {'if': {'column_id': 'DB Search Imported'},
                                'width': '15%'},
                                {'if': {'column_id': 'De Novo Imported'},
                                'width': '15%'}],
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
