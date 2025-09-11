from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from metapepview.constants import StyleConstants


# Data page
data_visual = [
    dbc.Tabs(
        [            
            dbc.Tab(
                [
                    html.Div(
                        id="sample_table",
                        style={"margin": "1rem", 
                               "padding": "1rem 1rem", 
                               "height": "25rem", 
                               "background-color": StyleConstants.plot_color, 
                               "border-radius": "1rem"}
                    ),
                    html.Div(
                        [
                            dbc.Button('Export peptide data', color="secondary", id='start_download_button',
                                n_clicks=0),
                            dcc.Download(id="download_peptides_dataset_button")
                        ],
                        style={'float': 'bottom', 'margin-left': '1rem', 'width': '15rem'} 
                    ),
                ],
                label="Samples"
            ),
            dbc.Tab(
                [
                    html.Div(
                        id="db_search_psm_table",
                        style={"margin": "1rem", 
                               "padding": "1rem 1rem", 
                               "height": "25rem", 
                               "background-color": StyleConstants.plot_color, 
                               "border-radius": "1rem"}
                    ),
                    html.Div(
                        dcc.Dropdown(
                            placeholder='Sample',     
                            id='psm_table_selector'
                        ),
                        style={'float': 'bottom', 'margin-left': '1rem', 'width': '15rem'}, 
                    )
                ],
                label="DB Search PSM"
            ),
            dbc.Tab(
                html.Div(
                    [
                        html.H5("Protein DB properties")    
                    ],
                    id="taxonomy_db_table",
                    style={"margin": "1rem", 
                           "padding": "1rem 1rem", 
                           "height": "25rem", 
                           "background-color": StyleConstants.plot_color, 
                           "border-radius": "1rem"}
                ),
                label="Protein DB"
            ),
            dbc.Tab(
                html.Div(
                    [
                        html.H5("Functional Annotation table")    
                    ],
                    id="func_annot_table",
                    style={"margin": "1rem", 
                           "padding": "1rem 1rem", 
                           "height": "25rem", 
                           "background-color": StyleConstants.plot_color, 
                           "border-radius": "1rem"}
                ),
                label="Functional Annotation"
            ),
        ]
    )
]