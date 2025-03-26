from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
import dash_bootstrap_components as dbc

from constants import StyleConstants


# Elements for the header bar
content_header = html.Div(
    [
    dbc.Navbar(
        dbc.Container(
            [
                dbc.NavbarBrand("Taxonomic and functional analysis of microbial communities", class_name="ms-2"),
                dbc.Row(# tab menu items
                    [
                        dbc.Col(html.P("Processing", 
                                       id='processing_button', 
                                       style=StyleConstants.header_button_style)),
                        dbc.Col(html.P("Data", 
                                       id='data_button', 
                                       style=StyleConstants.header_button_style)),
                        dbc.Col(html.P("Method Performance", 
                                       id='validation_button', 
                                       style=StyleConstants.header_button_style)),
                        dbc.Col(html.P("Taxonomy", 
                                       id='taxonomy_button', 
                                       style=StyleConstants.header_button_style)),
                        dbc.Col(html.P("Functional", 
                                       id='functional_button', 
                                       style=StyleConstants.header_button_style)),
                    ],
                    style={"float": "right", "margin-right": "1rem", "margin-top": "0.5rem"} 
                )
            ],
            fluid=True
        ),
        color=StyleConstants.header_color,
        dark=True,
    ),
    ],
    style=StyleConstants.content_header_style
)