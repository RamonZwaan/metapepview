from dash import html, dcc

# import layout elements
from constants import StyleConstants
from MetaPepView.layout.sidebar import new_sidebar

from backend import import_ref_statistics


app_layout = html.Div([html.Div(new_sidebar, style=StyleConstants.sidebar_style, id="sidebar"),
                       # content_header,
                       html.Div(
                           style=StyleConstants.content_style,
                           id="content_div",
                       ),
                       # parameter values stores
                       
                       dcc.Store(id="experiment_name"),
                       
                       # dataset for community analysis
                       dcc.Store(id="peptides"),
                       dcc.Store(id="peptides_metadata"),
                       
                       # dataset for spectral analysis
                       dcc.Store(id="mzml_data"),
                       dcc.Store(id="mzml_peaks_data"),
                       dcc.Store(id="mzml_metadata"),
                       dcc.Store(id="ref_statistics", 
                                 data=import_ref_statistics()),

                       # kegg ko map data
                       dcc.Store(id="kegg_ko_map_data"),
                       dcc.Store(id="kegg_db_class_data"),
                       ],
                       className="overflow-hidden",
                       style={"fontSize": "1rem"})
