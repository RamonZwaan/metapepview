from dash import html, dcc, Output, Input, ctx
from server import app
import plotly.io as pio

from constants import *

# import layout elements
from layout.style_constants import *
from layout.sidebar import new_sidebar
from layout.validation_page import ms_performance
from layout.taxonomy_page import taxonomy_sample_analysis, taxonomy_de_novo_analysis
from layout.func_annot_page import *
from layout.data_page import data_visual
from layout.annotation_page import import_block

# import callbacks
from callbacks.data_callbacks import *
from callbacks.taxonomy_callbacks import *
from callbacks.func_annot_callbacks import *
from callbacks.denovo_callbacks import *
from callbacks.ms_callbacks import *
from callbacks.annotation_callbacks import *

from backend import *


################################################################################
# Layout of Dashboard
################################################################################

# assign template theme for plotly figures
pio.templates.default = GraphConstants.default_template


# Build the app layout from the components
app.layout = html.Div([html.Div(new_sidebar, style=SIDEBAR_STYLE, id="sidebar"),
                       # content_header,
                       html.Div(
                           style=CONTENT_STYLE,
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



################################################################################
# Callback Functions
################################################################################


page_options = {
    'sidebar_data_button': [import_block],
    'sidebar_validation_button': [ms_performance],
    'sidebar_taxonomy_button': [taxonomy_sample_analysis],
    'sidebar_taxonomy_de_novo_button': [taxonomy_de_novo_analysis],
    'sidebar_functional_button': [functional_annotation_page],
}

@app.callback(
    [Output(component_id='content_div', component_property='children')] +\
    [Input(component_id=i, component_property='active') for i in page_options.keys()],

)
def update_tab(*args):
    page_blocks = list(page_options.values())
    for i, btn in enumerate(args):
        if btn is True:
            return page_blocks[i]
    
    return data_visual#, sidebar
        
callback_elems = [Output(component_id=i, component_property='active') for i in page_options.keys()] +\
                 [Input(component_id=i, component_property='n_clicks') for i in page_options.keys()]

@app.callback(callback_elems)
def update_sidebar_active(*args):
    return_list = [False] * len(page_options)
    
    if ctx.triggered_id is None:
        true_id = 0
    else:       
        button_ids = list(page_options.keys())
        true_id = button_ids.index(ctx.triggered_id)
    
    return_list[true_id] = True
    return return_list


if __name__ == '__main__':
    app.run(host="0.0.0.0", debug=True)
    