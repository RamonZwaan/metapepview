from dash import Output, Input, ctx

from metapepview.server import app

from metapepview.layout.annotation_page import import_block
from metapepview.layout.quality_control_page import ms_performance
from metapepview.layout.taxonomy_page import taxonomy_sample_analysis, \
    taxonomy_de_novo_analysis
from metapepview.layout.func_annot_page import functional_annotation_page
from metapepview.layout.data_page import data_visual



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
    
    return data_visual
        
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