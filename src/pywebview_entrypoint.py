import sys
import webview
from waitress import serve
from threading import Thread

from metapepview.server import app
import plotly.io as pio

from metapepview.constants import GraphConstants

# import layout elements
from metapepview.layout.app_layout import app_layout

# import callbacks
from metapepview.callbacks.page_select_callbacks import *
from metapepview.callbacks.data_callbacks import *
from metapepview.callbacks.taxonomy_callbacks import *
from metapepview.callbacks.denovo_callbacks import *
from metapepview.callbacks.func_annot_callbacks import *
from metapepview.callbacks.denovo_callbacks import *
from metapepview.callbacks.quality_control_callbacks import *
from metapepview.callbacks.annotation_callbacks import *
from metapepview.callbacks.sidebar_callbacks import *


def mpv_server():
    serve(app.server, host="127.0.0.1", port=8050, threads=64)#type:ignore

def main():    
    # set downcasting behavior to manage FutureWarning in `replace` function
    pd.set_option('future.no_silent_downcasting', True)

    # assign template theme for plotly figures
    pio.templates.default = GraphConstants.default_template
    # Build the app layout from the components
    app.layout = app_layout
    Thread(target=mpv_server,
           daemon=True).start()
    
    webview.create_window("MetaPepView", 
                          "http://127.0.0.1:8050", 
                          maximized=True) #type:ignore
    webview.start()




if __name__ == '__main__':
    sys.exit(main())