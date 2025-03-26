import webview

from MetaPepView.server import app
import plotly.io as pio

from constants import GraphConstants

# import layout elements
from MetaPepView.layout.app_layout import app_layout

# import callbacks
from MetaPepView.callbacks.page_select_callbacks import *
from MetaPepView.callbacks.data_callbacks import *
from MetaPepView.callbacks.taxonomy_callbacks import *
from MetaPepView.callbacks.denovo_callbacks import *
from MetaPepView.callbacks.func_annot_callbacks import *
from MetaPepView.callbacks.denovo_callbacks import *
from MetaPepView.callbacks.ms_callbacks import *
from MetaPepView.callbacks.annotation_callbacks import *
from MetaPepView.callbacks.sidebar_callbacks import *


# assign template theme for plotly figures
pio.templates.default = GraphConstants.default_template
# Build the app layout from the components
app.layout = app_layout

if __name__ == '__main__':
    webview.create_window("MetaPepView", app.server, maximized=True)
    webview.start()
    # app.run(host="0.0.0.0", debug=True)