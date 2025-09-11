import sys
from textwrap import dedent
from waitress import serve

from metapepview.server import app

import plotly.io as pio
import pandas as pd

from metapepview.constants import GraphConstants, GlobalConstants

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


def main():
    """Entrypoint to start MetaPepView dashboard.
    """
    # set downcasting behavior to manage FutureWarning in `replace` function
    pd.set_option('future.no_silent_downcasting', True)

    # assign template theme for plotly figures
    pio.templates.default = GraphConstants.default_template
    
    # Build the app layout from the components
    app.layout = app_layout
    
    # dash development server
    # app.run(host="0.0.0.0", debug=True)

    # production server
    print(
        dedent(f"""
        Starting MetaPepView server...
        Access dashboard in: http://127.0.0.1:{GlobalConstants.port}
        """)
    )
    serve(app.server, host="0.0.0.0", port=GlobalConstants.port, threads=64) #type:ignore


if __name__ == '__main__':
    sys.exit(main())
    