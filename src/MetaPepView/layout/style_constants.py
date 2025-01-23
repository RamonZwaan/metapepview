import plotly.express as px


# Color styles
PLOT_COLOR = "#f6f6f6"
IMPORT_BLOCK_COLOR = "#"
SIDEBAR_COLOR = "#eaeaea"
HEADER_COLOR = "#483d45"
TAB_COLOR = "#dddddd"

IMPORT_FAILED_COLOR = "#ffe5e5"
IMPORT_SUCCESS_COLOR = "#ecffe6"

# Common layout settings
SIDEBAR_WIDTH = "17%"

# define style arguments for the sidebar
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": SIDEBAR_WIDTH,
    "padding": "2rem 1rem",
    "background-color": SIDEBAR_COLOR,
    "overflow-y": "auto",
    "overflow-x": "hidden",
    "display": "flex",
    "flex-direction": "column",
    "justify-content": "space-between"
}

# define style arguments for the header bar
CONTENT_HEADER_STYLE = {
    "position": "fixed",
    "left": SIDEBAR_WIDTH,
    "right": 0,
    "top": 0,
    "height": "4rem",
    "zIndex": 1, 
}

# define style arguments for the content box
CONTENT_STYLE = {
    "margin-top": "1rem",
    "margin-left": SIDEBAR_WIDTH,
    "margin-right": "1rem",
    "padding": "0rem 0rem",
    "zIndex": 5
}


HEADER_BUTTON_STYLE = {
    "color": TAB_COLOR,
    "padding": "0 1rem",
    "cursor": "pointer",
    "display": "inline-block",
    "white-space": "nowrap"
}


# import container style
# regular import container regular colour
qa_import_box_style = {
    "margin": ".5rem", "padding": "1rem 1rem",
    "background-color": PLOT_COLOR,
    "border-radius": "1rem"
    }



# regular import container success colour, if validation passed
success_box_style = {"background-color": IMPORT_SUCCESS_COLOR}

# regular import container failed colour, if validation failed
failed_box_style = {"background-color": IMPORT_FAILED_COLOR}
