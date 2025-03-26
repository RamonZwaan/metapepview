"""
Contains functions that create html blocks for sections that are consistent
in format.
"""

from dash import html, dcc
import dash_bootstrap_components as dbc

from typing import Any, Sequence, List, Optional, Tuple, Dict

from backend.io import *
from constants import StyleConstants


def importer_block(
    title: str,
    upload_id: str,
    return_data_id: str,
    valid_state_id: str,
    multiple=False,
    title_id=None) -> Any:
    
    if multiple is True:
        filenames_feedback = html.Div(
            [
            html.P("No file...", style={"padding": "0rem .5rem",
                                        "text-align": "left"})
            ],
            id=return_data_id,
            className="overflow-auto", 
            style={"height": "5rem"}
        )
    else:
        filenames_feedback = html.Div(
            [
            html.P("No file...", style={"padding": "0rem .5rem",
                                        "text-align": "left"})
            ],
            id=return_data_id, 
        )
        
    if title_id is None:
        header_block = html.H4(title)
    else:
        header_block = html.H4(title, id=title_id)
    
    return [
        dbc.Row(
            [
                dbc.Col(
                    header_block,
                    className="align-self-center",
                    style={"width": "20rem"}
                ),
                dbc.Col(
                    dcc.Upload(
                        [
                            dbc.Button("Import", outline=True, color="primary",
                                    className="px-4")
                        ],
                        multiple=multiple,
                        id=upload_id,
                        className="justify-content-between p-auto float-end"
                    ),
                    width=3
                ),
            ],
            className="justify-content-between my-0"
        ),

        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    dbc.Spinner(
                        filenames_feedback
                    )
                )
            ],
            style={'display': "flex", "align-items": "center"}
        ),
        dcc.Store(id=valid_state_id, data=None)
    ]


def annotation_mini_importer_block(
    upload_id: str,
    return_data_id: str,
    valid_state_id: str,
    format_options: List[str] | List[Dict[str, str]] | None=None,
    format_id: str | None=None,
    allow_multiple: bool=False) -> html.Div:
    """Create dash component template for importer blocks at small size.

    Args:
        upload_id (str): Id of dash upload component.
        return_data_id (str): Id for the data name to be returned.
        valid_state_id (str): Component Id to check if block is valid.
        multiple (bool, optional): Allow multiple files. Defaults to False.
        format_options(List[str] | None, optional): Give option values for
            the potential data formats to upload in a dropdown menu.
            Defaults to None.
        format_id (str | None, optional): Dash component id for format dropdown
            menu. Defaults to None.
        
    Returns:
        Any: List[Any]: Dash component block.
    """
    header_block = html.B("Format", className="me-5")

    if format_options is not None:
        initial_val = format_options[0] if isinstance(format_options[0], str) else format_options[0]['value']
        format_block = dbc.Col(
            dcc.Dropdown(
                format_options,
                value=initial_val,
                clearable=False,
                id=format_id,
                style={'height': "30px"}
            ),
            className="w-auto"
            #width=6
        )
    else:
        format_block = dbc.Col(None)
        
    upload_style = {
        'borderWidth': '1px',
        'borderStyle': 'dashed',
        'borderRadius': '5px',
        'textAlign': 'center',
    }
    if allow_multiple == True:
        filenames_feedback = html.Div([
                'Drag and Drop or ',
                html.A("Select Files", style={'cursor': 'pointer', 'fontWeight': 'bold'})
            ],
            id=return_data_id, 
        )
    else:
        filenames_feedback = html.Div([
                'Drag and Drop or ',
                html.A("Select File", style={'cursor': 'pointer', 'fontWeight': 'bold'})
            ],
            id=return_data_id, 
        )
        

    return html.Div([
        html.Div(
            [
                header_block,
                format_block,
            ],
            className="d-flex align-items-center my-3"
        ),
        dcc.Upload(
            id=upload_id,
            children=dbc.Spinner(
                filenames_feedback,
                size="sm"
            ),
            className="my-2 mx-1 py-1",
            style=upload_style,
            multiple=allow_multiple
        ),
        dcc.Store(id=valid_state_id, data=None)
    ],
    )


def qa_importer_block(
    title: str,
    upload_id: str,
    return_data_id: str,
    valid_state_id: str,
    title_id=None,
    format_options: List[str] | List[Dict[str, str]] | None=None,
    format_id: str | None=None) -> List[Any]:
    """Create dash component template for importer blocks at small size.

    Args:
        title (str): Header text of block.
        upload_id (str): Id of dash upload component.
        return_data_id (str): Id for the data name to be returned.
        valid_state_id (str): Component Id to check if block is valid.
        multiple (bool, optional): Allow multiple files. Defaults to False.
        title_id (_type_, optional): Add id to title component. Defaults to None.
        format_options(List[str] | None, optional): Give option values for
            the potential data formats to upload in a dropdown menu.
            Defaults to None.
        format_id (str | None, optional): Dash component id for format dropdown
            menu. Defaults to None.
        
    Returns:
        Any: List[Any]: Dash component block.
    """
    filenames_feedback = html.Div([
            'Drag and Drop or ',
            html.A("Select File", style={'cursor': 'pointer', 'fontWeight': 'bold'})
        ],
        id=return_data_id, 
    )

    if title_id is None:
        header_block = html.H5(title)
    else:
        header_block = html.H5(title, id=title_id)

    if format_options is not None:
        initial_val = format_options[0] if isinstance(format_options[0], str) else format_options[0]['value']
        format_block = dbc.Col(
            dcc.Dropdown(
                format_options,
                value=initial_val,
                clearable=False,
                id=format_id
            ),
            width=4
        )
    else:
        format_block = dbc.Col(None, width=4)
    
    return [
        dbc.Row(
            [
                dbc.Col(
                    header_block,
                    className="align-self-center",
                    style={"width": "20rem"}
                ),
                format_block
            ],
            className="justify-content-between my-0"
        ),
        dcc.Upload(
            id=upload_id,
            children=dbc.Spinner(
                filenames_feedback,
                size="sm"
            ),
            className="my-2 mx-1 py-1",
            style={
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                #'margin': '1rem'
            },
        ),
        dcc.Store(id=valid_state_id, data=None)
    ]


def import_single_file(names: Optional[str],
                       dates: Any,
                       max_name_len: int=17,
                       drag_and_drop: bool = False) -> List[Any] | html.P:
    """Standard function to display name for single file on dashboard.

    Args:
        names (Optional[str]): File name.
        dates (Any): Date of file.
        max_name_len (int, optional): Maximum length of name to display.
        Defaults to 17.
        drag_and_drop (bool, optional): Return element text made for display inside
            drag and drop field. Defaults to False.

    Returns:
        List[Any] | html.P: filename element to display.
    """
    if names is None:
        # mini import blocks provide drag and drop fields
        if drag_and_drop is True:
            return [
                'Drag and Drop or ',
                html.A("Select File", style={'cursor': 'pointer', 'fontWeight': 'bold'})
            ]
        else:
            return ["No file..."]
    # update name in sidebar
    else:
        if len(names) > max_name_len:
            names = names[:max_name_len-3] + '...'
        return [html.P(names, className="fw-bold")]
    

def import_multiple_files(names: Optional[Sequence[str]],
                          dates: Optional[Sequence[str]],
                          max_name_len: int=17,
                          max_rows: int=0) -> List[html.P] | html.P:
    """Standard function to display names for multiple uploaded files
    on dashboard.

    Args:
        names (Optional[Sequence[str]]): List of file names.
        dates (Optional[Sequence[str]]): List of dates from file names
        max_name_len (int, optional): maximal length of name to display.
            Defaults to 17.
        max_rows (int, optional): Maximal number of rows to display.
            Defaults to 0.

    Returns:
        List[html.P] | html.P: List of html p elements to display.
    """
    
    if names is None or dates is None:
        return html.P("No file...")
    
    # limit displayed rows if desired
    if max_rows > 0 and max_rows < len(names):
        names, dates = names[:max_rows], dates[:max_rows]

    # get all names in the output list
    names_constr = []
    for name in names:
        if len(name) > max_name_len:
            names_constr.append(name[:max_name_len-3] + '...')
        else:
            names_constr.append(name)

    # create list of P elements to put in div
    return [html.P(name, className="fw-bold") for name in names_constr]


def validate_multiple_files(contents: List[str] | None,
                            names: List[str] | None,
                            dates: Any,
                            content_validator: Callable[[str, str | None], Tuple[bool, str | None]]) -> Tuple[
                            bool | None,
                            List[html.P] | html.P,
                            List[str] | None,
                            str | None,
                            bool,
                            Dict[str, str]]:
    """Display filename of denovo dataset import.
    """
    msg = None
    box_style = {}
    
    # validate file contents, if invalid format, notify user
    if contents is not None and names is not None:
        # allow import of single archive, extract all files
        archive_format = determine_archive_format(names[0])
        if len(contents) == 1 and archive_format is not None:
            file_data, file_names = archive_to_file_list(contents[0],
                                                         archive_format)
        else:
            file_data, file_names = contents, names
        
        current_name = None
        try:
            for idx, file in enumerate(file_data):
                if file is None:
                    raise ValueError("Non-file-like encountered, likely directory")
                
                # read buffer to keep content validator input consistent
                file = file.read() if isinstance(file, IO) else file

                current_name = file_names[idx]
                current_archive_format = determine_archive_format(current_name)
                success, msg = content_validator(file, current_archive_format)
                if success is False:
                    break
        except:
            success, msg = False, f"failed to read '{current_name}'"

        # if all files pass validation, update display names
        if success is True:
            valid_data = True
            box_style = StyleConstants.success_box_style
        else:
            name_list = import_multiple_files(None, None)
            return (False, name_list, None, msg, True,
                    StyleConstants.failed_box_style)
    else:
        valid_data = None
        file_names= None
        
    name_list = import_multiple_files(file_names,
                                      dates,
                                      max_name_len=30)
    
    return (valid_data, name_list, contents, msg, False, box_style)


def validate_single_file(contents: str | None,
                         name: str | None,
                         dates: Any,
                         content_validator: Callable[[str, str | None], Tuple[bool, str | None]],
                         drag_and_drop: bool = False) -> Tuple[
                            bool | None,
                            List[Any] | html.P,
                            str | None,
                            str | None,
                            bool,
                            Dict[str, str]]:
    """Perform initial validation and file processing to any data
    uploaded to a single file upload component. Returns
    file metadata, success status, potential error message, and
    styling of importer block based on success.

    Args:
        contents (str): File contents,
        name (str): File name
        dates (Any): Last modified date
        content_validator (Callable[[str, str  |  None], Tuple[bool, str  |  None]]): 
            Function that takes content file, as well as compression format as input
            and returns success status and error message in case of fail as output.
        drag_and_drop (bool, optional): Return element text made for display inside
            drag and drop field. Defaults to False.

    Returns:
        Tuple[bool | None, List[Any] | html.P, str | None, str | None, bool, Dict[str, str]]: 
            Validation status as well as html elements: 
            {valid data status, filename block, file contents, error message, success status,
            style element for import box}
    """
    msg, valid_data, box_style, success = None, None, {}, True
    
    if contents is not None and name is not None:
        # validate format if data given
        archive_format = determine_archive_format(name)
        success, msg = content_validator(contents, archive_format)

        # assign box colour based on success of failed
        if success:
            valid_data = True
            box_style = StyleConstants.success_box_style
        else:
            valid_data = False
            box_style = StyleConstants.failed_box_style
            contents = None

    return (valid_data,
            import_single_file(name, dates, max_name_len=30, drag_and_drop=drag_and_drop),
            contents,
            msg,
            not success,
            box_style)


def hidden_graph_with_text(graph_id: str, text_overlay: str, div_encapsulate: bool=False) -> List[Any]:
    """Create block that contains hidden dash graph and a text overlay.
    Used to give information to the user if graph cannot be created yet
    while having a valid id for further callback processing

    Args:
        graph_id (str): Id of graph object.
        text_overlay (cstr): Text to display in graph place
        div_encapsulate (bool, optional): _description_. Defaults to False.
    """
    block_element = [
        html.Div(dcc.Graph(id=graph_id), style={'display': 'None'}),
        html.P(text_overlay)
    ]
    return block_element


def sample_color_table_block(sample_color_map: Dict[str, str]) -> List[Any]:
    
    def row_block(sample: str,
                  color: str):
        return dbc.Row(
            [
                dbc.Col(
                    html.P(sample, 
                           className="ms-5"), 
                    width=8),
                dbc.Col(
                    html.Div(style={'backgroundColor': color, 
                                    'height': '1rem', 
                                    'width': '1rem'}),
                    width=4
                )
            ],
            className="d-flex align-items-center"
        )
    
    div_children = []
    
    for sample, color in sample_color_map.items():
        div_children.append(row_block(sample, color))
    
    return div_children