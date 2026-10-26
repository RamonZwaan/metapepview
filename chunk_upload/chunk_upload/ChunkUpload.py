# AUTO GENERATED FILE - DO NOT EDIT

import typing  # noqa: F401
from typing_extensions import TypedDict, NotRequired, Literal # noqa: F401
from dash.development.base_component import Component, _explicitize_args

ComponentSingleType = typing.Union[str, int, float, Component, None]
ComponentType = typing.Union[
    ComponentSingleType,
    typing.Sequence[ComponentSingleType],
]

NumberType = typing.Union[
    typing.SupportsFloat, typing.SupportsInt, typing.SupportsComplex
]


class ChunkUpload(Component):
    """A ChunkUpload component.


Keyword arguments:

- id (string; optional)

- appearance (string; default "default")

- chunk_size (number; optional)

- className (string; optional)

- file_infos (list; optional)

- multiple (boolean; default False)

- progress (number; optional)

- uploadText (string; optional)"""
    _children_props: typing.List[str] = []
    _base_nodes = ['children']
    _namespace = 'chunk_upload'
    _type = 'ChunkUpload'


    def __init__(
        self,
        id: typing.Optional[typing.Union[str, dict]] = None,
        file_infos: typing.Optional[typing.Sequence] = None,
        progress: typing.Optional[NumberType] = None,
        chunk_size: typing.Optional[NumberType] = None,
        multiple: typing.Optional[bool] = None,
        uploadText: typing.Optional[str] = None,
        className: typing.Optional[str] = None,
        appearance: typing.Optional[str] = None,
        style: typing.Optional[typing.Any] = None,
        **kwargs
    ):
        self._prop_names = ['id', 'appearance', 'chunk_size', 'className', 'file_infos', 'multiple', 'progress', 'style', 'uploadText']
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'appearance', 'chunk_size', 'className', 'file_infos', 'multiple', 'progress', 'style', 'uploadText']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs and excess named props
        args = {k: _locals[k] for k in _explicit_args}

        super(ChunkUpload, self).__init__(**args)

setattr(ChunkUpload, "__init__", _explicitize_args(ChunkUpload.__init__))
