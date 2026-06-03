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

- chunk_size (number; optional)

- className (string; optional)

- file_info (dict; optional)

- progress (number; optional)"""
    _children_props: typing.List[str] = []
    _base_nodes = ['children']
    _namespace = 'chunk_upload'
    _type = 'ChunkUpload'


    def __init__(
        self,
        id: typing.Optional[typing.Union[str, dict]] = None,
        file_info: typing.Optional[dict] = None,
        progress: typing.Optional[NumberType] = None,
        chunk_size: typing.Optional[NumberType] = None,
        className: typing.Optional[str] = None,
        style: typing.Optional[typing.Any] = None,
        **kwargs
    ):
        self._prop_names = ['id', 'chunk_size', 'className', 'file_info', 'progress', 'style']
        self._valid_wildcard_attributes =            []
        self.available_properties = ['id', 'chunk_size', 'className', 'file_info', 'progress', 'style']
        self.available_wildcard_properties =            []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs and excess named props
        args = {k: _locals[k] for k in _explicit_args}

        super(ChunkUpload, self).__init__(**args)

setattr(ChunkUpload, "__init__", _explicitize_args(ChunkUpload.__init__))
