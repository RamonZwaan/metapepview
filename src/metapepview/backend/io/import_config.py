from typing import Any

import tomli_w
import tomllib
import base64

from metapepview.backend.utils.io_utils import memory_to_str


def config_to_toml(config: dict[str, Any]) -> str:
    """Serialize configuration fields to TOML.

    Supported values are booleans, strings, integers, and floats.
    Fields with None values are omitted.
    """
    supported_types = (bool, str, int, float)

    filtered_config = {}

    for key, value in config.items():
        if value is None:
            continue

        if type(value) not in supported_types:
            raise TypeError(
                f"Unsupported value type for field {key!r}: "
                f"{type(value).__name__}"
            )

        filtered_config[key] = value

    return tomli_w.dumps(filtered_config)


def parse_toml_upload(contents: str) -> dict[str, Any]:
    """Decode and parse the base64 contents supplied by dcc.Upload."""
    if not contents:
        raise ValueError("No file contents received.")

    try:
        toml_text = memory_to_str(contents)
        data: dict[str, Any] = tomllib.loads(toml_text)
    except (ValueError, UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"Invalid TOML file: {exc}") from exc

    return data
