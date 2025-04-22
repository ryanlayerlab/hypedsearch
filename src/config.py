from abc import ABC, abstractmethod
from dataclasses import MISSING, fields
from typing import Callable, Generic, TypeVar, get_args

import click

T = TypeVar("T")
DEFAULT_HELP_STRING = "Command-line interface for app."


class Required(Generic[T]):
    """Type annotation for required command-line arguments."""

    def __init__(self, value: T):
        self.value = value

    def __repr__(self) -> str:
        return f"Required({self.value})"


class AppConfig(ABC):
    """Template class for application configuration."""

    @staticmethod
    @abstractmethod
    def get_variable_help_strings():
        pass

    @staticmethod
    @abstractmethod
    def get_variable_flags():
        pass


def process_dataclass_field_type(field_type: type) -> click.types:
    """Convert Python types to Click types."""
    required = False
    # Support required types
    if hasattr(field_type, "__origin__") and field_type.__origin__ is Required:
        field_type = get_args(field_type)[0]
        required = True
    if field_type is bool:
        field_type = click.BOOL
    elif field_type is int:
        field_type = click.INT
    elif field_type is float:
        field_type = click.FLOAT
    elif field_type is str:
        field_type = click.STRING
    else:
        raise ValueError(f"Unsupported type: {field_type}")

    return field_type, required


def generate_cli(
    cls: AppConfig,
    main_fcn: Callable,
    help_str: str = DEFAULT_HELP_STRING,
    name: str = "MyApp",
):
    context_settings = {"help_option_names": ["-h", "--help"], "max_content_width": 200}

    @click.command(name=name, help=help_str, context_settings=context_settings)
    @click.option(
        "--config_path",
        "-c",
        type=click.Path(exists=True),
        help="Path to a JSON config file. ",
    )
    def cli(config_path, **kwargs):
        """Function that receives CLI arguments."""
        # Get values from CLI options or defaults
        config = cls(**kwargs)

        # Validate the config

        # Run the main logic of the application
        main_fcn(config)

    # Dynamically add click options based on dataclass fields
    for field in fields(cls):
        field_type, required = process_dataclass_field_type(field_type=field.type)
        field_name = field.name
        default = None if field.default is MISSING else field.default
        if field_name in cls.get_variable_flags():
            cli = click.option(
                f"--{field_name}",
                f"-{cls.get_variable_flags()[field_name]}",
                required=required,
                default=default,
                show_default=True,
                type=field_type,
                help=cls.get_variable_help_strings().get(field_name, ""),
            )(cli)
        else:
            cli = click.option(
                f"--{field.name}",
                required=required,
                default=default,
                show_default=True,
                type=field_type,
                help=cls.get_variable_help_strings().get(field.name, ""),
            )(cli)

    return cli
