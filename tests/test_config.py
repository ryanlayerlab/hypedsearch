import re
from dataclasses import dataclass, fields

import click
from click.testing import CliRunner

from config import AppConfig, Required, generate_cli, process_dataclass_field_type


class Test_generate_cli:
    @staticmethod
    def test_smoke():
        @dataclass
        class TestConfig(AppConfig):
            a: Required[int]
            b: str
            c: int = 1

            @staticmethod
            def get_variable_help_strings():
                return {}

            @staticmethod
            def get_variable_flags():
                return {}

        def process_config(config):
            """Use the config values in a later function."""
            click.echo(f"Configuration: {config}")

        # def test_main(config: TestConfig):
        #     print("a=")
        cli = generate_cli(cls=TestConfig, main_fcn=process_config)
        runner = CliRunner()
        result = runner.invoke(cli, ["-h"])
        assert "Usage: MyApp" in result.output
        assert re.search(r"--a\s+INTEGER\s*\[required\]", result.output)


class Test_process_dataclass_field_type:
    @staticmethod
    def test_required_type():
        @dataclass
        class Object:
            a: Required[int]
            b: str

        f = fields(Object)[0]
        field_type, required = process_dataclass_field_type(field_type=f.type)

        assert field_type == click.INT
        assert required is True

        f = fields(Object)[1]
        field_type, required = process_dataclass_field_type(field_type=f.type)

        assert field_type == click.STRING
        assert required is False
