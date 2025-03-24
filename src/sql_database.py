import logging
import sqlite3
import time
from dataclasses import Field, asdict, dataclass, field, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Generic, List, Optional, TypeVar, get_args, get_origin

from pydantic import BaseModel

from src.constants import MEMORY
from src.utils import get_time_in_diff_units

logger = logging.getLogger(__name__)


SQL_PRIMARY_KEY = "PRIMARY KEY"
T = TypeVar("T")


class PrimaryKey(Generic[T]):
    """Type annotation for primary keys."""

    def __init__(self, value: T):
        self.value = value

    def __repr__(self) -> str:
        return f"PrimaryKey({self.value})"


class SqlColumnTypes(Enum):
    SQL_INT = "INTEGER"
    SQL_TEXT = "TEXT"
    SQL_REAL = "REAL"


@dataclass
class SqlColumn:
    name: str
    py_type: type
    primary_key: bool = field(init=False)
    sql_type: SqlColumnTypes = field(init=False)
    allowed_types = [int, float, str]

    def __post_init__(self):
        # Parse column type
        # Support
        field_type = self.py_type
        primary_key = False

        # Support PrimaryKey types
        if hasattr(field_type, "__origin__") and field_type.__origin__ is PrimaryKey:
            field_type = get_args(field_type)[0]
            primary_key = True

        if field_type == int:
            sql_type = SqlColumnTypes.SQL_INT
        elif field_type == str:
            sql_type = SqlColumnTypes.SQL_TEXT
        elif field_type == float:
            sql_type = SqlColumnTypes.SQL_REAL
        else:
            raise RuntimeError(
                f"Field type {field_type} is not supported. Allowed types are {cls.allowed_types}"
            )
        self.sql_type = sql_type
        self.primary_key = primary_key

    def sql_column_definition(self):
        val = f"{self.name} {self.sql_type.value}"
        if self.primary_key is True:
            val += f" {SQL_PRIMARY_KEY}"
        return val

    @classmethod
    def from_datclass_field(cls, field: Field) -> "SqlColumn":
        return cls(
            name=field.name,
            py_type=field.type,
        )


@dataclass
class SqlTableRow:
    @classmethod
    def sql_columns(cls) -> List[SqlColumn]:
        sql_colms = []
        for f in fields(cls):
            # Skip hidden fields
            if f.name[0] == "_":
                continue
            sql_colm = SqlColumn.from_datclass_field(field=f)
            sql_colms.append(sql_colm)
        return sql_colms

    @classmethod
    def sql_create_table_columns_str(cls):
        sql_colms = cls.sql_columns()
        return ", ".join([colm.sql_column_definition() for colm in sql_colms])

    @classmethod
    def sql_insert_columns_str(cls):
        sql_colms = cls.sql_columns()
        colm_str = ""
        for colm in sql_colms:
            colm_str += f":{colm.name}, "
        colm_str = colm_str[:-2]  # remove traling ', '
        return colm_str


class Sqlite3Database:
    def __init__(self, db_path: Optional[str] = MEMORY, overwrite: bool = False):
        self.db_path = db_path
        if self.db_path != MEMORY:
            self.db_path = Path(db_path).absolute()
            if self.db_path.exists():
                logger.info(f"The database already exists at {self.db_path}.")
                if overwrite:
                    logger.info(
                        "'overwrite' is True so deleting the existing databse and starting anew."
                    )
                    self.db_path.unlink()
                else:
                    logger.info(
                        "'overwrite' is False so connecting to the existing one."
                    )
        self.connection = sqlite3.connect(self.db_path)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()

    def execute_query(self, query: str) -> None:
        self.cursor.execute(query)
        self.connection.commit()

    def read_query(self, query: str) -> List[Dict]:
        results = self.cursor.execute(query).fetchall()
        return [dict(x) for x in results]

    def tables(self):
        name = "name"
        query = f"SELECT {name} FROM sqlite_master WHERE type = 'table';"
        rows = self.read_query(query=query)
        return [dict(row)[name] for row in rows]

    def indices(self):
        name = "name"
        query = f"SELECT {name} FROM sqlite_master WHERE type='index';"
        rows = self.read_query(query=query)
        return [row[name] for row in rows]

    def all_table_rows(self, table_name: str):
        query = f"SELECT * FROM {table_name}"
        return self.read_query(query=query)

    def table_info(self, table_name: str):
        query = f"PRAGMA table_info({table_name})"
        return self.read_query(query=query)

    def row_count(self, table_name: str):
        query = f"SELECT COUNT(*) FROM {table_name}"
        return self.read_query(query=query)

    def add_index(self, table_name: str, index_name: str, colms_to_index: List[str]):
        start_time = time.time()
        colms = ", ".join(colms_to_index)
        query = f"CREATE INDEX {index_name} ON {table_name}({colms})"
        self.execute_query(query=query)
        logger.info(
            f"Indexing took {get_time_in_diff_units(time_sec=(time.time() - start_time))}"
        )

    def create_table_from_dataclass(self, table_name: str, obj: type[SqlTableRow]):
        colm_str = obj.sql_create_table_columns_str()
        query = f"CREATE TABLE IF NOT EXISTS {table_name} ({colm_str})"
        self.execute_query(query=query)
        logger.info(f"Created table {table_name}")

    def insert_dataclasses(
        self, table_name: str, data_classes: List[SqlTableRow]
    ) -> None:

        colm_str = data_classes[0].sql_insert_columns_str()
        # Insert dataclasses into the table
        dataclasses_as_dicts = [asdict(x) for x in data_classes]
        self.cursor.executemany(
            f"INSERT INTO {table_name} VALUES({colm_str})", dataclasses_as_dicts
        )
        self.connection.commit()


def get_dataclass_fields_as_sql_columns():
    pass


def column_string_for_create_table_query_from_dataclass():
    pass


def add_data_classes_to_table(
    dataclasses: List[Any], table_name: str, db: Sqlite3Database
) -> None:
    pass
