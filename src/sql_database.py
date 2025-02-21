import sqlite3
from dataclasses import asdict, fields, is_dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel

SQL_PRIMARY_KEY = "PRIMARY KEY"


class SqlColumnTypes(Enum):
    SQL_INT = "INTEGER"
    SQL_TEXT = "TEXT"
    SQL_REAL = "REAL"


class SqlColumn(BaseModel):
    name: str
    type: SqlColumnTypes
    primary_key: Optional[bool] = False

    def __str__(self):
        val = f"{self.name} {self.type.value}"
        if self.primary_key is True:
            val += f" {SQL_PRIMARY_KEY}"
        return val


def get_sql_field_type(field_type: type):
    if field_type == int:
        return SqlColumnTypes.SQL_INT.value
    elif field_type == str:
        return SqlColumnTypes.SQL_TEXT.value
    elif field_type == float:
        return SqlColumnTypes.SQL_REAL.value
    else:
        raise RuntimeError(f"Field type {field_type} is not supported")


def get_sql_columns_from_dataclass(
    obj: Any, primary_key_colm: Optional[str] = None
) -> List[SqlColumn]:
    colms = []
    for field in fields(obj):
        primary_key = False
        if field.name == primary_key_colm:
            primary_key = True
        colms.append(
            SqlColumn(
                name=field.name,
                type=get_sql_field_type(field_type=field.type),
                primary_key=primary_key,
            )
        )
    return colms


def get_create_table_query_for_dataclass(
    obj: Any, table_name: str, primary_key_colm: Optional[str] = None
) -> str:
    assert is_dataclass(obj), f"{obj} isn't a dataclass"
    colms = get_sql_columns_from_dataclass(obj=obj, primary_key_colm=primary_key_colm)
    colms_for_query = ", ".join([str(colm) for colm in colms])
    return f"CREATE TABLE IF NOT EXISTS {table_name} ({colms_for_query})"


class Sqlite3Database:
    def __init__(self, db_path: str):
        self.connection = sqlite3.connect(db_path)
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


def add_data_classes_to_table(
    dataclasses: List[Any], table_name: str, db: Sqlite3Database
):
    # Create string ':<column 1>, :<column 2>, ...' that's used in the INSERT SQL query
    sql_colms = get_sql_columns_from_dataclass(obj=dataclasses[0])
    values = ""
    for colm in sql_colms:
        values += f":{colm.name}, "
    values = values[:-2]  # remove traling ', '

    # Insert dataclasses into the table
    dataclasses_as_dicts = [asdict(x) for x in dataclasses]
    _ = db.cursor.executemany(
        f"INSERT INTO {table_name} VALUES({values})", dataclasses_as_dicts
    )
    _ = db.connection.commit()
