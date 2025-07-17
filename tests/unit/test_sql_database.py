from dataclasses import dataclass, fields

from src.sql_database import (
    SQL_PRIMARY_KEY,
    PrimaryKey,
    SqlColumn,
    SqlColumnTypes,
    Sqlite3Database,
    SqlTableRow,
)


def get_columns_for_create_table_query_from_dataclass():
    pass


class Test_SqlColumn_initialization:
    @staticmethod
    def test_direct_not_primary_key():
        name, py_type = "a", int
        sql_type = SqlColumnTypes.SQL_INT
        actual = SqlColumn(name=name, py_type=py_type)
        assert actual.primary_key is False
        assert actual.py_type is py_type
        assert actual.name == name
        assert actual.sql_type == sql_type

    @staticmethod
    def test_direct_primary_key():
        name, py_type = "a", PrimaryKey[str]
        sql_type = SqlColumnTypes.SQL_TEXT
        sql_colm = SqlColumn(name=name, py_type=py_type)
        sql_colm
        assert sql_colm.primary_key is True
        assert sql_colm.py_type is py_type
        assert sql_colm.name == name
        assert sql_colm.sql_type == sql_type

    @staticmethod
    def test_from_dataclass_field():
        @dataclass
        class Product:
            id: PrimaryKey[int]
            val: float

        sql_colms = [
            SqlColumn.from_datclass_field(field=field) for field in fields(Product)
        ]

        sql_colm = sql_colms[0]
        assert sql_colm.primary_key is True
        assert sql_colm.name == "id"
        assert sql_colm.sql_type == SqlColumnTypes.SQL_INT

        sql_colm = sql_colms[1]
        assert sql_colm.primary_key is False
        assert sql_colm.name == "val"
        assert sql_colm.sql_type == SqlColumnTypes.SQL_REAL


class Test_SqlColumn_sql_column_definition:
    @staticmethod
    def test_no_primary_key():
        name, py_type = "a", int
        sql_type = SqlColumnTypes.SQL_INT
        expected_defn = f"{name} {sql_type.value}"
        sql_colm = SqlColumn(name=name, py_type=py_type)
        colm_defn = sql_colm.sql_column_definition()
        assert colm_defn == expected_defn

    @staticmethod
    def test_primary_key():
        name, py_type = "a", PrimaryKey[float]
        sql_type = SqlColumnTypes.SQL_REAL
        expected_defn = f"{name} {sql_type.value} {SQL_PRIMARY_KEY}"
        sql_colm = SqlColumn(name=name, py_type=py_type)
        colm_defn = sql_colm.sql_column_definition()
        assert colm_defn == expected_defn


class Test_SqlTableRow_sql_columns:
    @staticmethod
    def test_smoke():
        @dataclass
        class ProteinTest(SqlTableRow):
            id: PrimaryKey[int]
            name: str

        sql_colms = ProteinTest.sql_columns()
        # assert protein._sql_colms is not None
        assert len(sql_colms) == 2


class Test_SqlTableRow_sql_create_table_columns_str:
    @staticmethod
    def test_smoke():
        @dataclass
        class ProteinTest(SqlTableRow):
            id: PrimaryKey[int]
            name: str

        colm_str = ProteinTest.sql_create_table_columns_str()
        expected_str = "id INTEGER PRIMARY KEY, name TEXT"
        assert colm_str == expected_str


class Test_SqlTableRow_sql_insert_columns_str:
    @staticmethod
    def test_smoke():
        @dataclass
        class ProteinTest(SqlTableRow):
            id: PrimaryKey[int]
            name: str

        colm_str = ProteinTest.sql_insert_columns_str()
        expected_str = ":id, :name"
        assert colm_str == expected_str


class Test_Sqlite3Database_create_table_from_dataclass:
    @staticmethod
    def test_smoke():
        @dataclass
        class ProteinTest(SqlTableRow):
            id: PrimaryKey[int]
            seq: str

        table_name, obj = "proteins", ProteinTest

        db = Sqlite3Database()
        db.create_table_from_dataclass(table_name=table_name, obj=obj)
        assert db.tables() == [table_name]
        table_colms = db.table_info(table_name=table_name)
        assert table_colms[0]["name"] == "id"
        assert table_colms[0]["type"] == "INTEGER"
        assert table_colms[0]["pk"] == 1
        assert table_colms[1]["name"] == "seq"
        assert table_colms[1]["type"] == "TEXT"
        assert table_colms[1]["pk"] == 0


class Test_Sqlite3Database_insert_dataclasses:
    @staticmethod
    def test_smoke():
        @dataclass
        class ProteinTest(SqlTableRow):
            id: PrimaryKey[int]
            seq: str
            mass: float

        table_name, obj = "proteins", ProteinTest
        proteins = [
            ProteinTest(id=0, seq="ABC", mass=1.2),
            ProteinTest(id=1, seq="EFG", mass=0.5),
        ]
        db = Sqlite3Database()
        db.create_table_from_dataclass(table_name=table_name, obj=obj)
        db.insert_dataclasses(table_name=table_name, data_classes=proteins)
        rows = db.all_table_rows(table_name=table_name)
        rows = [ProteinTest(**row) for row in rows]
        assert rows == proteins
