from dataclasses import dataclass

from src.sql_database import (
    SqlColumn,
    SqlColumnTypes,
    Sqlite3Database,
    add_data_classes_to_table,
    get_create_table_query_for_dataclass,
    get_sql_columns_from_dataclass,
)


class TestGetCreateTableQueryForDataclass:
    @staticmethod
    def test_with_primary_key():
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        name = "products"
        actual = get_create_table_query_for_dataclass(
            obj=Product, table_name=name, primary_key_colm="name"
        )
        expected = "CREATE TABLE IF NOT EXISTS products (id INTEGER, name TEXT PRIMARY KEY, price REAL)"
        assert actual == expected

    @staticmethod
    def test_make_table_from_query():
        # Arrange
        @dataclass
        class Name:
            name: str

        table_name, db_path = "names", ":memory:"

        # Act
        db = Sqlite3Database(db_path=db_path)
        query = get_create_table_query_for_dataclass(obj=Name, table_name=table_name)
        db.execute_query(query=query)
        assert db.tables() == [table_name]

    @staticmethod
    def test_without_primary_key():
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        name = "products"
        actual = get_create_table_query_for_dataclass(obj=Product, table_name=name)
        expected = (
            "CREATE TABLE IF NOT EXISTS products (id INTEGER, name TEXT, price REAL)"
        )
        assert actual == expected


class TestGetSqlColumnsFromDataclass:
    @staticmethod
    def product_dataclass():
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        return Product

    @staticmethod
    def test_no_primary_key():
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        expected = [
            SqlColumn(name="id", type=SqlColumnTypes.SQL_INT),
            SqlColumn(name="name", type=SqlColumnTypes.SQL_TEXT),
            SqlColumn(name="price", type=SqlColumnTypes.SQL_REAL),
        ]
        actual = get_sql_columns_from_dataclass(obj=Product)
        assert actual == expected

    @staticmethod
    def test_with_primary_key():
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        expected = [
            SqlColumn(name="id", type=SqlColumnTypes.SQL_INT),
            SqlColumn(name="name", type=SqlColumnTypes.SQL_TEXT, primary_key=True),
            SqlColumn(name="price", type=SqlColumnTypes.SQL_REAL),
        ]
        actual = get_sql_columns_from_dataclass(obj=Product, primary_key_colm="name")
        assert actual == expected


class TestAddDataClassesToTable:
    @staticmethod
    def test_smoke():
        # Arrange
        @dataclass
        class Product:
            id: int
            name: str
            price: float

        table_name, db_path = "products", ":memory:"
        dataclasses = [
            Product(id=0, name="1", price=2),
            Product(id=1, name="2", price=3),
        ]

        # Act
        # Create table and add dataclasses
        db = Sqlite3Database(db_path=db_path)
        query = get_create_table_query_for_dataclass(
            obj=dataclasses[0], table_name=table_name, primary_key_colm="name"
        )
        db.execute_query(query=query)
        add_data_classes_to_table(dataclasses=dataclasses, table_name=table_name, db=db)

        # Assert
        rows = db.all_table_rows(table_name=table_name)
        rows = [Product(**row) for row in rows]
        assert dataclasses == rows
