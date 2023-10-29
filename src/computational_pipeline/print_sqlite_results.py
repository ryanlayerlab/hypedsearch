
import sqlite3
connectionObject = sqlite3.connect("kmers.db")
cursorObject = connectionObject.cursor()
cursorObject.execute("select * from SQLite_master")
tables = cursorObject.fetchall()
print("Listing tables and indices from main database:")

for table in tables:
        print("Type of database object: %s"%(table[0]))
        print("Name of the database object: %s"%(table[1]))
        print("Table Name: %s"%(table[2]))
        print("Root page: %s"%(table[3]))
        print("SQL statement: %s"%(table[4]))

connectionObject.close()