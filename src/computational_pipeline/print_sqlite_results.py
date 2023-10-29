
# Import the python module - sqlite3
import sqlite3

# Create database connection to a main database which is stored in a file
connectionObject = sqlite3.connect("kmers.db")

# Obtain a cursor object
cursorObject = connectionObject.cursor()

# Print the tables and indices present in the SQLite main database
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