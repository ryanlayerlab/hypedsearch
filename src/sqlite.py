import sqlite3
import time

from matplotlib.pyplot import connect

class database_file:
    def __init__(self, max_len, reset=True):
        self.connection = sqlite3.connect("kmers.db")
        self.cursor = self.connection.cursor()
        self.max_len = max_len
        if reset:
            self.cursor.execute("DROP TABLE IF EXISTS kmers")
            self.cursor.execute('''CREATE TABLE kmers (
                                    mass REAL,
                                    location_start INTEGER,
                                    location_end INTEGER,
                                    ion INTEGER,
                                    charge INTEGER,
                                    protein INTEGER
                                )''')
        
    def insert(self,data):
        self.cursor.executemany('INSERT INTO kmers VALUES(?, ?, ?, ?, ?, ?)', data)
        self.connection.commit()
        
    def read(self):
        rows = self.cursor.execute("SELECT * FROM kmers").fetchall()
        print(rows)
    
    def query_ion_mass(self, mass, tol, ion):
        upper = mass + tol
        lower = mass - tol
        rows = self.cursor.execute("SELECT * FROM kmers where ion = ? and mass between ? and ?", (ion, lower, upper)).fetchall()
        return rows
            
    def index_ion_mass(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX mass_ion_idx ON kmers(ion, mass)")
        print("time to index", time.time() - ctime)
        
    def check_sizes(self):
        size = self.cursor.execute("SELECT count(*) FROM kmers").fetchall()
        print(size)
        
    def run_query(self, query):
        results = self.cursor.execute(query).fetchall()
        print(results)
        
    def count_ion_mass(self, mass, tol, ion):
        upper = mass + tol
        lower = mass - tol
        qtime = time.time()
        rows = self.cursor.execute("SELECT count(*) FROM kmers where ion = ? and mass between ? and ?", (ion, lower, upper)).fetchall()
        print("Query time:", time.time() - qtime)
        print(mass, ion, rows)
        return rows

    def check_indices():
        pass