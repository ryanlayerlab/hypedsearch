import sqlite3
import time

class database_file:
    def __init__(self, max_len, reset=True):
        self.connection = sqlite3.connect("kmers.db")
        self.cursor = self.connection.cursor()
        self.max_len = max_len
        self.query_protein_average = 0
        self.query_mass_average = 0
        self.protein_count = 0
        self.mass_count = 0
        self.fetchall_protein_average = 0
        self.fetchall_mass_average = 0
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
    
    def query_mass(self, mass, tol):
        upper = mass + tol
        lower = mass - tol
        self.cursor.execute("CREATE TABLE temp.mass AS SELECT * FROM kmers where mass between ? and ? order by protein, location_start", (lower, upper))
        b_rows = self.cursor.execute("SELECT * FROM temp.mass where ion = 0 order by protein, location_start, location_end").fetchall()
        y_rows = self.cursor.execute("SELECT * FROM temp.mass where ion = 1 order by protein, location_start, location_end").fetchall()
        self.cursor.execute("DROP TABLE temp.mass")
        return b_rows, y_rows
    
    # def query_unique_mass(self, mass, tol):
        # upper = mass + tol
        # lower = mass - tol
        # self.cursor.execute("CREATE TABLE temp.mass AS SELECT * FROM kmers where mass between ? and ? order by protein, location_start", (lower, upper))
        # unique_b_rows = self.cursor.execute("SELECT DISTINCT location_end - location_start FROM temp.mass where ion = 0 order by protein, location_start, location_end").fetchall()
        # 
        
    def query_sequence(self, pid, start, end):
        rows = self.cursor.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end = ?", (pid, start, end)).fetchall()
        return rows
    
    def query_extensions_to_length(self, target_mass, pid, start, end, dist):
        query_start = time.time()
        rows_cursor = self.cursor.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end = ? and mass <= ? and ion = 0 and charge = 2", (pid, start, end+dist-1, target_mass)) #and ion = 0 and charge = 1
        query_time = time.time() - query_start
        fetchall_start = time.time()
        rows = rows_cursor.fetchall()
        fetchall_time = time.time() - fetchall_start
        self.query_protein_average = (self.query_protein_average * self.protein_count + query_time)/ (self.protein_count + 1)
        self.fetchall_protein_average = (self.fetchall_protein_average * self.protein_count + fetchall_time)/ (self.protein_count + 1)
        self.protein_count += 1
        return rows
        
    def query_extensions_b(self, target_mass, pid, start, end, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end >= ? and mass < ? and ion = ? and charge = 2 order by location_end", (pid, start, end, target_mass, ion))
        # rows_cursor = self.cursor.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end > ? and mass < ? and ion = ? and charge = 2 order by location_end", (pid, start, end, target_mass, ion))
        query_time = time.time() - query_start
        self.query_protein_average = (self.query_protein_average * self.protein_count + query_time)/ (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor
    
    def query_extensions_y(self, target_mass, pid, end, start, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute("SELECT * FROM kmers where protein = ? and location_end = ? and location_start <= ? and mass < ? and ion = ? and charge = 2 order by location_start desc", (pid, end, start, target_mass, ion))
        # rows_cursor = self.cursor.execute("SELECT * FROM kmers where protein = ? and location_end = ? and location_start < ? and mass < ? and ion = ? and charge = 2 order by location_start desc", (pid, end, start, target_mass, ion))
        query_time = time.time() - query_start
        self.query_protein_average = (self.query_protein_average * self.protein_count + query_time)/ (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor
    
    def query_fetchall(self, cursor):
        return cursor.fetchall()
    
    def index_ion_mass_b(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX b_mass_ion_idx ON kmers(protein, location_start)")
        print("time to index by protein", time.time() - ctime)
    
    def index_ion_mass_y(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX y_mass_ion_idx ON kmers(protein, location_end)")
        print("time to index by protein", time.time() - ctime)
        
    def index_ion_mass(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX mass_ion_idx ON kmers(mass, protein, location_start)")
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