import sqlite3
from preprocessing import merge_search
import time

from matplotlib.pyplot import connect

class database_file:
    def __init__(self, max_len): #This is like the named tuple
        self.connection = sqlite3.connect("kmers.db")
        self.cursor = self.connection.cursor()
        self.max_len = max_len
        self.cursor.execute("DROP TABLE IF EXISTS kmers")
        self.cursor.execute('''CREATE TABLE kmers (
                                mass REAL,
                                kmer TEXT,
                                location_start INTEGER,
                                location_end INTEGER,
                                ion INTEGER,
                                charge INTEGER,
                                protein INTEGER
                            )''')
        
    def insert(self,data):
        # print(data)
        self.cursor.executemany('INSERT INTO kmers VALUES(?, ?, ?, ?, ?, ?, ?)', data)
        # print(self.connection.total_changes)
        self.connection.commit()
        
    def read(self):
        rows = self.cursor.execute("SELECT * FROM kmers").fetchall()
        print(rows)
    
    def query_mass_ion(self, mass, tol, ion):
        upper = mass + tol
        lower = mass - tol
        rows = self.cursor.execute("SELECT * FROM kmers where mass between ? and ? and ion = ?", (lower, upper, ion)).fetchall()
        # print(mass, ion, rows)
        return rows
            
    def index_mass_ion(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX mass_ion_idx ON kmers(mass, ion)")
        print("time", time.time() - ctime)
    
    def __call__(self, proteins):
        # mass [(kmer string, location int, ion boolean, charge boolean), ... ]
        # (mass float, kmer string, location_start int, location_end int, ion boolean, charge boolean), mass is indexed
        # 1. Create a database table with name and columns (mass float, kmer string, location_start int, location_end int, ion boolean, charge boolean)
        # 
            
            
        # data = [
        #     ('1.2', 'ASD', 0, 1000, 0, 1, 45),
        #     ('6.7', 'HJK', 85, 682, 1, 0, 35),
        # ]
            seq_len = len(prot.sequence)
            all_data = []
            for size in range(1, max_len + 1):
                # size -> [1, max_len]
                for start in range(0, seq_len - size + 1):
                    end = start + size
                    kmer = prot.sequence[start:end]
                    # last_index = seq - size 6, end = start + size - 1 = 7
                    # [data.append(x) for x in get_data(kmer, start, end)]
                    print(kmer)
                    data = merge_search.get_data(kmer, start, end, i)
                    all_data.append(data)
                    
            self.connection.commit()
    