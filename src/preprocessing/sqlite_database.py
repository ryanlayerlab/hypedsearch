import sqlite3
import time
import os
import sys
import shutil
import computational_pipeline

class Sqllite_Database:
    def __init__(self, max_len, reset=True):
        self.connection = sqlite3.connect("hypedsearch.db")
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
            self.cursor.execute("DROP TABLE IF EXISTS proteins")
            self.cursor.execute('''CREATE TABLE proteins (
                                    id INTEGER,
                                    description TEXT,
                                    sequence TEXT
                                )''')            
        
    def insert_kmers(self,data):
        self.cursor.executemany('INSERT INTO kmers VALUES(?, ?, ?, ?, ?, ?)', data)
        self.connection.commit()
        
    def insert_proteins(self,data):
        self.cursor.executemany('INSERT INTO proteins VALUES(?, ?, ?)', data)
        self.connection.commit()

    def read_kmers(self):
        rows = self.cursor.execute("SELECT * FROM kmers").fetchall()
        print(rows)
    
    def query_mass_kmers(self, mass, tol):
        upper = mass + tol
        lower = mass - tol
        self.cursor.execute("CREATE TABLE temp.mass AS SELECT * FROM kmers where mass between ? and ? order by protein, location_start", (lower, upper))
        b_rows = self.cursor.execute("SELECT * FROM temp.mass where ion = 0 order by protein, location_start, location_end").fetchall()
        y_rows = self.cursor.execute("SELECT * FROM temp.mass where ion = 1 order by protein, location_start, location_end").fetchall()
        self.cursor.execute("DROP TABLE temp.mass")
        return b_rows, y_rows
      
    def query_sequence_kmers(self, pid, start, end):
        rows = self.cursor.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end = ?", (pid, start, end)).fetchall()
        return rows
    
    def query_extensions_to_length_kmers(self, target_mass, pid, start, end, dist):
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
        
    def query_extensions_b_kmers(self, target_mass, pid, start, end, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute("SELECT * FROM kmers where protein = ? and location_start = ? and location_end >= ? and mass < ? and ion = ? and charge = 2 order by location_end", (pid, start, end, target_mass, ion))
        query_time = time.time() - query_start
        self.query_protein_average = (self.query_protein_average * self.protein_count + query_time)/ (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor
    
    def query_extensions_y_kmers(self, target_mass, pid, end, start, ion):
        query_start = time.time()
        rows_cursor = self.connection.execute("SELECT * FROM kmers where protein = ? and location_end = ? and location_start <= ? and mass < ? and ion = ? and charge = 2 order by location_start desc", (pid, end, start, target_mass, ion))
        query_time = time.time() - query_start
        self.query_protein_average = (self.query_protein_average * self.protein_count + query_time)/ (self.protein_count + 1)
        self.protein_count += 1
        return rows_cursor
    
    def query_fetchall(self, cursor):
        return cursor.fetchall()
    
    def index_ion_mass_b_kmers(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX b_mass_ion_idx ON kmers(protein, location_start)")
        print("time to index by protein", time.time() - ctime)
    
    def index_ion_mass_y_kmers(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX y_mass_ion_idx ON kmers(protein, location_end)")
        print("time to index by protein", time.time() - ctime)
        
    def index_ion_mass_kmers(self):
        ctime = time.time()
        self.cursor.execute("CREATE INDEX mass_ion_idx ON kmers(mass, protein, location_start)")
        print("time to index", time.time() - ctime)
        
    def check_sizes_kmers(self):
        size = self.cursor.execute("SELECT count(*) FROM kmers").fetchall()
        print(size)
        
    def run_query(self, query):
        results = self.cursor.execute(query).fetchall()
        print(results)
        
    def count_ion_mass_kmers(self, mass, tol, ion):
        upper = mass + tol
        lower = mass - tol
        qtime = time.time()
        rows = self.cursor.execute("SELECT count(*) FROM kmers where ion = ? and mass between ? and ?", (ion, lower, upper)).fetchall()
        print("Query time:", time.time() - qtime)
        print(mass, ion, rows)
        return rows
    
    def query_proteins(self):
        rows = self.cursor.execute("SELECT * FROM proteins").fetchall()
        return rows
    
    def get_protein(self, pid):
        row = self.cursor.execute("SELECT * FROM proteins WHERE id = ?", (pid,)).fetchone()        
        return row       
    
#     def get_proteins_with_subsequence(db: Database, sequence: str) -> list:
#         return list(set(db.kmers[sequence]))

#     def get_proteins_with_subsequence_ion(db: Database, sequence: str, ion: str) -> list:
#         hits = []
#         subseq = sequence
#         while len(hits) == 0 and len(subseq) > 0:
#             hs = get_proteins_with_subsequence(db, subseq)
#             for h in hs:
#                 for entry in get_entry_by_name(db, h):
#                     if sequence in entry.sequence:
#                         hits.append(h)
#             if len(hits) == 0:
#                 subseq = subseq[1:] if ion == 'y' else subseq[:-1]
#         return hits

#     def get_entry_by_name(db: Database, name: str) -> namedtuple:
#         return db.proteins[name]    
    

    def get_kmers_for_protein(self, kmer, start, end, protein_id, ion):
        data_list = []
        for charge in [1,2]:
            mass = computational_pipeline.gen_spectra.get_max_mass(kmer, ion=ion, charge=charge)
            ion_int = 0 if ion == 'b' else 1
            input_tuple = (mass, start, end, ion_int, charge, protein_id)
            data_list.append(input_tuple)
        return data_list
    
    def db_make_set_for_protein_digest(self,protein_id,protein,max_peptide_length, digest_left, digest_right):
        data = []
        seq_len = len(protein)
        count_max = 1000000
        for size in range(2, max_peptide_length + 1):
            for start in range(0, seq_len - size + 1):
                end = start + size
                kmer = protein[start:end]
                if kmer[0] in digest_left or digest_left == ['-'] or (start > 0 and protein[start-1] in digest_right):
                    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                    if not any (x in bad_chars for x in kmer):
                        data_list = self.get_kmers_for_protein(kmer, start, end, protein_id, 'b')
                        data.extend(data_list)
                        if len(data) > count_max:
                            self.insert(data)
                            data.clear()
                if kmer[-1] in digest_right or digest_right == ['-'] or (end < seq_len and protein[end] in digest_left):
                    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                    if not any (x in bad_chars for x in kmer):
                        data_list = self.get_kmers_for_protein(kmer, start, end, protein_id, 'y')
                        data.extend(data_list)
                        if len(data) > count_max:
                            self.insert(data)
                            data.clear()
        return data
    
    def check_for_enough_disk_space(self, protein_id, plen, last_percent):
        percent = int((protein_id + 1) * 100 / plen)
        print(f'\rOn protein {protein_id + 1}/{plen} [{percent}%]', end='')
        if percent != last_percent:
            last_percent = percent
            free = shutil.disk_usage('/')[2]
            free = free / (1024 ** 3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB")
                return False, last_percent
            else:
                return True, last_percent
        return True, last_percent

    def insert_prepped_kmers(self, kv_proteins, max_peptide_length, digest_left, digest_right):
        plen = len(kv_proteins)
        last_percent = 0
        all_data = []
        
        for protein_id, (_, protein) in enumerate(kv_proteins):
            enough_space, last_percent = self.check_for_enough_disk_space(protein_id, plen, last_percent)
            if enough_space:
                data = self.db_make_set_for_protein_digest(protein_id,protein,max_peptide_length, digest_left, digest_right)
                all_data.extend(data)
            else:
                sys.exit(1)

        if len(all_data) != 0:
            self.insert_kmers(all_data)

    def insert_prepped_proteins(self,kv_proteins):
        all_data = []
        protein_id = 0
        for kv_protein in kv_proteins:
            (description,amino_acids) = kv_protein
            input_tuple = (protein_id,description,amino_acids)
            all_data.append(input_tuple)
            protein_id = protein_id + 1
        self.insert_proteins(all_data)

    def populate_database(self, kv_proteins, max_peptide_length, digest_left, digest_right):
        self.insert_prepped_proteins(kv_proteins)
        self.insert_prepped_kmers(kv_proteins, max_peptide_length, digest_left, digest_right)
        self.index_ion_mass_kmers()
        self.index_ion_mass_b_kmers()
        self.index_ion_mass_y_kmers()
