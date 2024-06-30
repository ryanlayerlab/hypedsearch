import computational_pipeline.gen_spectra
import shutil
import sys

def get_data(kmer, start, end, protein_num, ion):
    data_list = []
    for charge in [1,2]:
        mass = computational_pipeline.gen_spectra.get_max_mass(kmer, ion=ion, charge=charge)
        ion_int = 0 if ion == 'b' else 1
        input_tuple = (mass, start, end, ion_int, charge, protein_num)
        data_list.append(input_tuple)
    return data_list

def db_make_set_for_protein_digest(i,prot,max_len, dbf, data, digest):
    seq_len = len(prot)
    count_max = 1000000
    for size in range(2, max_len + 1):
        for start in range(0, seq_len - size + 1):
            end = start + size
            kmer = prot[start:end]
            if kmer[0] in digest[0] or digest[0] == ['-'] or (start > 0 and prot[start-1] in digest[1]):
                bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'b')
                    data.extend(data_list)
                    if len(data) > count_max:
                        dbf.insert(data)
                        data.clear()
            if kmer[-1] in digest[1] or digest[1] == ['-'] or (end < seq_len and prot[end] in digest[0]):
                bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
                if not any (x in bad_chars for x in kmer):
                    data_list = get_data(kmer, start, end, i, 'y')
                    data.extend(data_list)
                    if len(data) > count_max:
                        dbf.insert(data)
                        data.clear()
            
    return

def db_make_database_set_for_proteins(proteins,max_len,dbf,digest):
    plen = len(proteins)
    last_percent = 0
    data = []
    
    for i, (_, prot_entry) in enumerate(proteins):
        percent = int((i+1) * 100 / plen)
        print(f'\rOn protein {i+1}/{plen} [{int((i+1) * 100 / plen)}%]', end='')
        if percent != last_percent:
            last_percent = percent
            free = shutil.disk_usage('/')[2]
            free = free/(1024**3)
            if free < 10:
                print("\nUsed too much space, Space available =", free, "GB" )
                sys.exit(1)
        db_make_set_for_protein_digest(i,prot_entry,max_len, dbf, data, digest)
        
    if len(data) != 0:
        dbf.insert(data)
        
def modified_make_database_set(proteins: list, max_len: int, dbf, digest):
    db_make_database_set_for_proteins(proteins,max_len,dbf,digest)
    dbf.index_ion_mass()
    dbf.index_ion_mass_b()
    dbf.index_ion_mass_y()
    return


