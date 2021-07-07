from utils import make_valid_fasta_file, file_exists

'''write

DESC:
    write a fasta file
Inputs:
    output_name: str name of file to write to 
    sequences: list of dictionaries of form {'name': str, 'sequence': str}
Outputs:
    name of the output file written to
'''
def write(output_name, sequences):
    output_name = make_valid_fasta_file(output_name)
    with open(output_name, 'w') as o:
        for i, seq in enumerate(sequences):
            o.write('>sp|{}|{}\n{}\n'.format('id{}'.format(i), seq['name'], seq['sequence']))
    return output_name

def read(fasta_file: str, is_uniprot=False) -> list:
    '''read
    Read proteins into memory from fasta file
    
    Inputs:
        fasta_file: str path to fasta file
    kwargs:
        is_uniprot: bool adds attribute 'human_readable_name' to dictionary if True. Default=False
    Outputs:
        list        list of dictionaries of form {'name': str, 'sequence': str, 'identifier': str}
    '''
    if not file_exists(fasta_file):
        raise Exception('File {} does not exist'.format(fasta_file))
    prots = []
    with open(fasta_file, 'r') as i:
        name = None 
        seq = '' 
        identifier = ''
        hmn_rdble_name = ''
        for line in i:
            if '>' in line: #name line

                # add the last thing to the list
                if not ((name is None or name == '') and (seq is None or seq == '')):
                    entry = {
                        'name': name,
                        'sequence': seq,
                        'identifier': identifier
                    }
                    if is_uniprot:
                        entry['human_readable_name'] = hmn_rdble_name
                    prots.append(entry)

                seq = '' 
                name = str(str(line.split('|')[2]).split(' ')[0]).replace('\n', '')
                identifier = str(line.split('|')[1])
                if is_uniprot:
                    after_bar = str(line.split('|')[2])
                    hmn_rdble_name = str(' '.join(after_bar.split(' ')[1:]).split('OS=')[0]).strip()
            else:
                seq += line.replace('\n', '')
        # add the last one
        entry = {
            'name': name,
            'sequence': seq,
            'identifier': identifier
        }
        if is_uniprot:
            entry['human_readable_name'] = hmn_rdble_name
        prots.append(entry)
    return prots