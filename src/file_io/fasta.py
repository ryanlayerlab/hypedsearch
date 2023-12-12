from lookups.utils import make_valid_fasta_file, file_exists

def write(output_name, sequences):
    output_name = make_valid_fasta_file(output_name)
    with open(output_name, 'w') as o:
        for i, seq in enumerate(sequences):
            o.write('>sp|{}|{}\n{}\n'.format('id{}'.format(i), seq['name'], seq['sequence']))
    return output_name

def read(fasta_filename: str, is_uniprot=False) -> list:
    if not file_exists(fasta_filename):
        raise Exception('File {} does not exist'.format(fasta_filename))
    prots = []
    with open(fasta_filename, 'r') as fasta_file:
        name = None 
        seq = '' 
        identifier = ''
        human_readable_name = ''
        for line in fasta_file:
            if '>' in line:
                if not ((name is None or name == '') and (seq is None or seq == '')):
                    entry = {
                        'name': name,
                        'sequence': seq,
                        'identifier': identifier
                    }
                    if is_uniprot:
                        entry['human_readable_name'] = human_readable_name
                    prots.append(entry)

                seq = '' 
                name = str(str(line.split('|')[2]).split(' ')[0]).replace('\n', '')
                identifier = str(line.split('|')[1])
                if is_uniprot:
                    after_bar = str(line.split('|', 3)[2])
                    human_readable_name = str(' '.join(after_bar.split(' ')[1:]).split('OS=')[0]).strip()
            else:
                seq += line.replace('\n', '')
        entry = {
            'name': name,
            'sequence': seq,
            'identifier': identifier
        }
        if is_uniprot:
            entry['human_readable_name'] = human_readable_name
        prots.append(entry)
    return prots