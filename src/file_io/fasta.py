from src.lookups.utils import file_exists, make_valid_fasta_file


def write(output_name, sequences):
    output_name = make_valid_fasta_file(output_name)
    with open(output_name, "w") as o:
        for i, seq in enumerate(sequences):
            o.write(
                ">sp|{}|{}\n{}\n".format("id{}".format(i), seq["name"], seq["sequence"])
            )
    return output_name


def read(fasta_file: str, is_uniprot=False) -> list:
    if not file_exists(fasta_file):
        raise Exception("File {} does not exist".format(fasta_file))
    prots = []
    with open(fasta_file, "r") as i:
        name = None
        seq = ""
        identifier = ""
        hmn_rdble_name = ""
        for line in i:
            if ">" in line:
                if not ((name is None or name == "") and (seq is None or seq == "")):
                    entry = {"name": name, "sequence": seq, "identifier": identifier}
                    if is_uniprot:
                        entry["human_readable_name"] = hmn_rdble_name
                    prots.append(entry)

                seq = ""
                name = str(str(line.split("|")[2]).split(" ")[0]).replace("\n", "")
                identifier = str(line.split("|")[1])
                if is_uniprot:
                    after_bar = str(line.split("|")[2])
                    hmn_rdble_name = str(
                        " ".join(after_bar.split(" ")[1:]).split("OS=")[0]
                    ).strip()
            else:
                seq += line.replace("\n", "")
        entry = {"name": name, "sequence": seq, "identifier": identifier}
        if is_uniprot:
            entry["human_readable_name"] = hmn_rdble_name
        prots.append(entry)
    return prots
