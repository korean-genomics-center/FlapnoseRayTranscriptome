def read_fasta_as_dict(path_fasta):
    dict_fasta_to_lines = dict()

    header = None
    with open(path_fasta, 'r') as fr:
        for line in fr:
            if line.startswith('>'):
                header = line.strip('\n')[1:]
                assert dict_fasta_to_lines.get(header) == None
                dict_fasta_to_lines[header] = list()
            else:
                dict_fasta_to_lines[header].append(line.strip('\n'))
    dict_fasta = dict()
    for name, lines in dict_fasta_to_lines.items():
        dict_fasta[name] = ''.join(lines)

    return dict_fasta

def save_dict_fasta_into_file(dict_fasta, path_save):
    import os
    os.makedirs(os.path.dirname(path_save), exist_ok = True)
    with open(path_save, 'w') as fw:
        for name in dict_fasta.keys():
            fw.write(f">{name}\n")
            fw.write(dict_fasta[name] + '\n')
            