from Bio.Seq import Seq

def read_sequence_file(self, filepath, is_rc=False):
    '''
    '''
    filetype = filepath.split('.')[-1]
    match filetype:
        case "fasta":
            dict_seqs, max_seq_len = read_fasta(filepath, is_rc)
        case "csv":
            dict_seqs, max_seq_len = read_csv(filepath, is_rc)
        case "fna": # FASTA Nucleic Acids
            dict_seqs, max_seq_len = read_fna(filepath, is_rc)
        case "ffn": # FASTA Nucleotides of Gene Regions
            dict_seqs, max_seq_len = read_ffn(filepath, is_rc)
        case "faa": # FASTA Amino Acids
            dict_seqs, max_seq_len = read_faa(filepath, is_rc)
        case _:
            print("Filetype not supported")
    return dict_seqs, max_seq_len

def read_fasta(fasta_filepath, is_rc):
    fasta_reader = open(fasta_filepath, 'r')
    dict_seqs = {}
    max_seq_len = 0
    for line in fasta_reader:
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            print(seq_id)
        else:
            seq = line.strip('"').strip()
            if seq == '':
                continue
            else:
                if is_rc:
                    seq = str(Seq(seq).reverse_complement())
            print('sequence length = ' + str(len(seq)))
            if max_seq_len < len(seq):
                max_seq_len = len(seq)
            print(seq)
            dict_seqs[seq_id] = seq
    fasta_reader.close()
    return dict_seqs, max_seq_len

def read_csv(csv_filepath, is_rc):
    csv_reader = open(csv_filepath, 'r')
    dict_seqs = {}
    max_seq_len = 0
    for line in csv_reader:
        seq_id, seq = line.strip().split(',')
        if is_rc:
            seq = str(Seq(seq).reverse_complement())
        print('sequence length = ' + str(len(seq)))
        if max_seq_len < len(seq):
            max_seq_len = len(seq)
        print(seq)
        dict_seqs[seq_id] = seq
    csv_reader.close()
    return dict_seqs, max_seq_len

def read_fna(fna_filepath, is_rc):
    fna_reader = open(fna_filepath, 'r')
    dict_seqs = {}
    max_seq_len = 0
    for line in fna_reader:
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            print(seq_id)
        else:
            seq = line.strip('"').strip()
            if seq == '':
                continue
            else:
                if is_rc:
                    seq = str(Seq(seq).reverse_complement())
            print('sequence length = ' + str(len(seq)))
            if max_seq_len < len(seq):
                max_seq_len = len(seq)
            print(seq)
            dict_seqs[seq_id] = seq
    fna_reader.close()
    return dict_seqs, max_seq_len

def read_ffn(ffn_filepath, is_rc):
    ffn_reader = open(ffn_filepath, 'r')
    dict_seqs = {}
    max_seq_len = 0
    for line in ffn_reader:
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            print(seq_id)
        else:
            seq = line.strip('"').strip()
            if seq == '':
                continue
            else:
                if is_rc:
                    seq = str(Seq(seq).reverse_complement())
            print('sequence length = ' + str(len(seq)))
            if max_seq_len < len(seq):
                max_seq_len = len(seq)
            print(seq)
            dict_seqs[seq_id] = seq
    ffn_reader.close()
    return dict_seqs, max_seq_len

def read_faa(faa_filepath, is_rc):
    faa_reader = open(faa_filepath, 'r')
    dict_seqs = {}
    max_seq_len = 0
    for line in faa_reader:
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            print(seq_id)
        else:
            seq = line.strip('"').strip()
            if seq == '':
                continue
            else:
                if is_rc:
                    seq = str(Seq(seq).reverse_complement())
            print('sequence length = ' + str(len(seq)))
            if max_seq_len < len(seq):
                max_seq_len = len(seq)
            print(seq)
            dict_seqs[seq_id] = seq
    faa_reader.close()
    return dict_seqs, max_seq_len
