from Bio import SeqIO
import sys

def read_gene_transcript_mapping(mapping_file, match_column=3, map_column=1):
    gene_transcripts = {}
    try:
        from commands.read_data import open_file
        with open_file(mapping_file) as f:
            next(f)
            for line in f:
                name_list = line.strip().split('\t')
                if len(name_list) < max(match_column, map_column):
                    continue
                gene_name, match_name, map_name = name_list[0], name_list[match_column - 1], name_list[map_column - 1]
                if gene_name not in gene_transcripts:
                    gene_transcripts[gene_name] = []
                gene_transcripts[gene_name].append({'match_name': match_name, 'map_name': map_name})
    except FileNotFoundError:
        raise FileNotFoundError(f"Mapping file {mapping_file} not found.")
    except Exception as e:
        raise RuntimeError(f"Error reading mapping file: {e}")
    return gene_transcripts

def read_transcript_sequences(transcript_file):
    transcript_sequences = {}
    try:
        from commands.read_data import open_file
        with open_file(transcript_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                match_name = record.id
                sequence = str(record.seq)
                transcript_sequences[match_name] = sequence
    except FileNotFoundError:
        raise FileNotFoundError(f"Transcript file {transcript_file} not found.")
    except Exception as e:
        raise RuntimeError(f"Error reading transcript file: {e}")
    return transcript_sequences

def find_longest_transcripts(gene_transcripts, transcript_sequences, do_not_extract_longest=False):
    """
    筛选最长转录本或仅改名转录本序列。
    
    Parameters:
    - gene_transcripts: dict, 基因和转录本的对应关系
    - transcript_sequences: dict, 转录本序列
    - do_not_extract_longest: bool, 是否仅改名（默认 False）

    Returns:
    - dict: map_name 和序列的对应关系
    """
    renamed_transcripts = {}
    count_need_alter = 0
    count_altered = 0

    for gene_name, transcripts in gene_transcripts.items():
        if do_not_extract_longest:
            # 仅改名
            for transcript in transcripts:
                count_need_alter += 1
                match_name = transcript['match_name']
                map_name = transcript['map_name']
                if match_name in transcript_sequences:
                    if map_name in renamed_transcripts:
                        print(f"Warning: Sequence {map_name} already exists.", file=sys.stdout)
                        continue
                    else:
                        renamed_transcripts[map_name] = transcript_sequences[match_name]
                        count_altered += 1
                else:
                    print(f"ERROR: Sequence {match_name} does not exist in the id_mapping file.", file=sys.stderr)
        else:
            # 筛选最长转录本
            longest_transcript = None
            longest_length = 0
            for transcript in transcripts:
                count_need_alter += 1
                match_name = transcript['match_name']
                map_name = transcript['map_name']
                if match_name not in transcript_sequences:
                    continue
                seq_length = len(transcript_sequences[match_name])
                if seq_length > longest_length:
                    longest_transcript = map_name
                    longest_length = seq_length
            if longest_transcript:
                if longest_transcript in renamed_transcripts:
                    print(f"Warning: Sequence {map_name} already exists.", file=sys.stdout)
                    continue
                else:
                    renamed_transcripts[longest_transcript] = transcript_sequences[transcripts[0]['match_name']]
                    count_altered += 1
    print(f"The number of sequences inputed: {count_need_alter}.", file=sys.stdout)
    print(f"The number of sequences generated: {count_altered}.", file=sys.stdout)

    return renamed_transcripts


def write_longest_transcript_sequences(longest_transcripts, output_file):
    with open(output_file, 'w') as f:
        for map_name, sequence in longest_transcripts.items():
            f.write(f'>{map_name}\n{sequence}\n')

def setup_parser(parser):
    alter_parser = parser.add_parser(
        'alter', 
        help='Get longest isoform or rename isoforms according to ID_MAP and FASTA (cds/pep) files'
    )
    alter_parser.add_argument('-i', '--idmapping_file', required=True, help='Path to the ID_MAP file')
    alter_parser.add_argument('-f', '--fasta_file', required=True, help='Path to the transcript sequences file')
    alter_parser.add_argument('-o', '--output_file', default='output.fa', help='Path to the output file. [output.fa]')
    alter_parser.add_argument('-m', '--match_column',  type=int, default=3, help='Column for match_name in mapping file. [3]')
    alter_parser.add_argument('-n', '--map_column',  type=int, default=1, help='Column for map_name in mapping file. [1]')
    alter_parser.add_argument('-d', '--do_not_extract_longest', action='store_true', help='Do not extract longest transcript, just rename. [False]')
    return alter_parser


def run(args):
    """
    主函数，根据参数执行逻辑：
    - 筛选最长转录本并改名
    - 或仅改名
    """
    gene_transcripts = read_gene_transcript_mapping(
        args.idmapping_file, 
        match_column=args.match_column, 
        map_column=args.map_column
    )
    transcript_sequences = read_transcript_sequences(args.fasta_file)
    
    longest_transcripts = find_longest_transcripts(
        gene_transcripts, 
        transcript_sequences, 
        do_not_extract_longest=args.do_not_extract_longest
    )
    
    write_longest_transcript_sequences(longest_transcripts, args.output_file)