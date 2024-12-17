#!/usr/bin/env python3

def parse_extra_columns(extra_columns_str):
    extra_columns = extra_columns_str.split(';')
    column_mapping = {}
    for column in extra_columns:
        key, value = column.split('::')
        if key not in column_mapping:
            column_mapping[key] = []
        column_mapping[key].append(value)
    return extra_columns, column_mapping


def parse_gff(gff_file, mRNA_Type='mRNA', extra_columns=''):
    gene_id_mapping = {}
    mrna_id_mapping = {}

    if extra_columns:
        extra_columns, extra_columns_mapping = parse_extra_columns(extra_columns)
    else:
        extra_columns = ''
        extra_columns_mapping = {}
    from commands.read_data import open_file
    with open_file(gff_file) as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            SeqID, Source, Type, Start, End, Score, Strand, Phase, Attributes = fields
            Start, End = int(Start), int(End)

            attr_dict = {}
            for attr in Attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=')
                    attr_dict[key] = value

            if Type == 'gene':
                gene_id = attr_dict.get('ID', '')
                gene_name = attr_dict.get('Name', gene_id)
                gene_id_mapping[gene_id] = {
                    'gene_name': gene_name,
                    'longest_mRNA': '',
                    'max_length': 0
                }

                if Type in extra_columns_mapping.keys():
                    for col in extra_columns_mapping[Type]:
                        gene_id_mapping[gene_id][f'Extra_gene_{col}'] = attr_dict.get(col, '')

            elif Type == mRNA_Type:
                rna_id = attr_dict.get('ID', '')
                # gene_id = attr_dict.get('Parent', '')
                gene_id = attr_dict.get('Parent', rna_id)
                transcript_name = attr_dict.get('Name', rna_id)
                mrna_id_mapping[rna_id] = {
                    'gene_id': gene_id,
                    'transcript_name': transcript_name,
                    'SeqID': SeqID,
                    'Start': Start,
                    'End': End,
                    'Strand': Strand,
                    'CDS_names': [],
                    'CDS_starts': [],
                    'CDS_ends': [],
                    'CDS_length': 0
                }

                if Type in extra_columns_mapping.keys():
                    for col in extra_columns_mapping[Type]:
                        mrna_id_mapping[rna_id][f'Extra_mRNA_{col}'] = attr_dict.get(col, '')

            elif Type == "CDS":
                parent_id = attr_dict.get('Parent', '')
                cds_name = attr_dict.get('ID', '')
                cds_length = End - Start + 1

                if parent_id in mrna_id_mapping:
                    mrna_id_mapping[parent_id]['CDS_names'].append(cds_name)
                    mrna_id_mapping[parent_id]['CDS_starts'].append(Start)
                    mrna_id_mapping[parent_id]['CDS_ends'].append(End)
                    mrna_id_mapping[parent_id]['CDS_length'] += cds_length

                if Type in extra_columns_mapping.keys():
                    for col in extra_columns_mapping[Type]:
                        mrna_id_mapping[rna_id][f'Extra_CDS_{col}'] = attr_dict.get(col, '')

    # Determine the longest mRNA for each gene
    for rna_id, rna_attr in mrna_id_mapping.items():
        gene_id = rna_attr['gene_id']
        if gene_id in gene_id_mapping:
            if rna_attr['CDS_length'] > gene_id_mapping[gene_id]['max_length']:
                gene_id_mapping[gene_id]['max_length'] = rna_attr['CDS_length']
                gene_id_mapping[gene_id]['longest_mRNA'] = rna_id
        else:
            gene_id_mapping[gene_id] = {
                'gene_name': rna_id,
                'longest_mRNA': '',
                'max_length': 0
            }
            gene_id_mapping[gene_id]['max_length'] = rna_attr['CDS_length']
            gene_id_mapping[gene_id]['longest_mRNA'] = rna_id
            print(f"Warning: mRNA {rna_id}'s gene parent not found in gene id line. Map it's gene id/name to itself ({rna_id}).")


    return gene_id_mapping, mrna_id_mapping, extra_columns


def write_idmapping_file(gene_id_mapping, mrna_id_mapping, extra_columns, output_file):
    with open(output_file, 'w') as f:
        header_columns = [
            'gene_id', 'gene_name', 'transcript_id', 'transcript_name',
            'SeqID', 'Start', 'End', 'Strand', 'CDS_names', 'CDS_starts', 'CDS_ends',
            'CDS_length', 'is_longest_mRNA'
        ]
        if extra_columns:
            for col in extra_columns:
                type, col = col.split("::")
                header_columns.append(f'Extra::{type}::{col}')
        f.write('\t'.join(header_columns) + '\n')

        for rna_id, rna_attr in mrna_id_mapping.items():
            try:
                gene_id = rna_attr['gene_id']
                is_longest = (gene_id_mapping[gene_id]['longest_mRNA'] == rna_id)
                output_data = [
                    gene_id,
                    gene_id_mapping[gene_id]['gene_name'],
                    rna_id,
                    rna_attr['transcript_name'],
                    rna_attr['SeqID'],
                    rna_attr['Start'],
                    rna_attr['End'],
                    rna_attr['Strand'],
                    ','.join(rna_attr['CDS_names']),
                    ','.join(map(str, rna_attr['CDS_starts'])),
                    ','.join(map(str, rna_attr['CDS_ends'])),
                    rna_attr['CDS_length'],
                    'Yes' if is_longest else 'No'
                ]
                if extra_columns:
                    for col in extra_columns:
                        type = col.split("::")[0]
                        col = col.split("::")[1]
                        if type == "gene":
                            output_data.append(gene_id_mapping[gene_id].get(f'Extra_gene_{col}', '-'))
                        elif type == "mRNA":
                            output_data.append(rna_attr.get(f'Extra_mRNA_{col}', '-'))
                        elif type == "CDS":
                            output_data.append(rna_attr.get(f'Extra_CDS_{col}', '-'))
                f.write('\t'.join(map(str, output_data)) + '\n')
            except KeyError:
                print(f"Error: mRNA {rna_id} update failed. Please check.")


def setup_parser(parser):
    idmap_parser = parser.add_parser('map', help='Convert GFF file to ID_MAP format')
    idmap_parser.add_argument('-g', '--gff_file', required=True, help='Path to gff file')
    idmap_parser.add_argument('-o', '--output_file', default='id_mapping.txt', help='Path to the output file. [id_mapping.txt]')
    idmap_parser.add_argument('-t', '--trans_mRNA_info_to', default='mRNA', help='Transcript or mRNA. [mRNA]')
    idmap_parser.add_argument('-e', '--extra_info', help='Extra information that you need, for example: -e "mRNA::Dbxref;gene::gbkey". [NULL]')
    return idmap_parser
    '''
    如果-e有除了gene和mRNA的情况, 需要同时改变-t参数, 如果改变了-t参数, 就不能再得到mRNA的其他信息
    '''

def run(args):
    gene_id_mapping, mrna_id_mapping, extra_columns = parse_gff(args.gff_file, mRNA_Type = args.trans_mRNA_info_to, extra_columns=args.extra_info)
    write_idmapping_file(gene_id_mapping, mrna_id_mapping, extra_columns, args.output_file)

