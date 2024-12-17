# gf_mapper

## 1. Description
GFF and FASTA file id mapping. Get id mapping information (ID_MAP) from GFF files. And change FASTA files according to the ID_MAP file. This tool is inspired by [zhangrengang](https://github.com/zhangrengang), but an extension.

## 2. Installation
```
git clone https://github.com/shwzhao/gf_mapper
export PATH=$PATH:$PWD/gf_mapper/bin
```

dependences:
- biopython
- argparse

## 3. Usage Example

```
gf_mapper -h
usage: df_mapper [-h] {gff2idmap,alter} ...

Comparative genomics analysis toolkit

options:
  -h, --help         show this help message and exit

subcommands:
  {gff2idmap,alter}
    gff2idmap        Convert GFF file to ID_MAP format
    alter            Get longest isoform or rename isoforms according to ID_MAP and FASTA (cds/pep) files
```


### 3.1 download example data
```
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_genomic.gff.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_protein.faa.gz
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna.gz

gunzip *gz

grep -c ">" GCF_034140825.1_ASM3414082v1_protein.faa GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna
## GCF_034140825.1_ASM3414082v1_protein.faa:42933
## GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna:42933

grep -v "#" GCF_034140825.1_ASM3414082v1_genomic.gff | cut -f 3 | sort | uniq -c
242528 CDS
1637 cDNA_match
341091 exon
37483 gene
   1 intron
   2 inverted_repeat
6551 lnc_RNA
42933 mRNA
   1 match
1910 pseudogene
   3 pseudogenic_tRNA
1756 rRNA
  15 region
 489 sequence_feature
  71 snRNA
 611 snoRNA
 708 tRNA
3927 transcript
```

### 3.2 id mapping
```
gf_mapper map -h                                            
usage: df_mapper map [-h] -g GFF_FILE [-o OUTPUT_FILE] [-t TRANS_MRNA_INFO_TO] [-e EXTRA_INFO]

options:
  -h, --help            show this help message and exit
  -g GFF_FILE, --gff_file GFF_FILE
                        Path to gff file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output file. [id_mapping.txt]
  -t TRANS_MRNA_INFO_TO, --trans_mRNA_info_to TRANS_MRNA_INFO_TO
                        Transcript or mRNA. [mRNA]
  -e EXTRA_INFO, --extra_info EXTRA_INFO
                        Extra information that you need, for example: -e "mRNA::Dbxref;gene::gbkey". [NULL]
```

```
gf_mapper map \
  -g GCF_034140825.1_ASM3414082v1_genomic.gff
## Warning: rna-NC_001751.1:1732..2019 not found in gene_id_mapping. Skipping update.
## Warning: rna-NC_001751.1:1793..2035 not found in gene_id_mapping. Skipping update.

head -5 id_mapping.txt | column -t
## gene_id          gene_name   transcript_id       transcript_name  SeqID        Start  End    Strand  CDS_names                                                                                                                                                                                      CDS_starts                                                CDS_ends                                                  CDS_length  is_longest_mRNA  Extra::gene::Dbxref  Extra::CDS::Name
## gene-LOC4326813  LOC4326813  rna-XM_066304510.1  XM_066304510.1   NC_089035.1  7328   15219  +       cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1,cds-XP_066160607.1                     7851,8759,9859,11538,12430,12634,13612,14506,14676        8018,8857,9962,12346,12552,13010,14019,14589,14699        2196        Yes              GeneID:4326813       XP_066160607.1
## gene-LOC4326813  LOC4326813  rna-XM_015766610.3  XM_015766610.3   NC_089035.1  7328   15219  +       cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1,cds-XP_015622096.1  7851,8759,9859,11538,12430,12634,12810,13612,14506,14676  8018,8857,9962,12346,12552,12722,13010,14019,14589,14699  2109        No               GeneID:4326813       XP_015622096.1
## gene-LOC4326813  LOC4326813  rna-XM_066304523.1  XM_066304523.1   NC_089035.1  7328   15219  +       cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1,cds-XP_066160620.1  7851,8759,9859,11538,12430,12685,12810,13612,14506,14676  8018,8857,9962,12346,12552,12722,13010,14019,14589,14699  2058        No               GeneID:4326813       XP_066160620.1
## gene-LOC4326813  LOC4326813  rna-XM_066304514.1  XM_066304514.1   NC_089035.1  7328   15219  +       cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1,cds-XP_066160611.1                     7851,8759,9859,11538,12430,12634,13612,14506,14676        8018,8857,9962,12346,12552,13010,14019,14589,14699        2196        No               GeneID:4326813       XP_066160611.1
```

### 3.3 name changing

```
gf_mapper alter -h
## usage: df_mapper alter [-h] -i IDMAPPING_FILE -f FASTA_FILE [-o OUTPUT_FILE] [-m MATCH_COLUMN] [-n MAP_COLUMN] [-d]
## 
## options:
##   -h, --help            show this help message and exit
##   -i IDMAPPING_FILE, --idmapping_file IDMAPPING_FILE
##                         Path to the ID_MAP file
##   -f FASTA_FILE, --fasta_file FASTA_FILE
##                         Path to the transcript sequences file
##   -o OUTPUT_FILE, --output_file OUTPUT_FILE
##                         Path to the output file. [output.fa]
##   -m MATCH_COLUMN, --match_column MATCH_COLUMN
##                         Column for match_name in mapping file. [3]
##   -n MAP_COLUMN, --map_column MAP_COLUMN
##                         Column for map_name in mapping file. [1]
##   -d, --do_not_extract_longest
##                         Do not extract longest transcript, just rename. [False]
```


- Change CDS name to gene ID. And get the longest isoforms.
```
gf_mapper alter \
  -i id_mapping.txt -f GCF_034140825.1_ASM3414082v1_protein.faa \
  -m 15 \
  -n 1
## The number of sequences inputed: 42933.
## The number of sequences generated: 29429.

grep ">" output.fa | head -5
## >gene-LOC4326813
## >gene-LOC112163590
## >gene-LOC4326455
## >gene-LOC112163591
## >gene-LOC4326456
```

- Change CDS name to mRNA ID. And get the longest isoforms.
```
gf_mapper alter -i id_mapping.txt -f GCF_034140825.1_ASM3414082v1_protein.faa -m 15 -n 3
## The number of sequences inputed: 42933.
## The number of sequences generated: 29429.

grep ">" output.fa | head -5
## >rna-XM_066304510.1
## >rna-NM_001361325.1
## >rna-NM_001361324.1
## >rna-XM_026027108.2
## >rna-XM_015765670.3
```

- Change CDS name to mRNA ID.
```
gf_mapper alter \
  -i id_mapping.txt \
  -f GCF_034140825.1_ASM3414082v1_protein.faa \
  -m 15 \
  -n 3 \
  -d
The number of sequences inputed: 42933.
The number of sequences generated: 42933.

grep ">" output.fa | head -5                                                                     
## >rna-XM_066304510.1
## >rna-XM_015766610.3
## >rna-XM_066304523.1
## >rna-XM_066304514.1
## >rna-XM_066304520.1
```
