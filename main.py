#!/usr/bin/env python3

import argparse
from commands import map, alter

def main():
    parser = argparse.ArgumentParser(prog='df_mapper', description='GFF and FASTA file id mapping.')
    subparsers = parser.add_subparsers(title='subcommands', dest='subcommand')

    map.setup_parser(subparsers)
    alter.setup_parser(subparsers)

    args = parser.parse_args()

    if args.subcommand:
        if args.subcommand == 'map':
            map.run(args)
        elif args.subcommand == 'alter':
            alter.run(args)
        else:
            parser.print_help()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
