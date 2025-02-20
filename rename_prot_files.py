#! /usr/bin/env python3
import argparse
import fileinput

def rename_sequences(input_file, new_pattern, inplace):
    with fileinput.FileInput(input_file, inplace=inplace, backup='.bak') as infile:
        for line in infile:
            if line.startswith('>'):
                header, rest = line.strip().split(' ', 1)
                new_header = f'{header.split("_")[0]}_{new_pattern}'
                print(new_header)
            else:
                print(line, end='')

def main():
    parser = argparse.ArgumentParser(description='Rename sequences in a FASTA file.')
    parser.add_argument('-i', '--input_file', required=True, help='Input FASTA file')
    parser.add_argument('-p', '--new_pattern', required=True, help='New pattern to append to sequence headers')
    parser.add_argument('-o', '--output_file', help='Output FASTA file with renamed sequences. If not provided, the input file will be modified in place.')
    
    args = parser.parse_args()

    if args.output_file:
        rename_sequences(args.input_file, args.new_pattern, False)
    else:
        rename_sequences(args.input_file, args.new_pattern, True)

if __name__ == '__main__':
    main()

print("New headers! \U0001F609 \U0001F60E \U0001F37A")

print("Done - Salud")