#!/usr/bin/env python3
import sys
import gzip
import argparse
import os

def complement_seq(seq: str) -> str:
    # DNA complement table handling standard IUPAC codes and case
    trans = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRWSMKHVDBNtgcaayrwsmkhvdbn")
    return seq.translate(trans)[::-1]

def open_file(path: str):
    if path == "-" or not path:
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def process_stream(infile, outfile):
    first_char = None
    buffer = []
    
    # Peek at the first non-whitespace character to determine type (> for FASTA, @ for FASTQ)
    while True:
        line = infile.readline()
        if not line:
            break
        stripped = line.strip()
        if not stripped:
            continue
        first_char = stripped[0]
        buffer.append(line)
        break
        
    if first_char == '>':
        # FASTA format: >header followed by sequence lines
        header = buffer[0]
        outfile.write(header)
        seq_lines = []
        for line in infile:
            if line.startswith('>'):
                # Process previous record
                full_seq = "".join(seq_lines).replace("\n", "").replace("\r", "")
                rc_seq = complement_seq(full_seq)
                outfile.write(rc_seq + "\n")
                seq_lines = []
                outfile.write(line)
            else:
                seq_lines.append(line)
        if seq_lines:
            full_seq = "".join(seq_lines).replace("\n", "").replace("\r", "")
            rc_seq = complement_seq(full_seq)
            outfile.write(rc_seq + "\n")
            
    elif first_char == '@':
        # FASTQ format: @header, sequence, +, quality
        while True:
            # We already have @header in buffer[0] (or read next)
            header = buffer[0] if buffer else infile.readline()
            if not header:
                break
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()
            if not qual:
                break
            
            rc_seq = complement_seq(seq.strip())
            rc_qual = qual.strip()[::-1] + "\n"
            
            outfile.write(header)
            outfile.write(rc_seq + "\n")
            outfile.write(plus)
            outfile.write(rc_qual)
            buffer = []
    else:
        print(f"Error: Unable to detect file type. First character was '{first_char}' (expected '>' for FASTA or '@' for FASTQ).", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Reverse complement FASTA or FASTQ files with automatic format detection.")
    parser.add_argument("-i", "--input", default="-", help="Input FASTA/FASTQ file (plain or .gz), default stdin")
    parser.add_argument("-o", "--output", default="-", help="Output file, default stdout")
    args = parser.parse_args()

    try:
        with open_file(args.input) as infile:
            if args.output == "-" or not args.output:
                process_stream(infile, sys.stdout)
            else:
                mode = "wt"
                with gzip.open(args.output, mode) if args.output.endswith(".gz") else open(args.output, mode) as outfile:
                    process_stream(infile, outfile)
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(0)

if __name__ == "__main__":
    main()
