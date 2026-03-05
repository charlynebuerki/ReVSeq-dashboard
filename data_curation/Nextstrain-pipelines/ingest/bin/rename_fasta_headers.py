import argparse
from Bio import SeqIO

def rename_fasta_headers(input_file, output_file, new_name):
    """
    Renames the headers in a FASTA file to the specified new_name.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id = new_name  # Update the header
            record.description = ""  # Remove any description
            SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename FASTA headers to match the file name")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--name", required=True, help="New name for the header")
    args = parser.parse_args()

    rename_fasta_headers(args.input, args.output, args.name)
