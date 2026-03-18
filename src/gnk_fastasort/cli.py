import sys
import logging
import pathlib
import argparse

from gnk_fastasort.sanger_file_organiser import main as sanger
from gnk_fastasort.ncbi_file_organiser import main as ncbi

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("gnk_fastasort.log"),  # logs to file
        logging.StreamHandler()  # logs to console
    ]
)
logger = logging.getLogger('gnk_fastasort_logger')

def parse_args():
    parser = argparse.ArgumentParser(description="Sort FASTA files in ")
    parser.add_argument("--index", help="Index file", type=pathlib.Path)
    parser.add_argument("--gca_accession", help="Genome assembly accession number")
    parser.add_argument("--style", choices=["names", "full"], default="names",
        help="""Output style of data:
        'names' is a \\n seperated list of scaffold names,
        'full' is a tsv of  'names', 'original name', 'parent molecule' and 'length'.""")
    parser.add_argument("-o", "--output", help="Output file prefix", type=str)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s: 0.1.2")
    return parser.parse_args()

def main():
    args = parse_args()

    if args.index and not args.gca_accession:
        sanger(args)

    if args.gca_accession:
        ncbi(args)

    if not args.gca_accession and not args.index:
        sys.exit("MUST have at least one of --index (local file) or --gca_accession (ncbi hosted)")
