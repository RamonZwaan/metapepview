from typing import IO, Generator, Tuple, Any
from pathlib import Path

import pandas as pd
import numpy as np

from backend.utils import digest_proteins




def import_fasta(fasta_file: Path) -> pd.Series:
    """Import fasta file and return protein contents as Series.

    Args:
        fasta_file (Path): Path to fasta file.

    Returns:
        pd.Series: Series of protein sequences
    """
    # list to store protein elements
    proteome = dict()
    
    # parse fasta file and store protein elements into output list
    with fasta_file.open() as fasta_content:
        for prot_el in _fasta_parser(fasta_content):
            proteome.update({prot_el[0]: prot_el[1]})

    return pd.Series(proteome, name="Sequence")


def _fasta_parser(fasta_content: IO[str]) -> Generator[Tuple[str, str], None, None]:
    """Return iterator that parses fasta files and extracts information per item.

    Args:
        fasta_content (TextIOWrapper): Fasta file contents.
        header_parser (FastaHeaderParser): Function extracting data from fasta header.

    Yields:
        Generator[Tuple[str, str]]: Generator yielding fasta header and protein sequence.
    """
    # protein information
    header = ""
    seq = ""
    
    # iterate over fasta file lines and collect protein data
    while True:
        line = fasta_content.readline()
        
        # end of file
        if not line:
            if seq:
                yield (header, seq)
                seq = ""
            break
        
        # header line
        if line.startswith(">"):
            # add previous protein data to list and reset local vars
            if seq:
                yield (header, seq)
                seq = ""

            # extract header information
            header = line.lstrip(">")

        # sequence line
        else:
            seq += line.replace("\n", "")


def fasta_to_peptides(fasta_file: Path, *args, **kwargs) -> pd.Series:
    """Import fasta file, digest proteins, and return peptide contents as
    Series.

    Args:
        fasta_file (Path): Path to fasta file.

    Returns:
        pd.Series: Series of peptide sequences
    """
    protein_series = import_fasta(fasta_file)
    
    # digest cRAP proteins into tryptic peptides
    return protein_series.transform(
        lambda x: digest_proteins(x, *args, **kwargs)).explode()