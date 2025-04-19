import os
import click
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Dict, Union, Tuple
from loguru import logger


class BiologicalSequence(ABC):
    """
    Abstract class for biological sequences.

    Provides basic functionality for working with biological sequences, such as
    length, indexing, string representation, and a method to validate the
    sequence.
    Methods:
        __len__(): Returns the length of the sequence.
        __getitem__(index): Returns the sequence element at the given index.
        __str__(): Returns the sequence as a string.
        __repr__(): Returns a string representation of the sequence object.
        _is_bio_sequence(): Abstract method to check if the sequence is valid.
    """
    def __init__(self, sequence: str):
        self.sequence = sequence
        if not self._is_bio_sequence():
            raise ValueError("Invalid sequence")

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def _is_bio_sequence(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Class for nucleic acid sequences.

    Provides functionality specific to nucleic acid sequences,
    such as complement, reverse, reverse complement, GC content,
    palindrome checking, and nucleotide counting.

    Attributes:
        _complement_map (Dict[str, str]): A map for complementing nucleotides.
        _valid_monomers (set): A set of valid nucleotides for the sequence.
    Methods:
        complement(): Returns the complementary sequence.
        reverse(): Returns the reversed sequence.
        reverse_complement(): Returns the reverse complement of the sequence.
        gc_content(): Returns the GC content percentage.
        is_palindrome(): Checks if the sequence is a palindrome.
        nucleotide_count():
            Returns a dictionary with counts of each nucleotide.
        _is_bio_sequence():
            Validates if the sequence contains valid nucleotides.
    """
    _complement_map: Dict[str, str] = {}
    _valid_monomers: set = set()

    def complement(self):
        return self.__class__(self.sequence.translate(
            str.maketrans(self._complement_map)
        ))

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()

    def gc_content(self) -> float:
        gc_count = sum(1 for n in self.sequence if n.upper() in "GC")
        return round((gc_count / len(self.sequence)) * 100, 3) \
            if self.sequence \
            else 0.0

    def is_palindrome(self) -> bool:
        return self.sequence == self.reverse_complement().sequence

    def nucleotide_count(self) -> Dict[str, int]:
        nucleotide_counter = {}
        for nucleotide in self.sequence.upper():
            if nucleotide in nucleotide_counter:
                nucleotide_counter[nucleotide] += 1
            else:
                nucleotide_counter[nucleotide] = 1

        return nucleotide_counter

    def _is_bio_sequence(self) -> bool:
        return set(self.sequence.upper()).issubset(self._valid_monomers)


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences.

    Inherits from NucleicAcidSequence and provides
    functionality specific to DNA sequences: transcription to RNA.

    Attributes:
        _complement_map (Dict[str, str]):
            A map for complementing DNA nucleotides.
        _valid_monomers (set):
            A set of valid nucleotides for DNA sequences.
    Methods:
        transcribe(): Transcribes the DNA sequence to RNA.
    """
    _complement_map = {"A": "T", "T": "A", "C": "G", "G": "C",
                       "a": "t", "t": "a", "c": "g", "g": "c"}
    _valid_monomers = {"A", "T", "C", "G", "a", "t", "c", "g"}

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U").replace("t", "u"))


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences.

    Inherits from NucleicAcidSequence.

    Attributes:
        _complement_map (Dict[str, str]):
            A map for complementing RNA nucleotides.
        _valid_monomers (set):
            A set of valid nucleotides for RNA sequences.
    """
    _complement_map = {"A": "U", "U": "A", "C": "G", "G": "C",
                       "a": "u", "u": "a", "c": "g", "g": "c"}
    _valid_monomers = {"A", "U", "C", "G", "a", "u", "c", "g"}


class AminoAcidSequence(BiologicalSequence):
    """
    Class for amino acid sequences.

    Provides functionality for molecular weight calculation.

    Attributes:
        _valid_monomers (set): A set of valid amino acids.
    Methods:
        molecular_weight():
            Calculates the molecular weight of the amino acid sequence.
        _is_bio_sequence():
            Validates if the sequence contains only valid amino acids.
    """
    _valid_monomers = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def _is_bio_sequence(self) -> bool:
        return set(self.sequence).issubset(self._valid_monomers)

    def molecular_weight(self) -> float:
        weights = {
            "A": 89.1, "C": 121.2, "D": 133.1, "E": 147.1, "F": 165.2,
            "G": 75.1, "H": 155.2, "I": 131.2, "K": 146.2, "L": 131.2,
            "M": 149.2, "N": 132.1, "P": 115.1, "Q": 146.1, "R": 174.2,
            "S": 105.1, "T": 119.1, "V": 117.1, "W": 204.2, "Y": 181.2
        }
        weights.update({k.lower(): v for k, v in weights.items()})
        return round(sum(weights[aa] for aa in self.sequence), 2)


def filter_fastq(
        input_fastq: str,
        output_fastq: str,
        gc_bounds: Union[Tuple, int, float] = (0, 100),
        length_bounds: Union[Tuple[int, int], int] = (0, 2 ** 32),
        quality_threshold: Union[int, float] = 0
) -> None:
    """
    Filters reads from a FASTQ file based on GC content,
    sequence length, and quality score, and writes the filtered reads
    to a new FASTQ file in the 'filtered' directory.
    Creates "./logs" directory and performs logging into
    "./logs/filter_fastq.log" file.
    The function is available from the terminal.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        output_fastq (str): Name of the output filtered FASTQ file.
        gc_bounds (Union[Tuple, int, float]): GC content bounds.
        length_bounds (Union[Tuple[int, int], int]): Sequence length bounds.
        quality_threshold (Union[int, float]): Minimum read quality score.

    Returns: None
    """
    output_directory = 'filtered'
    output_path = os.path.join(output_directory, output_fastq)

    if not getattr(logger, "_configured", False):
        os.makedirs("logs", exist_ok=True)
        logger.add(
            "logs/filter_fastq.log",
            rotation="500 KB",
            level="INFO",
            format="""
==================== {level} ====================
Time:    <green>{time:YYYY-MM-DD HH:mm:ss}</green>
Module:  <cyan>{name}</cyan>:{function}:{line}
Message: {message}
            """
        )

    logger.info(
        f"""Starting FASTQ filtering:
        Input: {input_fastq}
        Output: {output_path}
        GC bounds: {gc_bounds}
        Length bounds: {length_bounds}
        Quality ≥ {quality_threshold}
        """
    )

    if os.path.exists(output_path):
        logger.error(f"Output file '{output_path}' already exists.")
        raise FileExistsError(f"File {output_fastq} already exists.")

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    filtered_count = 0
    total_count = 0

    with (open(output_path, 'w') as out_file):
        for fastq_seq in SeqIO.parse(input_fastq, "fastq"):
            total_count += 1
            seq_len = len(fastq_seq.seq)
            gc_content = gc_fraction(fastq_seq.seq) * 100
            sum_quality = sum(fastq_seq.letter_annotations["phred_quality"])
            len_quality = len(fastq_seq.letter_annotations["phred_quality"])
            quality = sum_quality / len_quality

            if not (length_bounds[0] <= seq_len <= length_bounds[1]):
                continue
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue
            if not (quality_threshold <= quality):
                continue

            SeqIO.write(fastq_seq, out_file, "fastq")
            filtered_count += 1

    logger.info(
        f"""Filtering finished.
        Total reads: {total_count} 
        Passed filters: {filtered_count}
        """
    )


@click.group()
def cli():
    """SeqSmith: bioinformatics tool"""
    pass


@cli.command(name="filter-fastq")
@click.argument('input_fastq')
@click.argument('output_fastq')
@click.option('-g', '--gc-bounds', nargs=2, type=float, default=(0, 100), help='GC content bounds (e.g., 40 60).')
@click.option('-l', '--length-bounds', nargs=2, type=int, default=(0, 2**32), help='Sequence length bounds.')
@click.option('-q', '--quality-threshold', type=float, default=0, help='Minimum average quality threshold.')
def filter_fastq_cmd(
        input_fastq,
        output_fastq,
        gc_bounds,
        length_bounds,
        quality_threshold
):
    filter_fastq(
        input_fastq=input_fastq,
        output_fastq=output_fastq,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold
    )
    click.echo(f"File saved to ./filtered/{output_fastq}")


if __name__ == '__main__':
    cli()
