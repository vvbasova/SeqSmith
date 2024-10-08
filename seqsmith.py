from typing import Union, List, Dict, Tuple
import os

from utils.dna_rna_utils import (
    is_dna,
    is_rna,
    transcribe,
    reverse,
    reverse_complement,
    complement,
    gc_content,
    is_palindrome,
    nucleotide_count
)

from utils.filter_fastq_utils import (
    filter_gc,
    filter_length,
    filter_quality,
    fastq_to_dict_generator,
    convert_dict_to_fastq
)


def run_dna_rna_tools(
        *args: str
) -> Union[
    str,
    bool,
    Dict[str, int],
    List[Union[str, bool, Dict[str, int]]]
]:
    """
    Takes nucleic acid sequences and a procedure for transforming them.
    Checks whether the sequences and the procedure are valid
    and if not, raises a ValueError.
    The nucleic acid sequence may only contain the following characters:
    A, T, C, G, U, a, t, c, g, u
    The sequence cannot contain both T/t and U/u characters simultaneously.

    Args:
        *args: str
        The first arguments are one or more strings containing
        nucleic acid sequences.
        The last argument is a string containing the procedure:
        - transcribe: Transcribe the sequences
        - reverse: Reverse the sequences
        - complement: Build complementary sequences
        - reverse_complement: Build reversed complementary sequences
        - gc_content: Calculate the percentage of G and C content
        - is_palindrome: Check if the sequences are palindromes
        - nucleotide_count: Count the number of each nucleotide

    Returns:
        str | bool | dict | list of str | list of bool | list of dict
        The results of the procedure for one or more sequences.

    Raises:
        ValueError:
        If the sequence contains invalid characters or is an empty string.
    """
    procedure = args[-1]
    sequences = list(args[:-1])

    for sequence in sequences:
        if not (is_dna(sequence) or is_rna(sequence)):
            raise ValueError("Invalid nucleotide sequence")

    procedures = {
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement,
        'gc_content': gc_content,
        'is_palindrome': is_palindrome,
        'nucleotide_count': nucleotide_count
    }

    if procedure not in procedures:
        raise ValueError("Invalid procedure")

    result = procedures[procedure](sequences)

    return result[0] if len(result) == 1 else result


def filter_fastq(
        input_fastq: str,
        output_fastq: str,
        gc_bounds: Union[
            Tuple[Union[int, float], Union[int, float]],
            int,
            float
        ] = (0, 100),
        length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
        quality_threshold: Union[int, float] = 0
) -> None:
    """
    Filters reads from a fastq file based on GC content,
    sequence length, and quality score, and writes the filtered reads
    to a new fastq file in the 'filtered' directory.
    If the 'filtered' directory does not exist, creates it.

    Args:
        input_fastq (str):
            Path to the input fastq file.
        output_fastq (str):
            Name of the output filtered fastq file.
            The file will be created in the 'filtered' directory.
        gc_bounds (Union[Tuple, int, float]):
            GC content bounds. If a single value is provided, it is
            interpreted as the upper bound. Default is (0, 100).
        length_bounds (Union[Tuple[int, int], int], optional):
            Sequence length bounds.If a single value is provided, it is
            interpreted as the upper bound. Default is (0, 2**32).
        quality_threshold (Union[int, float]):
            Minimum read quality score. Default is 0.

    Returns: None

    Raises:
        FileExistsError:
            If the output file with the specified name already
            exists in the 'filtered' directory, the function raises
            an error to prevent accidental overwriting of existing data.
    """
    output_path = os.path.join('filtered', output_fastq)
    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_fastq} already exists.")

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    for seq_dict in fastq_to_dict_generator(input_fastq):
        if not filter_length(seq_dict, length_bounds):
            continue
        if not filter_gc(seq_dict, gc_bounds):
            continue
        if not filter_quality(seq_dict, quality_threshold):
            continue

        convert_dict_to_fastq(seq_dict, output_fastq)
