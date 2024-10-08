from typing import Union, Dict, Tuple
import os


def filter_length(seqs: Dict[str, Tuple[str, str]],
                  length_bounds: Tuple[int, int]
                  ) -> Dict[str, Tuple[str, str]]:
    filtered_length_seqs = {}

    for name, (sequence, quality) in seqs.items():
        if length_bounds[0] <= len(sequence) <= length_bounds[1]:
            filtered_length_seqs[name] = (sequence, quality)

    return filtered_length_seqs


def gc_counter(sequence: str) -> Union[int, float]:
    g_content = sequence.upper().count('G')
    c_content = sequence.upper().count('C')
    sequence_length = len(sequence)
    gc_percentage = (g_content + c_content) / sequence_length * 100

    return gc_percentage


def filter_gc(seqs: Dict[str, Tuple[str, str]],
              gc_bounds: Tuple[Union[int, float], Union[int, float]]
              ) -> Dict[str, Tuple[str, str]]:
    filtered_gc_seqs = {}

    for name, (sequence, quality) in seqs.items():
        if gc_bounds[0] <= gc_counter(sequence) <= gc_bounds[1]:
            filtered_gc_seqs[name] = (sequence, quality)

    return filtered_gc_seqs


def quality_counter(quality: str) -> Union[int, float]:
    q_score = [ord(symbol) - 33 for symbol in quality]
    average_q_score = sum(q_score) / len(q_score)

    return average_q_score


def filter_quality(seqs: Dict[str, Tuple[str, str]],
                   quality_threshold: Union[int, float]
                   ) -> Dict[str, Tuple[str, str]]:
    filtered_quality_seqs = {}

    for name, (sequence, quality) in seqs.items():
        if quality_threshold <= quality_counter(quality):
            filtered_quality_seqs[name] = (sequence, quality)

    return filtered_quality_seqs


def fastq_to_dict_generator(input_fastq: str) -> None:
    with open(input_fastq, 'r') as fastq_file:
        while True:
            try:
                name = next(fastq_file).strip()
                sequence = next(fastq_file).strip()
                next(fastq_file)
                quality = next(fastq_file).strip()

                yield {name: (sequence, quality)}

            except StopIteration:
                break


def convert_dict_to_fastq(
        seq_dict: Dict[str, Tuple[str, str]],
        output_fastq: str
) -> None:
    output_directory = 'filtered'

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_path = os.path.join(output_directory, output_fastq)

    with open(output_path, 'a') as fastq_file:
        for name, (sequence, quality) in seq_dict.items():
            fastq_file.write(f"{name}\n")
            fastq_file.write(f"{sequence}\n")
            fastq_file.write(f"+{name[1:]}\n")
            fastq_file.write(f"{quality}\n")
