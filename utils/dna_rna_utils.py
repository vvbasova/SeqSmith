from typing import Union, List, Dict


def is_dna(sequence: str) -> bool:
    if not sequence:
        return False
    return set(sequence).issubset({'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'})


def is_rna(sequence: str) -> bool:
    if not sequence:
        return False
    return set(sequence).issubset({'A', 'U', 'C', 'G', 'a', 'u', 'c', 'g'})


def transcribe(sequences: List[str]) -> List[str]:
    transcribed_sequences = []

    for sequence in sequences:
        if is_dna(sequence):
            transcribed_sequence = sequence.replace('T', 'U').replace('t', 'u')
        elif is_rna(sequence):
            transcribed_sequence = sequence.replace('U', 'T').replace('u', 't')
        transcribed_sequences.append(transcribed_sequence)

    return transcribed_sequences


def reverse(sequences: List[str]) -> List[str]:
    reversed_sequences = []

    for sequence in sequences:
        reversed_sequence = sequence[::-1]
        reversed_sequences.append(reversed_sequence)

    return reversed_sequences


def complement(sequences: List[str]) -> List[str]:
    complement_sequences = []

    for sequence in sequences:
        if is_dna(sequence):
            trans_table = sequence.maketrans("ATCGatcg", "TAGCtagc")
            complement_sequence = sequence.translate(trans_table)
            complement_sequences.append(complement_sequence)
        elif is_rna(sequence):
            trans_table = sequence.maketrans("AUCGaucg", "UAGCuagc")
            complement_sequence = sequence.translate(trans_table)
            complement_sequences.append(complement_sequence)

    return complement_sequences


def reverse_complement(sequences: List[str]) -> List[str]:
    return reverse(complement(sequences))


def gc_content(sequences: List[str]) -> List[Union[int, float]]:
    gc_contents = []

    for sequence in sequences:
        g_content = sequence.upper().count('G')
        c_content = sequence.upper().count('C')
        sequence_length = len(sequence)
        gc_percentage = (g_content + c_content) / sequence_length * 100
        gc_contents.append(round(gc_percentage, 3))

    return gc_contents


def is_palindrome(sequences: List[str]) -> List[bool]:
    complement_seq = reverse_complement(sequences)
    palindromes = []

    for sequences, complement_seq in zip(sequences, complement_seq):
        palindromes.append(sequences.upper() == complement_seq.upper())

    return palindromes


def nucleotide_count(sequences: List[str]) -> List[Dict[str, int]]:
    nucleotide_counts = []

    for sequence in sequences:
        nucleotide_counter = {}
        for nucleotide in sequence.upper():
            if nucleotide in nucleotide_counter:
                nucleotide_counter[nucleotide] += 1
            else:
                nucleotide_counter[nucleotide] = 1
        nucleotide_counts.append(nucleotide_counter)

    return nucleotide_counts
