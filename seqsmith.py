from typing import Union, List, Dict, Tuple

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
    filter_quality
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
    Функция принимает последовательности нуклеиновых кислот и
    процедуру для их преобразования.
    Функция выполняет проверку, являются ли последовательности и
    процедура корректными, и, если нет, поднимает ValueError.
    Последовательность нуклеиновых кислот может содержать только следующие символы:
    A, T, C, G, U, a, t, c, g, u
    Последовательность не может одновременно содержать символы T/t и U/u

    Args:
        *args: str
        Первыми предеаются строка/строки, содержащие последовательность нуклеиновых кислот.
        Последней передаётся строка, содержащая процедуру:
        - transcribe: Транскрибировать последовательности
        - reverse: Перевернуть последовательности
        - complement: Построить комплементарные последовательности
        - reverse_complement: Построить перевернутые комплементарные последовательности
        - gc_content: Вычислить процент содержания G и C
        - is_palindrome: Проверить, являются ли последовательности палиндромами
        - nucleotide_count: Подсчитать количество каждого нуклеотида

    Returns:
        str | bool | dict | list of str | list of bool | list of dict
        Результаты процедуры для одной или нескольких последовательностей.

    Raises:
        ValueError: Если последовательность содержит некорректные символы
        либо является пустой строкой.
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
        seqs: Dict[str, Tuple[str, str]],
        gc_bounds: Union[Tuple[Union[int, float], Union[int, float]], int, float] = (0, 100),
        length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
        quality_threshold: Union[int, float] = 0
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрует последовательности формата FASTQ по длине, содержанию GC и качеству.
    Если функции подается одно число в качестве границ содержания GC или длины
    последовательности, оно принимается за верхнюю границу.

    Args:
        seqs: dict
        Словарь, где ключ — ID последовательности, значение —
        кортеж с последовательностью и её качеством.
        gc_bounds: tuple | int | float
        Диапазон допустимого содержания GC. По умолчанию (0, 100).
        length_bounds: tuple | int
        Диапазон допустимых длин последовательностей. По умолчанию (0, 2**32).
        quality_threshold: int | float
        Минимальный порог качества. По умолчанию 0.

    Returns:
        dict
        Отфильтрованный словарь удовлетворяющих условиям последовательностей.

    """

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    filtered_seqs = filter_length(seqs, length_bounds)

    if gc_bounds != (0, 100):
        filtered_seqs = filter_gc(filtered_seqs, gc_bounds)
    if quality_threshold > 0:
        filtered_seqs = filter_quality(filtered_seqs, quality_threshold)

    return filtered_seqs
