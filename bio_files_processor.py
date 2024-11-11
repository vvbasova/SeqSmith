from typing import Optional, List
import os


def convert_multiline_fasta_to_oneline(
        input_fasta: str,
        output_fasta: Optional[str] = None
) -> None:
    """
    Takes fasta file where nucleic acid sequences are on
    multiple lines and creates new fasta file where every
    sequence is on a single line.

    Args:
        input_fasta: str
            Name of an input fasta file.
        output_fasta: str
            Name of a new fasta file which will be created.

    Returns: None

    Raises:
        FileExistsError: If the output file with the specified name already
        exists, the function raises an error to prevent accidental
        overwriting of existing data.

    """
    if output_fasta is None:
        output_fasta = input_fasta.rsplit(".", 1)[0] + "_one_line.fasta"

    if os.path.exists(output_fasta):
        raise FileExistsError(f"The file '{output_fasta}' already exists.")

    with (open(input_fasta, 'r') as input_fasta_file,
          open(output_fasta, 'w') as output_fasta_file):
        sequence = ""

        for line in input_fasta_file:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    output_fasta_file.write(sequence + "\n")
                output_fasta_file.write(line + "\n")
                sequence = ""
            else:
                sequence += line

        if sequence:
            output_fasta_file.write(sequence + "\n")


def parse_blast_output(
        input_file: str,
        output_file: Optional[str] = None
) -> None:
    """
    Takes the first protein name from each query in the input
    BLAST file, sorts the result and creates a new file with
    list of proteins.

    Args:
        input_file: str
            Name of an input BLAST file.
        output_file: str
            Name of a new text file which will be created.

    Returns: None

    Raises:
        FileExistsError: If the output file with the specified name already
        exists, the function raises an error to prevent accidental
        overwriting of existing data.
    """
    if output_file is None:
        output_file = input_file.rsplit(".", 1)[0] + "_parse.txt"

    if os.path.exists(output_file):
        raise FileExistsError(f"The file '{output_file}' already exists.")

    descriptions = []
    with open(input_file, 'r') as blast_input_file:
        for line in blast_input_file:
            line = line.strip()
            if line.startswith("Description"):
                descr = next(blast_input_file).strip()
                mod_descr = descr.split("  ")[0].strip()
                descriptions.append(mod_descr)
    descriptions.sort()
    with open(output_file, 'w') as blast_output_file:
        for description in descriptions:
            blast_output_file.write(description + '\n')


def select_genes_from_gbk_to_fasta(
        input_gbk: str,
        genes: List[str],
        n_before: int = 1,
        n_after: int = 1,
        output_fasta: Optional[str] = None,
) -> None:
    """
    Takes gbk file, adds selected genes from genes list,
    genes before and after selected genes and its translations
    to the new fasta file.

    Args:
        input_gbk: str
            Name of an input gbk file.
        output_fasta: str
            Name of a new fasta file which will be created.
        genes: List[str]
            List with names of selected genes
        n_before:
            Number of genes before the selected gene which
            will be added to a new fasta file. Default: 1
        n_after
            Number of genes after the selected gene which
            will be added to a new fasta file. Default: 1

    Returns: None

    Raises:
        FileExistsError: If the output file with the specified name already
        exists, the function raises an error to prevent accidental
        overwriting of existing data.
    """
    if output_fasta is None:
        output_fasta = input_gbk.rsplit(".", 1)[0] + "_selected.fasta"

    if os.path.exists(output_fasta):
        raise FileExistsError(f"The file '{output_fasta}' already exists.")

    fasta_list = []
    gene = None
    translation = None
    translation_is_in_gene = False

    with open(input_gbk, 'r') as input_gbk_file:
        for line in input_gbk_file:
            line = line.strip()

            if line.startswith("/gene="):

                if gene and translation:
                    fasta_list.append((gene, translation))

                gene = line.split("=")[-1].strip('"')
                translation = None
                translation_is_in_gene = True

            elif line.startswith("/translation=") and translation_is_in_gene:
                translation = line.split("=")[-1].strip().strip('"')

                while not line.endswith('"'):
                    line = next(input_gbk_file).strip()
                    translation += line.strip('"')

                translation_is_in_gene = False

        if gene and translation:
            fasta_list.append((gene, translation))

    genes_index = []
    for index, (gene, translation) in enumerate(fasta_list):
        if gene in genes:
            genes_index.append(index)

    with open(output_fasta, 'w') as output_fasta_file:
        for index in genes_index:
            start = max(0, index - n_before)
            end = min(len(fasta_list), index + n_after + 1)

            for gene, translation in fasta_list[start:end]:
                output_fasta_file.write(f">{gene}\n")
                output_fasta_file.write(f"{translation}\n")
