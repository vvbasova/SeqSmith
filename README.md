

# SeqSmith


---

**SeqSmith** is a simple program for working with nucleic acid sequences, filtering data in FASTQ format and processing bioinformatics files.


## Table of Contents

---
- [Main Features](#main-features)
- [How to install](#how-to-install)
- [How to use](#how-to-use)
- [System requirements](#system-requirements)
- [Project structure](#project-structure)

## Main Features

___

1. **run_dna_rna_tools**: proceeds various operations with DNA or RNA sequences.

    - Supported operations: 

       - `transcribe`: returns the result of transcribing DNA to RNA or reverse transcribing RNA to DNA.
       - `reverse`: returns the reversed sequence of DNA or RNA.
       - `complement`: returns a DNA sequence complementary to the input DNA, or an RNA sequence complementary to the input RNA.
       - `reverse complement`: returns the reversed complementary sequence of DNA or RNA.
       - `gc_content`: returns the GC content of the sequence as a percentage.
       - `is_palindrome`: checks if the sequence is a palindrome.
       - `nucleotide_count`: counts the number of each type of nucleotide in the sequence.

2. **filter_fastq**: filters FASTQ format sequences depending on set parameters.

    - Supported parameters:

       - Sequence length
       - GC content
       - Sequence quality

3. **convert_multiline_fasta_to_oneline**: takes fasta file where nucleic acid sequences are on multiple lines and creates new fasta file where every sequence is on a single line.

4. **parse_blast_output**: takes BLAST result file and adds the sorted list of the first protein names from each query to a new file.

5. **select_genes_from_gbk_to_fasta**: selects set genes and set number of genes before and after target gene and converts its names and translation sequences into new fasta file. 
      
## How to install

___

- Use Git to clone the repository to your local machine.

```bash
git clone git@github.com:vvbasova/SeqSmith.git
```
- Add the path to the repository in your script using sys:

```python
import sys

sys.path.append('/path/to/repository/SeqSmith')

import seqsmith
```

## How to use

---

1. **run_dna_rna_tools** function
   
   - Takes any number of nucleic acid sequences with the last argument corresponding to the selected procedure.

   - If multiple sequences are provided, the function will return a list with the result of the operation for each sequence. If a single sequence is provided, the result will be a string/dictionary/number (depending on the procedure).

**Example:**

```python
run_dna_rna_tools("ATG", "GTa", "transcribe") #['AUG', 'GUa']
run_dna_rna_tools("GTa", "reverse") #aTG
run_dna_rna_tools("GTaAaaaTTTc", "complement") #CAtTtttAAAg
run_dna_rna_tools("GTaAaaaTTTc", "reverse_complement") #gAAAtttTtAC
run_dna_rna_tools("GTaAaaaTTTc", "gc_content") #18.182
run_dna_rna_tools("ATGCAT", "GTa", "is_palindrome") #[True, False]
run_dna_rna_tools("GTaAaaaTTTc", "nucleotide_count") #{'G': 1, 'T': 4, 'A': 5, 'C': 1}
```

2. **filter_fastq** function
   
   - Argument `input_fastq` sets the name to the .fastq file that needs to be filtered.
   - Argument `output_fastq` sets the name for the new file containing filtered .fastq file content.
   - The `gc_bounds` argument sets the GC content range (default is `(0, 100)`), while `length_bounds` sets the acceptable length range for sequences (default is `(0, 2**32)`). If a single number is provided for these arguments, it will be treated as the upper limit.
   - The `quality_threshold` argument sets the minimum read quality threshold (default is 0).
   - Creates `<output_fastq>` file in the `filtered` directory located in a directory where script is executed.
   - If the directory `filtered` does not exist, it will be created.


**Example:**

```python
filter_fastq(input_fastq = "file.fastq", output_fastq = "new_file.fastq", gc_bounds=(40, 100), length_bounds=50, quality_threshold=30)
```

3. **convert_multiline_fasta_to_oneline** function
   
   - Creates new fasta file in working directory
   - `output_fasta` argument is optional: if it was not passed to the function, new file will be named `<file>_one_line.fasta`
   - If file with specified name already exists in working directory, raises FileExistsError

**Example:**

```python
convert_multiline_fasta_to_oneline(input_fasta='file.fasta', output_fasta='oneline_file.fasta')
```
4. **parse_blast_output** function

   - Creates new txt file in working directory
   - `output_file` argument is optional: if it was not passed to the function, new file will be named `<file>_parsed.fasta`
   - If file with specified name already exists in working directory, raises FileExistsError


**Example:**

```python
parse_blast_output(input_file='file.txt', output_file='new_file.txt')
```

5. **select_genes_from_gbk_to_fasta** function
    
   - Creates new fasta file in working directory
   - Arguments `n_before` and `n_after` are 1 by default.
   - `output_fasta` argument is optional: if it was not passed to the function, new file will be named `<file>_selected.fasta`
   - If file with specified name already exists in working directory, raises FileExistsError

**Example:**

```python
select_genes_from_gbk_to_fasta(input_gbk='file.gbk', output_fasta='new_file.fasta', genes=['yicJ_1', 'barA'], n_before=3, n_after=2)
```
## System Requirements

---

- Python 3.x
- `os` module installed

## Project Structure

```
SeqSmith/
├── seqsmith.py                # run_dna_rna_tools and filter_fastq functions
├── bio_files_processor.py     # Reading bioinf files functions          
├── utils/                     # Folder with additional functions
│   ├── dna_rna_utils.py       # Additional functions for run_dna_rna_tools          
│   └── filter_fastq_utils.py  # Additional functions for filter_fastq 
└── README.md                  # Project description

```