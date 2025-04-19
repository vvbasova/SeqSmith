

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

1. **Biological Sequence Classes**: Represents biological sequences (DNA, RNA, and amino acids) with support for various operations.

    - Supported operations: 

       - `transcribe()`: returns the result of transcribing DNA to RNA or reverse transcribing RNA to DNA.
       - `reverse()`: returns the reversed sequence of DNA or RNA.
       - `complement()`: returns a DNA sequence complementary to the input DNA, or an RNA sequence complementary to the input RNA.
       - `reverse complement()`: returns the reversed complementary sequence of DNA or RNA.
       - `gc_content()`: returns the GC content of the sequence as a percentage.
       - `is_palindrome()`: checks if the sequence is a palindrome.
       - `nucleotide_count()`: counts the number of each type of nucleotide in the sequence.
       - `molecular_weight()`: counts molecular weight of amino acid sequence
    
    - The following sequence types are supported:
       
       - `DNASequence`: Represents a DNA sequence, with transcription to RNA.
       - `RNASequence`: Represents an RNA sequence.
       - `AminoAcidSequence`: Represents an amino acid sequence, with functionality for calculating molecular weight.

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

1. **Biological Sequence Classes** 
   
   - `BiologicalSequence` class is the base class for all biological sequences. It provides basic functionality for sequence operations, such as length, indexing, and string representation.

**Example:**


```python
my_sequence = DNASequence("ACTG")
len(my_sequence) #4
my_sequence[2] #T
print(my_sequence) #ACTG
print(repr(my_sequence)) #DNASequence('ACTG')
```   

   - `NucleicAcidSequence` class inherits from BiologicalSequence and
includes methods for working with DNA and RNA sequences, such as complementing, reversing, and calculating GC content.


**Example:**


```python
my_sequence = DNASequence("GTaAaaaTTTc")
my_sequence.complement() #CAtTtttAAAg
my_sequence.reverse_complement() #gAAAtttTtAC
my_sequence.gc_content() #18.182
my_sequence.reverse() #cTTTaaaAaTG
my_sequence.nucleotide_count() #{'G': 1, 'T': 4, 'A': 5, 'C': 1}
my_sequence.is_palindrome() #False
```

   - `DNASequence` and `RNASequence` classes inherit from `NucleicAcidSequence` and allows to use its methods specifically for DNA and RNA sequences. `DNASequence` also has method for transcribing DNA to RNA. `transcribe` method returns an RNASequence class instance. Methods `complement()`, `reverse_complement()`, `reverce()` return an instance of corresponding class (`DNAsequence` or `RNASequence`).  

**Example:**

```python
my_dna_sequence = DNASequence("GTaAaaaTTTc")
my_rna = my_dna_sequence.transcribe() 
my_rna #RNASequence('GUaAaaaUUUc')

my_rna_sequence = RNASequence("UCGA")
my_rna_sequence.complement() #AGCU
```

   - The `AminoAcidSequence` class provides as calculating molecular weight operation.

**Example:**

```python
my_sequence = AminoAcidSequence("KEK")
my_sequence.molecular_weight() #439.5
```


2. **filter_fastq** function
   
   - Argument `input_fastq` sets the name to the .fastq file that needs to be filtered.
   - Argument `output_fastq` sets the name for the new file containing filtered .fastq file content.
   - The `gc_bounds` argument sets the GC content range (default is `(0, 100)`), while `length_bounds` sets the acceptable length range for sequences (default is `(0, 2**32)`). If a single number is provided for these arguments, it will be treated as the upper limit.
   - The `quality_threshold` argument sets the minimum read quality threshold (default is 0).
   - Creates `<output_fastq>` file in the `filtered` directory located in a directory where script is executed.
   - If the directory `filtered` does not exist, it will be created.
   - The function logs the execution process  in "./logs/filter_fastq.log" file. If the `log` directory doesn't exist, creates it.
   - The function is available via Terminal


**Example:**

```python
filter_fastq(input_fastq = "file.fastq", output_fastq = "new_file.fastq", gc_bounds=(40, 100), length_bounds=50, quality_threshold=30)
```

```Terminal
python seqsmith.py filter-fastq {input_file_name}.fastq {output_file_name}.fastq -g 50 100 -l 16 20 -q 30

python seqsmith.py filter-fastq {input_file_name}.fastq {output_file_name}.fastq --gc_bounds 50 100 --length-bounds 16 20 --quality-threshold 30
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
- Biopython>=1.57

## Project Structure

```
SeqSmith/
├── seqsmith.py                # Biological Sequence Classes and filter_fastq functions
├── bio_files_processor.py     # Reading bioinf files functions
├── tests/                     # Directory for tests
│   ├── data/
│   │   └── example.fastq      
│   ├── __init__.py
│   └── test_filter_fastq.py   # Tests for filter_fastq function
├── requirements.txt           # System requirements       
└── README.md                  # Project description

```