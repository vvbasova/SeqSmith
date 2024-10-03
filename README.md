

# SeqSmith

- [russian version](#russian-version)
- [english version](#english-version)
---

### russian version

---

**SeqSmith** — это простая программа для работы с последовательностями нуклеиновых кислот, а также фильтрации данных в формате FASTQ.



## Оглавление

---
- [Основные функции](#основные-функции)
- [Установка](#установка)
- [Использование](#использование)
- [Системные требования](#системные-требования)
- [Структура проекта](#структура-проекта)

## Основные функции

___

1. **run_dna_rna_tools**: позволяет выполнять различные операции с последовательностями ДНК или РНК.

    - Поддерживаемые операции: 

       - `transcribe`: возвращает результат транскрипции ДНК в РНК или обратной транскрипции РНК в ДНК
       - `reverse`: возвращает развёрнутую последовательность ДНК или РНК
       - `complement`: возвращает последовательность ДНК, комплементарную введённой ДНК, или последовательность РНК, комплементарную введённой РНК.
       - `reverse complement`: возвращает развёрную комплементарную последовательность ДНК или РНК
       - `gc_content`: возвращает содержание GC последовательности в процентах
       - `is_palindrome`: проверка, является ли последовательность палиндромом
       - `nucleotide_count`: считает количество каждого типа нуклеотидов в последовательности

2. **filter_fastq**: фильтрует последовательности формата FASTQ в зависимости от заданных параметров.

    - Поддерживаемые параметры:
       
       - Длина последовательностей
       - Содержание GC
       - Качество последовательностей
      
## Установка

___

- Используя Git, клонируйте репозиторий на ваш локальный компьютер.

```bash
git clone git@github.com:vvbasova/SeqSmith.git
```
- Добавьте путь к репозиторию в ваш скрипт, используя sys:

```python
import sys

sys.path.append('/path/to/repository/SeqSmith')

import seqsmith
```

## Использование

---

1. Функция **run_dna_rna_tools** 
   - Принимает любое количество последовательностей нуклеиновых кислот с последним аргументом, соответствующим выбранной процедуре. 

   - Если было введено несколько последовательностей, функция вернёт список с результатом операции для каждой. Если была введена одна последовательность, результатом будет строка/словарь/число (в зависимости от процедуры).

**Пример использования:**

```python
run_dna_rna_tools("ATG", "GTa", "transcribe") #['AUG', 'GUa']
run_dna_rna_tools("GTa", "reverse") #aTG
run_dna_rna_tools("GTaAaaaTTTc", "complement") #CAtTtttAAAg
run_dna_rna_tools("GTaAaaaTTTc", "reverse_complement") #gAAAtttTtAC
run_dna_rna_tools("GTaAaaaTTTc", "gc_content") #18.182
run_dna_rna_tools("ATGCAT", "GTa", "is_palindrome") #[True, False]
run_dna_rna_tools("GTaAaaaTTTc", "nucleotide_count") #{'G': 1, 'T': 4, 'A': 5, 'C': 1}
```

2. Функция **filter_fastq**
   
   - Принимает в качестве первого аргумента словарь, соответствующий формату файла FASTQ, где ключ — ID последовательности, значения — кортеж из последовательности и качества прочтения.
   - Аргумент `gc_bounds` устанавливает диапазон содержания GC (по умолчанию `(0, 100)`), `length_bounds` — диапазон допустимой длины последовательностей (по умолчанию `(0, 2**32)`). Если задать в качестве этих аргументов одно число, оно будет воспринято как верхняя граница.
   - Аргумент `quality_threshold` устанавливает минимальный порог качества прочтения (по умолчанию `0`)
   - Возвращает отфильтрованный словарь
   
**Пример использования:**

```python
example_fastq = {
    '@SRX079801': ('ACAGCAACATAAAC', 'FGGGFGGGFGGGF'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGAT', 'BFFFFFFFB@B@A<@D>BDDACDDDE'),
    '@SRX079803': ('GAACGACAGCAGC', 'DFFFEGDGGGGFGD'),
    '@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCC', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEF'),
}

filter_fastq(example_fastq, gc_bounds=(40, 100), length_bounds=(20, 100), quality_threshold=30)
# {'@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGAT', 'BFFFFFFFB@B@A<@D>BDDACDDDE'), '@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCC', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEF')}
```

## Системные требования

---

- Python 3.9+

## Структура проекта

```
SeqSmith/
├── seqsmith.py                # Основные функции
├── utils/                     # Папка с дополнительными функциями
│   ├── dna_rna_utils.py       # Дополнительные функции для run_dna_rna_tools          
│   └── filter_fastq_utils.py  # Дополнительные функции для filter_fastq 
└── README.md                  # Описание проекта
```
---

### english version

---


# SeqSmith



---

**SeqSmith** — is a simple program for working with nucleic acid sequences and filtering data in FASTQ format.


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
   
   - Takes a dictionary as the first argument that corresponds to the FASTQ file format, where the key is the sequence ID, and the values are tuples of the sequence and read quality.
   - The `gc_bounds` argument sets the GC content range (default is `(0, 100)`), while `length_bounds` sets the acceptable length range for sequences (default is `(0, 2**32)`). If a single number is provided for these arguments, it will be treated as the upper limit.
   - The `quality_threshold` argument sets the minimum read quality threshold (default is 0).
   - Returns a filtered dictionary.
   
**Example:**

```python
example_fastq = {
    '@SRX079801': ('ACAGCAACATAAAC', 'FGGGFGGGFGGGF'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGAT', 'BFFFFFFFB@B@A<@D>BDDACDDDE'),
    '@SRX079803': ('GAACGACAGCAGC', 'DFFFEGDGGGGFGD'),
    '@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCC', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEF'),
}

filter_fastq(example_fastq, gc_bounds=(40, 100), length_bounds=(20, 100), quality_threshold=30)
# {'@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGAT', 'BFFFFFFFB@B@A<@D>BDDACDDDE'), '@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCC', 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEF')}
```

## System Requirements

---

- Python 3.9+

## Project Structure

```
SeqSmith/
├── seqsmith.py                # Main functions
├── utils/                     # Folder with additional functions
│   ├── dna_rna_utils.py       # Additional functions for run_dna_rna_tools          
│   └── filter_fastq_utils.py  # Additional functions for filter_fastq 
└── README.md                  # Project description

```