import os
import pytest
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from seqsmith import filter_fastq


TEST_OUTPUT_FILES = [
    "created_output.fastq",
    "gc_filtered.fastq",
    "length_filtered.fastq",
    "quality_filtered.fastq",
    "gc_len.fastq",
    "combined.fastq",
    "gc_single_value.fastq",
    "len_single_value.fastq",
    "all_zero.fastq",
    "already_exists.fastq",
    "bad_types_1.fastq",
    "bad_types_2.fastq",
    "bad_types_3.fastq"
]


@pytest.fixture(autouse=True)
def cleanup_after():
    yield
    for filename in TEST_OUTPUT_FILES:
        filepath = os.path.join("filtered", filename)
        if os.path.exists(filepath):
            os.remove(filepath)


@pytest.fixture
def example_fastq_path():
    return os.path.join(os.path.dirname(__file__), "data", "example.fastq")


class TestExampleFastq:
    """
    Test suite for the `filter_fastq` function from SeqSmith.
    Tests use example.fastq file with FASTQ examples.
    All generated output files are cleaned up after each test run.

    - test_output_file_created:
        Verifies that the function creates an output
        file when only input and output are provided.
    - test_gc_filtering:
        Applies GC content filter and checks that all
        reads in the output meet the lower GC threshold.
    - test_length_filtering:
        Applies length filter and ensures all reads are
        within the specified length bounds.
    - test_quality_filtering:
        Filters reads by quality and checks that average
        quality is above a given threshold.
    - test_gc_and_length_only:
        Applies both GC and length filters and validates
        that all reads pass both criteria.
    - test_combined_filters:
        Applies GC, length, and quality filters simultaneously
        and checks combined correctness.
    - test_gc_bounds_as_single_value:
        Passes a single GC upper bound as an integer and
        validates the result.
    - test_length_bounds_as_single_value:
        Passes a single length upper bound as an integer and
        validates the result.
    - test_all_zero_parameters:
        Uses zero values for all filters to check if no
        reads are allowed through.
    - test_existing_file_raises_error:
        Ensures the function raises FileExistsError
        if the output file already exists.
    - test_invalid_argument_types:
        Verifies that inappropriate argument types
        (e.g., strings instead of numbers) raise TypeError
        or ValueError.
    """

    def test_output_file_created(self, example_fastq_path):
        output_name = "created_output.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(input_fastq=example_fastq_path, output_fastq=output_name)

        assert os.path.exists(out_path)
        records = list(SeqIO.parse(out_path, "fastq"))
        assert len(records) > 0

    def test_gc_filtering(self, example_fastq_path):
        output_name = "gc_filtered.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            gc_bounds=(60, 100)
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert len(records) > 0
        assert all(gc_fraction(r.seq) * 100 >= 60 for r in records)

    def test_length_filtering(self, example_fastq_path):
        output_name = "length_filtered.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            length_bounds=(80, 200)
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert len(records) > 0
        assert all(80 <= len(r.seq) <= 200 for r in records)

    def test_quality_filtering(self, example_fastq_path):
        output_name = "quality_filtered.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            quality_threshold=30
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert len(records) > 0
        for r in records:
            avg_q = sum(r.letter_annotations["phred_quality"]) / len(r)
            assert avg_q >= 30

    def test_gc_and_length_only(self, example_fastq_path):
        output_name = "gc_len.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(input_fastq=example_fastq_path, output_fastq=output_name,
                     gc_bounds=(40, 70), length_bounds=(50, 120))

        records = list(SeqIO.parse(out_path, "fastq"))
        assert all(40 <= gc_fraction(r.seq) * 100 <= 70 for r in records)
        assert all(50 <= len(r.seq) <= 120 for r in records)

    def test_combined_filters(self, example_fastq_path):
        output_name = "combined.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            gc_bounds=(40, 70),
            length_bounds=(80, 150),
            quality_threshold=20
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        for r in records:
            assert 40 <= gc_fraction(r.seq) * 100 <= 70
            assert 80 <= len(r.seq) <= 150
            avg_q = sum(r.letter_annotations["phred_quality"]) / len(r)
            assert avg_q >= 20

    def test_gc_bounds_as_single_value(self, example_fastq_path):
        output_name = "gc_single_value.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            gc_bounds=60
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert all(gc_fraction(r.seq) * 100 <= 60 for r in records)

    def test_length_bounds_as_single_value(self, example_fastq_path):
        output_name = "len_single_value.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            length_bounds=100
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert all(len(r.seq) <= 100 for r in records)

    def test_all_zero_parameters(self, example_fastq_path):
        output_name = "all_zero.fastq"
        out_path = os.path.join("filtered", output_name)

        filter_fastq(
            input_fastq=example_fastq_path,
            output_fastq=output_name,
            gc_bounds=(0, 0),
            length_bounds=(0, 0),
            quality_threshold=0
        )

        records = list(SeqIO.parse(out_path, "fastq"))
        assert len(records) == 0

    def test_existing_file_raises_error(self, example_fastq_path):
        output_name = "already_exists.fastq"
        out_path = os.path.join("filtered", output_name)
        os.makedirs("filtered", exist_ok=True)
        with open(out_path, "w") as f:
            f.write("already exists")

        with pytest.raises(FileExistsError):
            filter_fastq(
                input_fastq=example_fastq_path,
                output_fastq=output_name
            )

    def test_invalid_argument_types(self, example_fastq_path):
        with pytest.raises((TypeError, ValueError)):
            filter_fastq(
                input_fastq=example_fastq_path,
                output_fastq="bad_types_1.fastq",
                gc_bounds=("low", "high"),
                length_bounds=(0, 10),
                quality_threshold=0
            )

        with pytest.raises((TypeError, ValueError)):
            filter_fastq(
                input_fastq=example_fastq_path,
                output_fastq="bad_types_2.fastq",
                gc_bounds=(0, 100),
                length_bounds="short",
                quality_threshold=0
            )

        with pytest.raises((TypeError, ValueError)):
            filter_fastq(
                input_fastq=example_fastq_path,
                output_fastq="bad_types_3.fastq",
                gc_bounds=(0, 100),
                length_bounds=(0, 10),
                quality_threshold="high"
            )
