import json
import hashlib
from collections import Counter

# ============================================================
# CONFIGURATION
# ============================================================

INPUT_FILE = "binary_output.txt"

DNA_OUTPUT_FILE = "final_dna_sequence.txt"
MAPPING_OUTPUT_FILE = "final_mapping.json"

# Constraint parameters
MIN_GC = 0.40
MAX_GC = 0.60

MAX_HOMOPOLYMER = 3

# k-mer constraint
KMER_SIZE = 4
MAX_KMER_FREQUENCY = 3

# Maximum number of repeated consecutive patterns
MAX_REPEAT_COUNT = 2
MIN_REPEAT_LENGTH = 4

# ============================================================
# BASIC BINARY VALIDATION
# ============================================================

def validate_binary(binary):
    """
    Validates binary input and removes whitespace.

    Accepts:
        01010100 01101000 01101001

    Converts it to:
        010101000110100001101001
    """

    # Remove spaces, newlines, tabs, etc.
    binary = "".join(binary.split())

    if not binary:
        raise ValueError("Binary input is empty.")

    if any(bit not in "01" for bit in binary):
        raise ValueError(
            "Input contains characters other than 0 and 1."
        )

    return binary


# ============================================================
# BINARY PREPARATION
# ============================================================

def prepare_binary(binary):
    """
    Ensures the binary length is divisible by 2.

    2 bits -> 1 nucleotide
    """

    binary = validate_binary(binary)

    padding = 0

    if len(binary) % 2 != 0:
        binary += "0"
        padding = 1

    return binary, padding

# ============================================================
# DNA CONVERSION UTILITIES
# ============================================================

BASES = ["A", "C", "G", "T"]

TWO_BIT_VALUES = ["00", "01", "10", "11"]


def bits_to_pairs(binary):
    """
    Converts:

    00101100

    into:

    ['00', '10', '11', '00']
    """

    return [
        binary[i:i + 2]
        for i in range(0, len(binary), 2)
    ]


# ============================================================
# CONSTRAINT 1: GC CONTENT
# ============================================================

def check_gc_content(sequence):
    """
    Checks whether GC content is between MIN_GC and MAX_GC.
    """

    if not sequence:
        return False, 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    gc_content = gc_count / len(sequence)

    passed = MIN_GC <= gc_content <= MAX_GC

    return passed, gc_content

# ============================================================
# CONSTRAINT 2: HOMOPOLYMER
# ============================================================

def check_homopolymer(sequence):
    """
    Checks for excessive consecutive identical bases.

    Example:

    AAA  -> allowed
    AAAA -> violation
    """

    if not sequence:
        return True, 0

    longest = 1
    current = 1

    for i in range(1, len(sequence)):

        if sequence[i] == sequence[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1

    passed = longest <= MAX_HOMOPOLYMER

    return passed, longest

# ============================================================
# CONSTRAINT 3: SEQUENCE COMPLEXITY / k-MER FREQUENCY
# ============================================================

def check_kmer_frequency(sequence):
    """
    Checks whether any k-mer occurs too many times.
    """

    if len(sequence) < KMER_SIZE:
        return True, {}

    kmers = [
        sequence[i:i + KMER_SIZE]
        for i in range(len(sequence) - KMER_SIZE + 1)
    ]

    frequencies = Counter(kmers)

    violations = {
        kmer: count
        for kmer, count in frequencies.items()
        if count > MAX_KMER_FREQUENCY
    }

    passed = len(violations) == 0

    return passed, violations

# ============================================================
# CONSTRAINT 4: REPEATED SEQUENCES
# ============================================================

def check_repeated_sequences(sequence):
    """
    Checks for excessive TANDEM repeats.

    A tandem repeat occurs when the same block appears
    immediately twice:

        ACGTACGT
        ---- ----
         ACGT ACGT

    We only consider repeats of length 4 or greater.
    Short 2- and 3-base repetitions are not treated as
    constraint violations.
    """

    MIN_REPEAT_LENGTH = 4
    MAX_REPEAT_LENGTH = 8

    violations = []

    for repeat_length in range(4, 9):

        for i in range(
            len(sequence) - 2 * repeat_length + 1
        ):

            block1 = sequence[
                i:i + repeat_length
            ]

            block2 = sequence[
                i + repeat_length:
                i + 2 * repeat_length
            ]

            if block1 == block2:

                violations.append({
                    "position": i,
                    "repeat_length": repeat_length,
                    "sequence": block1
                })

    passed = len(violations) == 0

    return passed, violations

# ============================================================
# ALL CONSTRAINTS
# ============================================================

def check_all_constraints(sequence):
    """
    Executes all constraints.

    Returns:
        overall_pass
        detailed_results
    """

    gc_pass, gc_value = check_gc_content(sequence)

    homopolymer_pass, longest_homopolymer = \
        check_homopolymer(sequence)

    kmer_pass, kmer_violations = \
        check_kmer_frequency(sequence)

    repeat_pass, repeat_violations = \
        check_repeated_sequences(sequence)

    results = {

        "GC_Content": {
            "passed": gc_pass,
            "value": gc_value,
            "allowed_range": [
                MIN_GC,
                MAX_GC
            ]
        },

        "Homopolymer": {
            "passed": homopolymer_pass,
            "longest_run": longest_homopolymer,
            "maximum_allowed": MAX_HOMOPOLYMER
        },

        "kmer_frequency": {
            "passed": kmer_pass,
            "k": KMER_SIZE,
            "maximum_allowed": MAX_KMER_FREQUENCY,
            "violations": kmer_violations
        },

        "Repeated_sequences": {
            "passed": repeat_pass,
            "violation_count": len(repeat_violations),
            "maximum_allowed": MAX_REPEAT_COUNT
        }
    }

    overall_pass = (
        gc_pass
        and homopolymer_pass
        and kmer_pass
        and repeat_pass
    )

    return overall_pass, results

# ============================================================
# MAPPING 1: DIRECT MAPPING
# ============================================================

def direct_mapping(binary):
    """
    Direct mapping:

    00 -> A
    01 -> C
    10 -> G
    11 -> T
    """

    mapping = {
        "00": "A",
        "01": "C",
        "10": "G",
        "11": "T"
    }

    pairs = bits_to_pairs(binary)

    sequence = "".join(
        mapping[pair]
        for pair in pairs
    )

    metadata = {
        "technique": "Direct Mapping",
        "mapping": mapping
    }

    return sequence, metadata

# ============================================================
# MAPPING 2: PERMUTATION BASED MAPPING
# ============================================================

def permutation_mapping(binary):
    """
    Uses a different permutation of A,C,G,T.

    00 -> G
    01 -> T
    10 -> A
    11 -> C
    """

    mapping = {
        "00": "G",
        "01": "T",
        "10": "A",
        "11": "C"
    }

    pairs = bits_to_pairs(binary)

    sequence = "".join(
        mapping[pair]
        for pair in pairs
    )

    metadata = {
        "technique": "Permutation Based Mapping",
        "mapping": mapping
    }

    return sequence, metadata



# ============================================================
# PROGRAM ENTRY POINT
# ============================================================

if __name__ == "__main__":

    encode_file()
