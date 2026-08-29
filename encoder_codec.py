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


