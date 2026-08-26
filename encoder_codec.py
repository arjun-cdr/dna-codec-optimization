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
