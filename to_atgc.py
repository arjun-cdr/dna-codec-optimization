# ==========================================
# DNA ENCODER
# Binary -> ACGT
# ==========================================

MIN_GC = 40
MAX_GC = 60
MAX_HOMOPOLYMER = 4

BINARY_TO_DNA = {
    "00": "A",
    "01": "C",
    "10": "G",
    "11": "T"
}


# ---------- Read Binary File ----------

def read_binary_file(filename):
    with open(filename, "r") as file:
        binary = file.read().replace(" ", "").replace("\n", "")

    # Make sure the input contains only 0s and 1s
    if not all(bit in "01" for bit in binary):
        raise ValueError("Input file contains characters other than 0 and 1.")

    return binary


# ---------- Binary -> DNA ----------

def binary_to_dna(binary):
    dna = ""

    for i in range(0, len(binary), 2):
        bits = binary[i:i + 2]

        # Binary data should normally be divisible by 2
        if len(bits) != 2:
            raise ValueError("Binary data length must be even.")

        dna += BINARY_TO_DNA[bits]

    return dna


# ---------- GC Content Check ----------

def calculate_gc_content(dna):
    gc_count = dna.count("G") + dna.count("C")

    gc_percentage = (gc_count / len(dna)) * 100

    return gc_percentage


def check_gc_constraint(dna):
    gc_percentage = calculate_gc_content(dna)

    if MIN_GC <= gc_percentage <= MAX_GC:
        return True

    return False


# ---------- Homopolymer Check ----------

def check_homopolymer(dna):
    current_base = dna[0]
    current_length = 1

    for base in dna[1:]:

        if base == current_base:
            current_length += 1

            if current_length > MAX_HOMOPOLYMER:
                return False

        else:
            current_base = base
            current_length = 1

    return True


# ---------- Complete Constraint Check ----------

def validate_dna(dna):

    if not check_gc_constraint(dna):
        return False

    if not check_homopolymer(dna):
        return False

    return True


# ---------- Save DNA ----------

def save_dna_file(dna, filename):

    with open(filename, "w") as file:
        file.write(dna)


# ---------- Main Encoder ----------

def encode_file(binary_file, dna_file):

    # 1. Read binary
    binary = read_binary_file(binary_file)

    # 2. Convert binary -> DNA
    dna = binary_to_dna(binary)

    # 3. Check constraints
    if not validate_dna(dna):

        print("DNA sequence failed the constraints.")

        print(
            f"GC Content: "
            f"{calculate_gc_content(dna):.2f}%"
        )

        print(
            f"Maximum allowed homopolymer: "
            f"{MAX_HOMOPOLYMER}"
        )

        return False

    # 4. Save DNA
    save_dna_file(dna, dna_file)

    print("DNA encoding successful.")
    print(f"Binary length : {len(binary)} bits")
    print(f"DNA length    : {len(dna)} nt")
    print(
        f"Density       : "
        f"{len(binary) / len(dna):.2f} bits/nt"
    )
    print(
        f"GC Content    : "
        f"{calculate_gc_content(dna):.2f}%"
    )

    return True


# ---------- Run ----------

encode_file(
    "binary_output.txt",
    "dna_output.txt"
)