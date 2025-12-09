

CODON_TABLE = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",  
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C",
    "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "UAA": "STOP", "UAG": "STOP", "UGA": "STOP",
}

def clean_dna(seq: str) -> str:
    """Keep only A/T/G/C, uppercase."""
    return "".join(b for b in seq.upper() if b in "ATGC")

def transcribe_dna_to_rna(dna: str) -> str:
    dna = clean_dna(dna)
    return dna.replace("T", "U")

def translate_rna_to_protein(rna: str) -> str:
    """
    Translate RNA from the first AUG (start codon) until a stop codon or end.
    Returns a single protein string (simple ORF).
    """
    rna = rna.upper()

    
    start = rna.find("AUG")
    if start == -1:
        return ""  

    protein = []
    for i in range(start, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) < 3:
            break
        aa = CODON_TABLE.get(codon, "")
        if aa == "STOP":
            break
        if aa:
            protein.append(aa)
    return "".join(protein)
