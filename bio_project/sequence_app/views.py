import matplotlib.pyplot as plt
from io import BytesIO
import base64
from django.shortcuts import render, redirect, get_object_or_404
from .models import DNASequence
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect
from .models import SavedSequence




# =======================
#  HELPER FUNCTIONS
# =======================

# Simple function to calculate the GC content (returns ratio 0–1)
def calculate_gc_content(seq):
    seq = seq.upper()
    g_count = seq.count('G')
    c_count = seq.count('C')
    total_bases = len(seq)
    
    if total_bases == 0:
        return 0.0
        
    return (g_count + c_count) / total_bases


def reverse_complement(seq):
    seq = seq.upper()
    complement = str.maketrans("ATCG", "TAGC")
    return seq.translate(complement)[::-1]


# --- NEW: transcription + translation helpers ---

CODON_TABLE = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",  # start codon (Methionine)
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


def transcribe_dna_to_rna(dna: str) -> str:
    """
    DNA (A,T,G,C) -> RNA (A,U,G,C)
    Assumes input already validated to only contain A/C/G/T.
    """
    dna = (dna or "").upper()
    return dna.replace("T", "U")


def translate_rna_to_protein(rna: str) -> str:
    """
    Translate RNA from the first AUG (start codon) until a stop codon or end.
    Returns a simple single protein sequence (no frameshifts/complex ORFs).
    """
    rna = (rna or "").upper()
    start = rna.find("AUG")
    if start == -1:
        return ""  # no start codon found

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




def analyze_sequence(request):
    if request.method == 'POST':
        fasta_file = request.FILES.get('fasta_file')
        manual_sequence = request.POST.get('dna_sequence', '')

        raw_sequence = ""

        # PRIORITY 1: File upload
        if fasta_file:
            try:
                file_content = fasta_file.read().decode('utf-8')
                lines = file_content.splitlines()
                raw_sequence = "".join(
                    line.strip() for line in lines if not line.startswith(">")
                )
            except Exception as e:
                return render(request, 'sequence_app/sequence_form.html', {
                    'error': f'Error reading file: {e}'
                })

        # PRIORITY 2: Manual input
        elif manual_sequence:
            raw_sequence = manual_sequence

        # Clean up sequence
        raw_sequence = "".join(raw_sequence.split()).upper()

        if not raw_sequence:
            return render(request, 'sequence_app/sequence_form.html', {
                'error': 'Please enter a sequence or upload a FASTA file.'
            })

        allowed = set("ACGT")
        if not set(raw_sequence).issubset(allowed):
            return render(request, 'sequence_app/sequence_form.html', {
                'error': 'Sequence contains invalid characters. Only A, C, G, T allowed.'
            })

        # ========== PERFORM ANALYSIS ==========
        gc_content = calculate_gc_content(raw_sequence)
        seq_length = len(raw_sequence)
        rev_comp = reverse_complement(raw_sequence)
        rna_seq = transcribe_dna_to_rna(raw_sequence)
        protein_seq = translate_rna_to_protein(rna_seq)

        # Save extra info in session
        request.session['reverse_complement'] = rev_comp
        request.session['rna_sequence'] = rna_seq
        request.session['protein_sequence'] = protein_seq

        # ========== SAVE ONLY IF LOGGED IN ==========
        if request.user.is_authenticated:
            new_analysis = DNASequence.objects.create(
                user=request.user,
                sequence=raw_sequence,
                length=seq_length,
                gc_content=gc_content
            )
            return redirect('results_view', analysis_id=new_analysis.id)

        # ========== SHOW RESULTS WITHOUT SAVING ==========
        context = {
            'sequence': raw_sequence,
            'length': seq_length,
            'gc_percentage': f"{gc_content * 100:.2f}%",
            'submission_date': None,
            'reverse_complement': rev_comp,
            'rna_sequence': rna_seq,
            'protein_sequence': protein_seq,
            'base_chart': base_count_chart(raw_sequence),
        }
        return render(request, 'sequence_app/results.html', context)

    # GET → show form
    return render(request, 'sequence_app/sequence_form.html')


def results_view(request, analysis_id):
    # Retrieve the saved object from PostgreSQL, or return 404 if not found
    analysis = get_object_or_404(DNASequence, id=analysis_id)

    # Get reverse complement, RNA, and protein from session
    reverse_comp = request.session.get('reverse_complement', '')
    rna_seq = request.session.get('rna_sequence', '')
    protein_seq = request.session.get('protein_sequence', '')

    # Create base count chart for the saved sequence
    chart_data_uri = base_count_chart(analysis.sequence)
    
    context = {
        'sequence': analysis.sequence,
        'length': analysis.length,
        # Format the float ratio to a percentage string (e.g., 45.00%) for display
        'gc_percentage': f"{analysis.gc_content * 100:.2f}%", 
        'submission_date': analysis.submission_date,
        'reverse_complement': reverse_comp,
        'base_chart': chart_data_uri,
        # NEW:
        'rna_sequence': rna_seq,
        'protein_sequence': protein_seq,
    }
    return render(request, 'sequence_app/results.html', context)


def primer_designer(request):
    """
    Very simple primer designer:
    - Takes a DNA sequence (from ?sequence= or from form POST)
    - Uses first 20 bases as forward primer
    - Uses reverse complement of last 20 bases as reverse primer
    - Calculates a simple Tm for each primer
    """
    # Get sequence from query parameter or POST form
    seq = request.GET.get("sequence", "") or request.POST.get("sequence", "")
    seq = "".join((seq or "").split()).upper()  # remove spaces/newlines, uppercase

    context = {
        "sequence": seq,
        "forward_primer": None,
        "reverse_primer": None,
        "tm_forward": None,
        "tm_reverse": None,
        "error": None,
    }

    if request.method == "POST":
        if not seq:
            context["error"] = "Please provide a DNA sequence."
        else:
            # simple validation: only allow A,C,G,T
            allowed = set("ACGT")
            if not set(seq).issubset(allowed):
                context["error"] = "Sequence contains invalid characters. Only A, C, G, T are allowed."
            elif len(seq) < 40:
                context["error"] = "Sequence is too short for primer design (need at least 40 bases)."
            else:
                # super simple demo: first 20 and last 20
                fwd = seq[:20]
                rev = reverse_complement(seq[-20:])  # uses your existing helper

                # simple Wallace rule: Tm = 2*(A+T) + 4*(G+C)
                def simple_tm(primer: str) -> float:
                    primer = primer.upper()
                    a = primer.count("A")
                    t = primer.count("T")
                    g = primer.count("G")
                    c = primer.count("C")
                    return 2 * (a + t) + 4 * (g + c)

                context["forward_primer"] = fwd
                context["reverse_primer"] = rev
                context["tm_forward"] = simple_tm(fwd)
                context["tm_reverse"] = simple_tm(rev)

    return render(request, "sequence_app/primer_designer.html", context)


def motif_finder(request):
    sequence = request.POST.get("sequence", "").upper().replace(" ", "")
    motif = request.POST.get("motif", "").upper()

    positions = []
    error = None

    if request.method == "POST":
        if not sequence:
            error = "Please enter a DNA sequence."
        elif not motif:
            error = "Please enter a motif to search."
        else:
            # find all occurrences (including overlapping)
            idx = 0
            while True:
                idx = sequence.find(motif, idx)
                if idx == -1:
                    break
                positions.append(idx + 1)  # 1-based index
                idx += 1

            if not positions:
                error = "Motif not found in sequence."

    return render(request, "sequence_app/motif_finder.html", {
        "sequence": sequence,
        "motif": motif,
        "positions": positions,
        "error": error,
    })

STOP_CODONS = {"TAA", "TAG", "TGA"}

def find_orfs_internal(dna, min_aa_len=20):
    dna = dna.upper()
    dna = "".join(dna.split())        # remove whitespace
    n = len(dna)
    rc = reverse_complement(dna)
    orfs = []

    # Scan one sequence in one strand and frame
    def scan(seq, strand, frame):
        i = frame
        while i + 3 <= len(seq):
            codon = seq[i:i+3]
            if codon == "ATG":
                j = i + 3
                while j + 3 <= len(seq):
                    stop = seq[j:j+3]
                    if stop in STOP_CODONS:
                        nt_len = (j + 3) - i
                        aa_len = nt_len // 3

                        if aa_len >= min_aa_len:
                            if strand == "+":
                                start = i + 1
                                end = j + 3
                            else:
                                # reverse complement coordinate conversion
                                start = n - (j + 3) + 1
                                end = n - i

                            orfs.append({
                                "strand": strand,
                                "frame": frame,
                                "start": start,
                                "end": end,
                                "length_nt": nt_len,
                                "length_aa": aa_len,
                                "dna": seq[i:j+3],
                            })
                        break
                    j += 3
                i = j  # move to next scanning point
            else:
                i += 3

    # Forward frames
    for f in range(3):
        scan(dna, "+", f)

    # Reverse complement frames
    for f in range(3):
        scan(rc, "-", f)

    # Sort longest first
    orfs.sort(key=lambda x: x["length_aa"], reverse=True)
    return orfs


def orf_finder(request):
    raw_input = request.POST.get("sequence", "")
    print("RAW INPUT =", repr(raw_input))

    # Clean sequence fully
    sequence = (
        raw_input.upper()
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
    )
    print("CLEANED =", repr(sequence))

    orfs = []
    longest_orf = None
    protein = ""
    error = None

    if request.method == "POST":

        # Validate characters
        allowed = set("ACGT")
        if not set(sequence).issubset(allowed):
            error = "Invalid characters found. Only A, C, G, T allowed."
            return render(request, "sequence_app/orf_finder.html", {
                "sequence": raw_input,
                "orfs": [],
                "longest_orf": None,
                "protein": "",
                "error": error,
            })

        # Allow ANY length ORF ≥ 1 aa
        orfs = find_orfs_internal(sequence, min_aa_len=1)

        if not orfs:
            error = "No ORFs found."
        else:
            longest_orf = orfs[0]
            rna = transcribe_dna_to_rna(longest_orf["dna"])
            protein = translate_rna_to_protein(rna)

    return render(request, "sequence_app/orf_finder.html", {
        "sequence": raw_input,
        "orfs": orfs,
        "longest_orf": longest_orf,
        "protein": protein,
        "error": error,
    })


def restriction_finder(request):
    sequence = request.POST.get("sequence", "").upper().replace(" ", "")

    # Common restriction enzymes + their recognition sites
    enzymes = {
        "EcoRI": "GAATTC",
        "HindIII": "AAGCTT",
        "BamHI": "GGATCC",
        "NotI": "GCGGCCGC",
        "XhoI": "CTCGAG",
        "NheI": "GCTAGC",
        "PstI": "CTGCAG",
        "SmaI": "CCCGGG",
        "KpnI": "GGTACC",
    }

    results = []
    error = None
    highlighted_sequence = sequence  # default: plain (no color)

    if request.method == "POST":
        if not sequence:
            error = "Please enter a DNA sequence."
        else:
            # 1) Find positions for each enzyme (your original logic)
            for enzyme, site in enzymes.items():
                positions = []
                i = 0
                while True:
                    i = sequence.find(site, i)
                    if i == -1:
                        break
                    positions.append(i + 1)  # 1-based index
                    i += 1

                results.append({
                    "enzyme": enzyme,
                    "site": site,
                    "positions": positions,
                    "count": len(positions),
                })

            # 2) Build a colored version of the sequence
            colored = sequence

            # simple approach: replace each site string with a span
            for enzyme, site in enzymes.items():
                colored = colored.replace(
                    site,
                    f"<span style='background:#ffe08a; padding:2px 3px; border-radius:4px; font-weight:600;'>{site}</span>"
                )

            highlighted_sequence = colored

    return render(request, "sequence_app/restriction_finder.html", {
        "sequence": sequence,
        "results": results,
        "error": error,
        "highlighted_sequence": highlighted_sequence,
    })



def base_count_chart(sequence):
    """
    Create a matplotlib bar chart for base counts (A,T,G,C).
    Returns a base64-encoded PNG data URI (string) ready to use in <img src="...">.
    """
    seq = (sequence or "").upper()
    # Count bases
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')

    # Create bar chart (one plot only).
    plt.figure(figsize=(4, 3))
    plt.bar(['A', 'T', 'G', 'C'], [a, t, g, c], color=['#2563EB', '#FBBF24', '#10B981', '#EF4444'])
    plt.title('Base Counts Distribution', fontsize=14)
    plt.xlabel('Base', fontsize=12)
    plt.ylabel('Count', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    # Save to buffer
    buf = BytesIO()
    plt.tight_layout()
    plt.savefig(buf, format='png')
    buf.seek(0)

    # Encode to base64 for embedding
    img_b64 = base64.b64encode(buf.read()).decode('ascii')
    plt.close()
    return f"data:image/png;base64,{img_b64}"


# =======================
#  CRUD VIEWS
# =======================

# READ ALL sequences
@login_required
def sequence_list(request):
    sequences = DNASequence.objects.filter(user=request.user)
    return render(request, 'sequence_app/sequence_list.html', {'sequences': sequences})


# READ ONE sequence
@login_required
def sequence_detail(request, pk):
    sequence = get_object_or_404(DNASequence, pk=pk, user=request.user)
    return render(request, 'sequence_app/sequence_detail.html', {'sequence': sequence})




# CREATE new sequence
@login_required
def sequence_create(request):
    if request.method == 'POST':
        fasta_file = request.FILES.get('fasta_file')
        raw_sequence = ""
        if fasta_file:
            # PRIORITY 1: Handle file content
            try:
                file_content = fasta_file.read().decode('utf-8')
                lines = file_content.splitlines()
                # Clean FASTA headers
                raw_sequence = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
            except Exception as e:
                # If file reading fails
                return render(request, 'sequence_app/sequence_form.html', {"error": f"Error reading file: {e}"})
        
        # PRIORITY 2: If no file content, check manual input (name='sequence')
        if not raw_sequence:
            raw_sequence = request.POST.get('sequence', '').strip().upper()
        
        # Final validation and creation
        if raw_sequence:
            gc_content_val = calculate_gc_content(raw_sequence)
            DNASequence.objects.create(
                user=request.user,
                sequence=raw_sequence, 
                length=len(raw_sequence),
                gc_content=gc_content_val,
             
            )
            return redirect('sequence_list')

        # If we reach here, no data was provided
        return render(request, 'sequence_app/sequence_form.html', {"error": "Please enter a sequence or upload a file."})
    return render(request, 'sequence_app/sequence_form.html')
        

# UPDATE sequence
def sequence_update(request, pk):
    sequence = get_object_or_404(DNASequence, pk=pk, user=request.user)
    if request.method == 'POST':
        seq = request.POST.get('sequence', '').strip().upper()
        if seq:
            gc_content_val = calculate_gc_content(seq)
            sequence.sequence = seq
            sequence.length = len(seq)
            sequence.gc_content = gc_content_val
            sequence.save()
            return redirect('sequence_list')
    return render(request, 'sequence_app/sequence_form.html', {'sequence': sequence})


# DELETE sequence
def sequence_delete(request, pk):
    sequence = get_object_or_404(DNASequence, pk=pk, user=request.user)
    if request.method == 'POST':
        sequence.delete()
        return redirect('sequence_list')
    return render(request, 'sequence_app/sequence_confirm_delete.html', {'sequence': sequence})


@login_required
def save_sequence(request):
    if request.method == 'POST':
        sequence = request.POST.get('sequence')

        if not sequence:
            return redirect('analyze_sequence')  # fallback

        SavedSequence.objects.create(
            user=request.user,
            sequence=sequence
        )

        return redirect('saved_sequences')  # we will create this page later
    
@login_required
def saved_sequences(request):
    saved = SavedSequence.objects.filter(user=request.user).order_by('-created')
    return render(request, "sequence_app/saved_sequences.html", {"saved": saved})


