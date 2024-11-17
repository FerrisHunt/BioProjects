import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import random
from Bio import pairwise2

# Streamlit page configuration
st.set_page_config(page_title="DNA Sequence Analyzer", layout="wide")

# Display the title at the top of the page
st.title("DNA Sequence Analyzer")

# DNA sequence input field (always visible)
sequence = st.text_area("Enter DNA Sequence (ATCG)", "")

# Validate DNA sequence input (Remove spaces, newlines, and check for non-ATCG characters)
def clean_sequence(seq):
    seq = seq.replace(" ", "")  # Remove spaces from the sequence
    seq = seq.replace("\n", "")  # Remove newlines from the sequence
    valid_bases = "ATCGatcg"
    invalid_bases = [base for base in seq if base not in valid_bases]
    seq = seq.upper()  # Convert the sequence to uppercase for consistent processing
    return seq, invalid_bases

# Functions for analysis
def nucleotide_frequencies(seq):
    seq, invalid_bases = clean_sequence(seq)  # Clean the sequence first
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "C": seq.count("C"),
        "G": seq.count("G")
    }, invalid_bases

def wallace_tm(seq):
    seq = seq.upper()  # Ensure uppercase for Tm calculation
    A_T_count = seq.count('A') + seq.count('T')
    G_C_count = seq.count('G') + seq.count('C')
    tm_celsius = 2 * A_T_count + 4 * G_C_count
    return tm_celsius  # Return only in Celsius

# Display the sequence with the codon or peptide highlighted in red and position displayed
def display_alignment(seq, search_sequence=""):
    highlighted_sequence = ""
    positions = []
    start = 0
    while start < len(seq):
        pos = seq.find(search_sequence, start)
        if pos == -1:
            highlighted_sequence += seq[start:]
            break
        end = pos + len(search_sequence)
        highlighted_sequence += seq[start:pos]
        highlighted_sequence += f"<span style='color:red;cursor:pointer;' title='Position: {pos+1}-{end}'> {seq[pos:end]} </span>"
        positions.append(f"{pos+1}-{end}")
        start = end

    return highlighted_sequence, positions

# Initialize session state variables if they are not initialized
if "show_tm" not in st.session_state:
    st.session_state.show_tm = False
if "show_rna" not in st.session_state:
    st.session_state.show_rna = False
if "show_cDNA" not in st.session_state:
    st.session_state.show_cDNA = False
if "show_protein" not in st.session_state:
    st.session_state.show_protein = False
if "show_length" not in st.session_state:
    st.session_state.show_length = False
if "show_nucleotide_freq" not in st.session_state:
    st.session_state.show_nucleotide_freq = False
if "show_gc" not in st.session_state:
    st.session_state.show_gc = False
if "sequence_submitted" not in st.session_state:
    st.session_state.sequence_submitted = False

# Store the submit button container
submit_container = st.empty()

# Handle sequence submission
if not st.session_state.sequence_submitted:
    if submit_container.button("Submit"):
        if sequence:
            # Clean the sequence (remove spaces, newlines, and check for invalid characters)
            cleaned_sequence, invalid_bases = clean_sequence(sequence)
            st.session_state.sequence_submitted = True  # Set flag to hide the submit button
            submit_container.empty()  # Hide the submit button container

            # If there are invalid characters, display an error message
            if invalid_bases:
                st.error(f"Your sequence contains invalid nucleotide(s): {', '.join(set(invalid_bases))}. Please enter only A, T, C, and G.")
            else:
                # Proceed with analysis if valid
                sequence = cleaned_sequence
        else:
            # Show the message if the user clicked "Submit" without entering a sequence
            st.write("**Please enter a DNA sequence first!**")

# If sequence is already submitted, show analysis options
if st.session_state.sequence_submitted:
    # Sequence Length
    if st.button("Sequence Length"):
        st.session_state.show_length = not st.session_state.show_length
        if st.session_state.show_length:
            st.write(f"{len(sequence)} base(s)")

    # Nucleotide Frequencies
    if st.button("Nucleotide Frequencies"):
        st.session_state.show_nucleotide_freq = not st.session_state.show_nucleotide_freq
        if st.session_state.show_nucleotide_freq:
            freqs, _ = nucleotide_frequencies(sequence)
            st.write(f"**A**: {freqs['A']}")
            st.write(f"**T**: {freqs['T']}")
            st.write(f"**C**: {freqs['C']}")
            st.write(f"**G**: {freqs['G']}")

    # GC Content
    if st.button("GC Content"):
        st.session_state.show_gc = not st.session_state.show_gc
        if st.session_state.show_gc:
            gc_content = gc_fraction(sequence) * 100  # GC content in percentage
            st.write(f"{gc_content:.2f}%")

    # Melting Temperature (Tm) Calculation
    if st.button("Tm (Wallace's Rule)"):
        st.session_state.show_tm = not st.session_state.show_tm
        if st.session_state.show_tm:
            tm_celsius = wallace_tm(sequence)
            st.write(f"{tm_celsius} °C")  # Only show in Celsius

    # Convert to cDNA (Reverse Complement)
    if st.button("cDNA (Reverse Complement)"):
        st.session_state.show_cDNA = not st.session_state.show_cDNA
        if st.session_state.show_cDNA:
            cDNA = Seq(sequence).reverse_complement()
            st.write(f"{cDNA}")

    # Transcribe to RNA
    if st.button("RNA Transcription"):
        st.session_state.show_rna = not st.session_state.show_rna
        if st.session_state.show_rna:
            rna_seq = Seq(sequence).transcribe()
            st.write(f"{rna_seq}")

    # Translate to Protein (removing spaces and newlines, and checking sequence length)
    if st.button("Protein Translation"):
        st.session_state.show_protein = not st.session_state.show_protein
        if st.session_state.show_protein:
            # Clean the sequence by removing spaces and newlines
            cleaned_sequence_for_protein = sequence.replace(" ", "").replace("\n", "")
            
            # Check if the sequence length is a multiple of 3 (for correct translation)
            seq_length = len(cleaned_sequence_for_protein)
            if seq_length % 3 != 0:
                st.warning(f"Warning: Sequence length ({seq_length}) is not a multiple of 3, translation may be incomplete.")

            try:
                # Translate the cleaned sequence
                protein = Seq(cleaned_sequence_for_protein).translate()

                # **Stop Codon Indication**: The `*` symbol represents the stop codon
                st.write("The `*` represents a stop codon, marking the end of translation.")
                # Display the protein sequence as plain text using the <pre> tag to avoid formatting issues
                # Wrap the sequence for better readability
                wrapped_protein_sequence = "\n".join([str(protein)[i:i+60] for i in range(0, len(str(protein)), 60)])  # Wrap at 60 characters per line
                st.markdown(f"<pre>{wrapped_protein_sequence}</pre>", unsafe_allow_html=True)
            except Exception as e:
                st.write(f"Error in translation: {e}")

    # Sequence Alignment (Show sequence with highlighted codon)
    codon_to_highlight = st.text_input("Sequence Alignment (e.g., ATG or atg)")
    if codon_to_highlight:
        st.markdown("**Hover to show position range**")
        highlighted_sequence, positions_text = display_alignment(sequence, codon_to_highlight)
        st.markdown(f"<p style='white-space: pre-wrap; font-family: monospace'>{highlighted_sequence}</p>", unsafe_allow_html=True)

    # Amino Acid Alignment (Show protein sequence with highlighted amino acids)
    amino_acid_to_highlight = st.text_input("Amino Acid Alignment (e.g., MAQ)")
    if amino_acid_to_highlight:  
        st.markdown("**Hover to show position range**")
        protein_sequence = Seq(sequence).translate()  # Translate the sequence
        highlighted_protein_sequence, positions_text = display_alignment(str(protein_sequence), amino_acid_to_highlight)
        st.markdown(f"<p style='white-space: pre-wrap; font-family: monospace'>{highlighted_protein_sequence}</p>", unsafe_allow_html=True)
