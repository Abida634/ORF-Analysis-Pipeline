import io

import pandas as pd
import streamlit as st

from sequence_utils import prepare_dna, read_fasta
from translator import find_longest_orf, calculate_statistics
from output_utils import CSV_COLUMNS
from main import analyze_sequence

# ---------- Page setup (must be the first Streamlit call) ----------
st.set_page_config(page_title="ORF Finder", page_icon="🧬", layout="wide")

EXAMPLE_DNA = "ATGGCCAAATTTTAACCATGAAATAGGGCTATTTCAT"

# ---------- Custom styling (CSS) ----------
# Fonts come from Google Fonts: Poppins for text, JetBrains Mono for DNA/proteins.
CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Poppins:wght@400;500;600;700&family=JetBrains+Mono:wght@400;600&display=swap');

.stApp, .stApp p, .stApp label, .stApp h1, .stApp h2, .stApp h3 {
    font-family: 'Poppins', sans-serif;
}
.stApp code, .stApp pre {
    font-family: 'JetBrains Mono', monospace !important;
}

/* Hero banner */
.hero {
    background: linear-gradient(120deg, #4f46e5 0%, #06b6d4 100%);
    border-radius: 18px;
    padding: 2rem 2.2rem;
    margin-bottom: 1.5rem;
    box-shadow: 0 8px 24px rgba(79, 70, 229, 0.25);
}
.hero h1 { color: #fff; margin: 0; font-size: 2.2rem; font-weight: 700; }
.hero p  { color: #e0f2fe; margin: 0.4rem 0 0 0; font-size: 1.05rem; }

/* Metric cards */
.card {
    background: rgba(99, 102, 241, 0.08);
    border: 1px solid rgba(99, 102, 241, 0.25);
    border-left: 5px solid #6366f1;
    border-radius: 14px;
    padding: 0.9rem 1.1rem;
}
.card .label { font-size: 0.8rem; opacity: 0.75; text-transform: uppercase; letter-spacing: 0.06em; }
.card .value { font-size: 1.8rem; font-weight: 700; color: #6366f1; line-height: 1.2; }

/* ORF map */
.track-name { font-size: 0.8rem; font-weight: 600; margin: 0.6rem 0 0.2rem 0; opacity: 0.8; }
.track { position: relative; height: 30px; border-radius: 8px; background: rgba(148, 163, 184, 0.2); }
.orf-bar {
    position: absolute; top: 4px; height: 22px; border-radius: 6px;
    min-width: 4px; opacity: 0.9;
}
.fwd { background: #6366f1; }
.rev { background: #06b6d4; }
.scale { display: flex; justify-content: space-between; font-size: 0.75rem; opacity: 0.6; margin-top: 0.3rem; }

/* Sidebar + buttons */
section[data-testid="stSidebar"] { border-right: 1px solid rgba(99, 102, 241, 0.2); }
.stDownloadButton button {
    background: linear-gradient(120deg, #4f46e5, #06b6d4);
    color: white; border: none; border-radius: 10px; font-weight: 600;
}
.footer { text-align: center; opacity: 0.55; font-size: 0.8rem; margin-top: 2.5rem; }
</style>
"""
st.markdown(CSS, unsafe_allow_html=True)


# ---------- Small helper functions that build HTML ----------
def metric_card(label, value):
    return f'<div class="card"><div class="label">{label}</div><div class="value">{value}</div></div>'


def orf_map_html(orfs, sequence_length):
    """Draw two tracks (forward / reverse strand) with one bar per ORF."""
    tracks = ""
    for strand, name, css_class in [("+", "Forward strand (+)", "fwd"),
                                    ("-", "Reverse strand (-)", "rev")]:
        bars = ""
        for orf in orfs:
            if orf["strand"] != strand:
                continue
            left = (orf["start"] - 1) / sequence_length * 100
            width = (orf["end"] - orf["start"] + 1) / sequence_length * 100
            tip = f"{orf['orf_id']}: {orf['start']}-{orf['end']} ({orf['length']} nt)"
            bars += f'<div class="orf-bar {css_class}" style="left:{left}%; width:{width}%;" title="{tip}"></div>'
        tracks += f'<div class="track-name">{name}</div><div class="track">{bars}</div>'
    scale = f'<div class="scale"><span>1</span><span>{sequence_length} bp</span></div>'
    return tracks + scale


# ---------- Hero banner ----------
st.markdown(
    """
    <div class="hero">
        <h1>🧬 ORF Finder</h1>
        <p>Find candidate open reading frames in all six reading frames, translate them, and download the results.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# ---------- Sidebar: all settings in one place ----------
with st.sidebar:
    st.header("⚙️ Settings")
    input_type = st.radio("Input type", ["Paste DNA", "Upload FASTA file"])
    min_length = st.slider("Minimum ORF length (nt)", 9, 300, 9, step=3,
                           help="Shorter ORFs are filtered out. Real gene finding often uses 300 nt or more.")
    use_example = st.checkbox("Use example sequence")
    st.caption("ORF = ATG ... first in-frame stop (TAA, TAG, TGA).")

# ---------- Read the input ----------
sequences = []   # list of (sequence_id, clean_dna)

try:
    if input_type == "Paste DNA":
        if use_example:
            st.info("Using the built-in example sequence.")
            st.code(EXAMPLE_DNA, language=None)
            raw_dna = EXAMPLE_DNA
        else:
            raw_dna = st.text_area("DNA sequence", height=140,
                                   placeholder="Paste A, T, G, C letters here...")
        if raw_dna.strip() != "":
            sequences = [("manual_input", prepare_dna(raw_dna))]
    else:
        uploaded_file = st.file_uploader("FASTA file", type=["fasta", "fa", "fna", "txt"])
        if uploaded_file is not None:
            fasta_text = uploaded_file.getvalue().decode("utf-8")
            sequences = read_fasta(io.StringIO(fasta_text))
except ValueError as error:
    st.error(f"⚠️ Input problem: {error}")
    st.stop()

if len(sequences) == 0:
    st.info("👈 Paste a DNA sequence, upload a FASTA file, or tick 'Use example sequence' to begin.")
    st.stop()

# ---------- Analysis and display ----------
all_orfs = []

for sequence_id, dna in sequences:
    st.subheader(f"🔹 {sequence_id}  ·  {len(dna)} bases")
    orfs = analyze_sequence(sequence_id, dna, min_length)

    if len(orfs) == 0:
        st.warning(f"No ORFs of at least {min_length} nt found.")
        continue

    stats = calculate_statistics(orfs)
    longest = find_longest_orf(orfs)

    cols = st.columns(4)
    cols[0].markdown(metric_card("Total ORFs", stats["total_orfs"]), unsafe_allow_html=True)
    cols[1].markdown(metric_card("Forward / Reverse",
                                 f"{stats['forward_orfs']} / {stats['reverse_orfs']}"), unsafe_allow_html=True)
    cols[2].markdown(metric_card("Longest (nt)", longest["length"]), unsafe_allow_html=True)
    cols[3].markdown(metric_card("Average (nt)", stats["average_length"]), unsafe_allow_html=True)
    st.write("")

    tab_table, tab_map, tab_proteins = st.tabs(["📋 Table", "🗺️ ORF map", "🔬 Proteins"])

    with tab_table:
        table = pd.DataFrame(orfs)[CSV_COLUMNS]
        st.dataframe(table, hide_index=True, width="stretch")

    with tab_map:
        st.markdown(orf_map_html(orfs, len(dna)), unsafe_allow_html=True)
        st.caption("Hover over a bar to see its ORF. Coordinates are always on the forward strand.")

    with tab_proteins:
        for orf in orfs:
            title = f"{orf['orf_id']}  ·  {orf['strand']}{orf['frame']}  ·  {orf['protein_length']} aa"
            with st.expander(title):
                st.code(orf["protein"], language=None)

    all_orfs.extend(orfs)

# ---------- Download ----------
if len(all_orfs) > 0:
    st.divider()
    csv_text = pd.DataFrame(all_orfs)[CSV_COLUMNS].to_csv(index=False)
    st.download_button("⬇️ Download results as CSV", data=csv_text,
                       file_name="orf_results.csv", mime="text/csv")

st.markdown('<div class="footer">Candidate ORFs only, not confirmed genes · Built with Python, Biopython & Streamlit</div>',
            unsafe_allow_html=True)