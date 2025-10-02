import streamlit as st
from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.Data import CodonTable
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner
import matplotlib.pyplot as plt
import re

# ----- App config -----
st.set_page_config(page_title="GeneScope", layout="wide")
st.title("GeneScope")
st.write("Simple Biopython pipeline — fetch an accession or upload a FASTA. All sequence ops use Biopython.")

# ----- Sidebar -----
st.sidebar.header("Settings")
email = st.sidebar.text_input("NCBI email (required for Entrez/BLAST)", value="")
if email:
    Entrez.email = email

st.sidebar.write("Optional outputs (toggle on to compute)")
show_orfs = st.sidebar.checkbox("Show ORFs (6-frame, M->stop)", value=False)
show_codon = st.sidebar.checkbox("Show codon usage (frame 1)", value=False)
show_cds = st.sidebar.checkbox("Show annotated CDS (GenBank)", value=False)

BLAST_MAX_HITS = 5
ORF_MIN_AA = 30

# ----- Input -----
st.header("Input")
col1, col2 = st.columns([3,2])
with col1:
    acc = st.text_input("NCBI accession (e.g. NC_005816)")
    fetch_btn = st.button("Fetch from NCBI")
with col2:
    uploaded = st.file_uploader("Or upload FASTA", type=["fa","fasta","txt"])
    load_btn = st.button("Load uploaded FASTA")

if "rec" not in st.session_state:
    st.session_state.rec = None

# ----- Helper functions (short and clear) -----
def fetch_from_ncbi(accession):
    if not email:
        raise ValueError("Provide email in sidebar first.")
    # try genbank to get annotations; fallback to fasta
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as h:
            return SeqIO.read(h, "genbank")
    except Exception:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as h:
            return SeqIO.read(h, "fasta")

def seq_stats(rec):
    s = str(rec.seq).upper()
    L = len(s)
    gc = 100.0 * (s.count("G")+s.count("C"))/L if L else 0.0
    at = 100.0 * (s.count("A")+s.count("T"))/L if L else 0.0
    amb = sum(1 for ch in s if ch not in "ACGT")
    amb_pct = 100.0 * amb / L if L else 0.0
    counts = {b: s.count(b) for b in ("A","C","G","T","N")}
    mw = round(molecular_weight(s, seq_type="DNA"),2) if L else 0.0
    return {"id": rec.id, "len": L, "GC%": round(gc,3), "AT%": round(at,3), "amb%": round(amb_pct,3), "mw": mw, **counts}

def transcribe_translate(rec):
    rna = rec.seq.transcribe()
    prot = rec.seq.translate()
    return rna, prot

def find_orfs(rec, min_aa=ORF_MIN_AA):
    s = str(rec.seq).upper()
    L = len(s)
    out = []
    # forward
    for frame in range(3):
        prot = Seq(s[frame:]).translate()
        aa = str(prot)
        start = None
        for i,ch in enumerate(aa):
            if ch=="M" and start is None:
                start = i
            if ch=="*" and start is not None:
                end = i-1
                length = end - start + 1
                if length >= min_aa:
                    nt_start = frame + start*3 + 1
                    nt_end = frame + end*3 + 3
                    out.append(("+", frame, length, nt_start, nt_end))
                start = None
    # reverse
    rc = str(Seq(s).reverse_complement())
    for frame in range(3):
        prot = Seq(rc[frame:]).translate()
        aa = str(prot)
        start = None
        for i,ch in enumerate(aa):
            if ch=="M" and start is None:
                start = i
            if ch=="*" and start is not None:
                end = i-1
                length = end - start + 1
                if length >= min_aa:
                    nt_end_rc = frame + start*3 + 1
                    nt_start_rc = frame + end*3 + 3
                    nt_start = L - nt_start_rc + 1
                    nt_end = L - nt_end_rc + 1
                    out.append(("-", frame, length, nt_start, nt_end))
                start = None
    return out

def codon_usage(rec):
    s = str(rec.seq).upper()
    codons = [s[i:i+3] for i in range(0,len(s)-2,3) if len(s[i:i+3])==3]
    cnt = {}
    for c in codons:
        cnt[c] = cnt.get(c,0)+1
    table = CodonTable.unambiguous_dna_by_id[1]
    rows = [{"codon":k,"count":v,"aa":("*" if k in table.stop_codons else table.forward_table.get(k,"X"))}
            for k,v in sorted(cnt.items(), key=lambda x:-x[1])]
    return rows

def extract_cds(rec):
    cds = []
    try:
        for f in rec.features:
            if f.type=="CDS":
                quals = {k:(";".join(v) if isinstance(v,list) else str(v)) for k,v in f.qualifiers.items()}
                cds.append({"loc":str(f.location), **quals})
    except Exception:
        pass
    return cds

def run_blastn_top(seq_text, top=BLAST_MAX_HITS):
    res = NCBIWWW.qblast("blastn","nt", seq_text)
    hits = []
    for r in NCBIXML.parse(res):
        for a in r.alignments[:top]:
            h = a.hsps[0] if a.hsps else None
            hits.append({"title":a.title, "len":a.length, "score": h.score if h else None, "e": h.expect if h else None})
        break
    res.close()
    return hits

def find_pdb(prot_seq):
    res = NCBIWWW.qblast("blastp","pdb", str(prot_seq))
    pdb = None
    summary = []
    for r in NCBIXML.parse(res):
        for a in r.alignments[:1]:
            m = re.search(r"\b([0-9][A-Za-z0-9]{3})\b", a.title)
            if m:
                pdb = m.group(1).upper()
            h = a.hsps[0] if a.hsps else None
            summary.append({"title": a.title, "len": a.length, "score": h.score if h else None, "e": h.expect if h else None})
        break
    res.close()
    return pdb, summary

# ----- Load sequence -----
loaded = None
if fetch_btn and acc:
    if not email:
        st.error("Enter email in the sidebar before fetching.")
    else:
        try:
            st.session_state.rec = fetch_from_ncbi(acc)
            loaded = st.session_state.rec
            st.success(f"Fetched {loaded.id} ({len(loaded.seq)} nt)")
        except Exception as e:
            st.error(f"Fetch error: {e}")
elif load_btn and uploaded:
    try:
        text = uploaded.getvalue().decode("utf-8")
        rec = SeqIO.read(StringIO(text), "fasta")
        st.session_state.rec = rec
        loaded = rec
        st.success(f"Loaded {rec.id} from upload")
    except Exception as e:
        st.error(f"Upload error: {e}")
elif st.session_state.rec:
    loaded = st.session_state.rec

# ----- Main output -----
if loaded:
    rec = loaded
    stats = seq_stats(rec)

    st.header("Summary")
    c1,c2,c3,c4,c5 = st.columns(5)
    c1.metric("Length (nt)", stats["len"])
    c2.metric("GC %", f"{stats['GC%']}")
    c3.metric("AT %", f"{stats['AT%']}")
    c4.metric("Ambiguous %", f"{stats['amb%']}")
    c5.metric("Mol weight (Da)", stats["mw"])

    st.subheader("Base counts")
    cols = st.columns(5)
    for i,(base,val) in enumerate([("A",stats["A"]),("C",stats["C"]),("G",stats["G"]),("T",stats["T"]),("N",stats["N"])]):
        cols[i].metric(base, val)

    st.subheader("FASTA")
    s_io = StringIO()
    SeqIO.write(rec, s_io, "fasta")
    st.code(s_io.getvalue())

    # full transcription & translation (always shown, no truncation)
    st.subheader("Transcription & Translation")
    rna, prot = transcribe_translate(rec)
    st.markdown("**RNA**")
    st.code(str(rna))
    st.markdown("**Protein**")
    st.code(str(prot))

    # meaningful plots
    st.subheader("Plots")
    # base composition bar chart
    bases = ["A","C","G","T"]
    counts = [stats[b] for b in bases]
    fig1, ax1 = plt.subplots(figsize=(4,3))
    ax1.bar(bases, counts)
    ax1.set_title("Base composition")
    ax1.set_ylabel("Count")
    st.pyplot(fig1)

    # ORF histogram if enabled
    if show_orfs:
        orfs = find_orfs(rec)
        if orfs:
            lengths = [o[2] for o in orfs]  # aa lengths
            fig2, ax2 = plt.subplots(figsize=(4,3))
            ax2.hist(lengths, bins=10)
            ax2.set_title("ORF length distribution (aa)")
            ax2.set_xlabel("AA length")
            ax2.set_ylabel("Count")
            st.pyplot(fig2)
            st.write("ORFs (strand,frame,len,nt_start,nt_end)")
            st.table([{"strand":o[0],"frame":o[1],"aa_len":o[2],"nt_start":o[3],"nt_end":o[4]} for o in orfs])
        else:
            st.info("No ORFs found with current settings.")

    # codon usage plot if enabled
    if show_codon:
        codons = codon_usage(rec)
        if codons:
            top = codons[:15]
            labels = [r["codon"] for r in top]
            vals = [r["count"] for r in top]
            fig3, ax3 = plt.subplots(figsize=(6,3))
            ax3.bar(labels, vals)
            ax3.set_title("Top codons (frame 1)")
            ax3.set_ylabel("Count")
            plt.xticks(rotation=45)
            st.pyplot(fig3)
            st.table(top)
        else:
            st.info("Sequence too short for codon counts.")

    # CDS table if enabled
    if show_cds:
        cds = extract_cds(rec)
        if cds:
            st.subheader("Annotated CDS (from GenBank)")
            st.table(cds)
        else:
            st.info("No CDS features found in GenBank record.")

    # pairwise alignment section
    st.subheader("Pairwise alignment (optional)")
    align_up = st.file_uploader("Upload FASTA with targets for pairwise alignment", type=["fa","fasta","txt"], key="align2")
    if st.button("Run pairwise alignment"):
        if not align_up:
            st.error("Upload target FASTA first.")
        else:
            targets = list(SeqIO.parse(StringIO(align_up.getvalue().decode("utf-8")), "fasta"))
            if not targets:
                st.error("No sequences found in uploaded file.")
            else:
                aligner = PairwiseAligner()
                rows = []
                for t in targets:
                    score = aligner.score(rec.seq, t.seq)
                    rows.append({"target": t.id, "len": len(t.seq), "score": score})
                rows = sorted(rows, key=lambda x:-x["score"])
                st.table(rows)
                # preview top alignment
                top = targets[0]
                aln = aligner.align(rec.seq, top.seq)[0]
                st.code(str(aln))

    # PDB lookup (BLASTp -> PDB) and show viewer embedded if found
    st.subheader("Protein structure lookup")
    if st.button("Find PDB (BLASTp top hit)"):
        if not email:
            st.error("Enter email in sidebar for BLASTp.")
        else:
            # prefer CDS translation if present
            protein_seq = None
            try:
                for f in rec.features:
                    if f.type=="CDS":
                        if "translation" in f.qualifiers:
                            protein_seq = f.qualifiers["translation"][0]
                            break
                        else:
                            try:
                                cds_seq = f.location.extract(rec).seq
                                protein_seq = cds_seq.translate()
                                break
                            except Exception:
                                continue
            except Exception:
                protein_seq = None
            if protein_seq is None:
                protein_seq = rec.seq.translate()
            with st.spinner("Running BLASTp against PDB..."):
                try:
                    pdb_id, summ = find_pdb(protein_seq)
                    if summ:
                        st.table(summ)
                    else:
                        st.info("No BLASTp hits returned.")
                    if pdb_id:
                        st.success(f"PDB match: {pdb_id}")
                        # Try to embed Mol* viewer via RCSB. If embedding blocked, show link.
                        viewer = f"https://www.rcsb.org/3d-view/entry/{pdb_id}"
                        iframe = f'<iframe src="{viewer}" width="100%" height="650" frameborder="0"></iframe>'
                        try:
                            st.components.v1.html(iframe, height=660, scrolling=True)
                        except Exception:
                            st.markdown(f"Open PDB viewer: {viewer}")
                    else:
                        st.info("No PDB id could be extracted from BLASTp top hit.")
                except Exception as e:
                    st.error(f"PDB lookup failed: {e}")

    # BLASTn top 5
    st.subheader("BLASTn (top 5 hits)")
    if st.button("Run BLASTn (top 5)"):
        if not email:
            st.error("Enter email in sidebar for BLAST.")
        else:
            with st.spinner("Running BLASTn..."):
                try:
                    hits = run_blastn_top(str(rec.seq), top=BLAST_MAX_HITS)
                    if hits:
                        st.table(hits)
                    else:
                        st.info("No BLASTn hits.")
                except Exception as e:
                    st.error(f"BLASTn failed: {e}")

else:
    st.info("No sequence loaded yet. Provide an accession or upload a FASTA and click Fetch/Load.")