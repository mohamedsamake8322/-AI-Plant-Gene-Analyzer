"""
config.py
---------
Centralized configuration for the AI-Powered Plant Gene Analyzer.
Manages all constants, thresholds, and settings.
"""

import logging
from pathlib import Path

# ─── Project Paths ────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent
DATABASE_PATH = PROJECT_ROOT / "genes_database.json"
STYLE_CSS_PATH = PROJECT_ROOT / "style.css"
LOG_DIR = PROJECT_ROOT / "logs"
RESULTS_DIR = PROJECT_ROOT / "results"

# Create directories if they don't exist. Guarded because config.py is
# imported at module load time by every other module; a read-only deployment
# target (e.g. a locked-down container) would otherwise crash the whole app
# on import with an unhandled PermissionError instead of failing gracefully
# where logging/export actually happens.
for _dir in (LOG_DIR, RESULTS_DIR):
    try:
        _dir.mkdir(exist_ok=True)
    except OSError:
        pass

# ─── Logging Configuration ────────────────────────────────────────────────────
LOG_LEVEL = logging.INFO
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_FILE = LOG_DIR / "analyzer.log"

# ─── Bioinformatics Parameters ────────────────────────────────────────────────

# Sequence validation
MIN_SEQUENCE_LENGTH = 10  # bp, for DNA/RNA input
MIN_PROTEIN_LENGTH = 5  # aa, for protein input
MAX_SEQUENCE_LENGTH = 1_000_000  # 1 million bp — hard cap for O(n) stats (GC%, motifs, translation)

# Length-ratio prefilter for similarity search: a global alignment between
# sequences differing in length by more than this factor is dominated by
# forced gaps and essentially never ranks highly, so candidates outside this
# ratio are skipped before running a full O(n*m) alignment against them.
# Shared by similarityengine.compare_with_database (in-Python prefilter) and
# postgres_utils.find_candidate_genes_by_kmer (server-side candidate lookup)
# so both stay in sync -- previously each hardcoded its own copy of "3.0".
LENGTH_RATIO_PREFILTER = 3.0
# similarity search and mutation detection is O(n*m) in both time AND
# memory (two full DP matrices), unlike the O(n) basic statistics above.
# With int32 score + int8 traceback matrices (see alignment_engine.py),
# two 5,000x5,000 matrices are ~125 MB combined -- reasonable for a shared,
# publicly hosted instance. A larger cap here is not just slower but risks
# exhausting memory outright (e.g. 20,000x20,000 would be ~1.6 GB). This
# covers the vast majority of real single-gene/transcript-length plant
# sequences; true long-sequence alignment (whole loci, chromosomes) needs a
# different algorithm (BLAST-style seed-and-extend heuristic, already on
# the project roadmap) rather than a bigger cap on this exact DP.
MAX_ALIGNMENT_SEQUENCE_LENGTH = 5_000  # bp / aa

# Uploaded file size cap (bytes), checked before reading the file into memory.
MAX_UPLOAD_SIZE_BYTES = 5 * 1024 * 1024  # 5 MB

# GC content thresholds
GC_HIGH = 60.0
GC_LOW = 35.0
GC_OPTIMAL = 50.0

# Similarity thresholds
SIMILARITY_VERY_HIGH = 90.0
SIMILARITY_HIGH = 75.0
SIMILARITY_MODERATE = 55.0
SIMILARITY_LOW = 35.0

# Mutation thresholds
MUTATION_RATE_HIGH = 10.0
MUTATION_RATE_LOW = 2.0

# Window size for sliding window similarity
DEFAULT_WINDOW_SIZE = 30
MIN_WINDOW_SIZE = 5
MAX_WINDOW_SIZE = 60

# Reading frames
READING_FRAMES = [0, 1, 2]

# ─── Database Parameters ──────────────────────────────────────────────────────
DEFAULT_TOP_N_MATCHES = 3
MAX_TOP_N_MATCHES = 8

# ─── UI/UX Settings ────────────────────────────────────────────────────────────
PAGE_TITLE = "AI Plant Gene Analyzer"
PAGE_ICON = "🧬"
DEFAULT_LAYOUT = "wide"
DEFAULT_SIDEBAR_STATE = "expanded"

# Color palette for nucleotides
NUCLEOTIDE_COLORS = {
    "A": "#00c853",  # green
    "T": "#2979ff",  # blue
    "G": "#ffd600",  # yellow
    "C": "#ff3d00",  # red
    "N": "#9e9e9e",  # grey
}

# Light theme colors
CHART_BG = "#ffffff"
CHART_PAPER = "#ffffff"
CHART_FONT_COLOR = "#111111"
CHART_GRID_COLOR = "rgba(0,0,0,0.08)"

# ─── Export Settings ──────────────────────────────────────────────────────────
EXPORT_FORMATS = ["JSON", "CSV", "XLSX", "HTML"]
DEFAULT_EXPORT_FORMAT = "JSON"

# ─── Sequence input settings ────────────────────────────────────────────────
SUPPORTED_INPUT_TYPES = ["Auto detect", "DNA", "Protein"]
MAX_INPUT_SEQUENCES = 20

# ─── Demo Sequences ───────────────────────────────────────────────────────────
DEMO_SEQUENCES = {
    "Select a demo…": {"seq": "", "desc": ""},
    "DREB-like (Drought Resistance)": {
        "seq": "ATGCGTAGCTAGCGATCGATCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCG",
        "desc": "Dehydration-Responsive Element Binding protein fragment",
    },
    "HSP-like (Heat Stress)": {
        "seq": "GGCTAACCGTAGCTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT",
        "desc": "Heat Shock Protein gene fragment",
    },
    "RbcL (Photosynthesis)": {
        "seq": "ATGGCTTCAGAAACAGGTATTGATTTAGAGAAAGAAATGAAAGGTTATAAATTAGGTGGT",
        "desc": "RuBisCO large subunit — Calvin cycle enzyme",
    },
    "Novel / Unknown sequence": {
        "seq": "AATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT",
        "desc": "AT-rich sequence with low database similarity",
    },
}

# ─── Utility function to get logger ────────────────────────────────────────────
def get_logger(name: str) -> logging.Logger:
    """Get or create a configured logger."""
    logger = logging.getLogger(name)
    
    if not logger.handlers:
        # Create handlers
        file_handler = logging.FileHandler(LOG_FILE)
        console_handler = logging.StreamHandler()
        
        # Set level
        file_handler.setLevel(LOG_LEVEL)
        console_handler.setLevel(LOG_LEVEL)
        
        # Create formatter
        formatter = logging.Formatter(LOG_FORMAT)
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)
        
        # Add handlers
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        logger.setLevel(LOG_LEVEL)
    
    return logger