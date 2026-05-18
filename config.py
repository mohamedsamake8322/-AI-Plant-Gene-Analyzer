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

# Create directories if they don't exist
LOG_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True)

# ─── Logging Configuration ────────────────────────────────────────────────────
LOG_LEVEL = logging.INFO
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_FILE = LOG_DIR / "analyzer.log"

# ─── Bioinformatics Parameters ────────────────────────────────────────────────

# Sequence validation
MIN_SEQUENCE_LENGTH = 10
MAX_SEQUENCE_LENGTH = 1_000_000  # 1 million bp

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
