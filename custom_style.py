import streamlit as st

def inject_custom_css():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&family=IBM+Plex+Mono:wght@400;500&display=swap');

    :root {
        --sage: #7A8B5C;
        --sage-dim: #5f6e47;
        --amber: #B8873B;
        --border: #2E3323;
        --text-dim: #9BA08D;
        --mono-bg: #10130C;
    }

    /* Titres en serif */
    h1, h2, h3 {
        font-family: 'Source Serif 4', serif !important;
        font-weight: 600 !important;
        letter-spacing: 0.2px;
    }

    /* Enlever le fond à pois / motif décoratif éventuel sur le bloc principal */
    [data-testid="stAppViewContainer"] {
        background-image: none !important;
    }

    /* Boutons : couleur pleine, pas de dégradé */
    .stButton > button, button[kind="primary"] {
        background: var(--sage) !important;
        color: #12150C !important;
        border: none !important;
        border-radius: 6px !important;
        font-weight: 600 !important;
        transition: background 0.15s ease;
    }
    .stButton > button:hover, button[kind="primary"]:hover {
        background: #8b9d69 !important;
        color: #12150C !important;
    }

    /* Zones de texte / séquences en mono */
    textarea, .stTextArea textarea {
        font-family: 'IBM Plex Mono', monospace !important;
        background: var(--mono-bg) !important;
        border: 1px solid var(--border) !important;
        color: var(--text-dim) !important;
    }
    textarea:focus {
        border-color: var(--sage-dim) !important;
        box-shadow: none !important;
    }

    /* Onglets : soulignement amber sur l'onglet actif au lieu du rouge/rose par défaut */
    .stTabs [aria-selected="true"] {
        color: #EDEAE0 !important;
        border-bottom-color: var(--amber) !important;
    }

    /* Sidebar : bordure discrète pour la séparer du contenu */
    [data-testid="stSidebar"] {
        border-right: 1px solid var(--border);
    }

    /* Slider : accent sauge au lieu du rouge par défaut */
    [data-testid="stSlider"] div[role="slider"] {
        background-color: var(--sage) !important;
        border-color: var(--sage) !important;
    }
    </style>
    """, unsafe_allow_html=True)
