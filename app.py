"""
app.py — entrypoint
--------------------
Declares the app's pages with clean, explicit titles, icons and URL slugs
via st.Page / st.navigation, instead of relying on filenames to derive
them. Page content lives in views/home.py and views/independent_tools.py.
"""

import streamlit as st
import config

st.set_page_config(
    page_title=config.PAGE_TITLE,
    page_icon=config.PAGE_ICON,
    layout=config.DEFAULT_LAYOUT,
    initial_sidebar_state=config.DEFAULT_SIDEBAR_STATE,
)

home_page = st.Page(
    "views/home.py",
    title="Analyse de séquence",
    icon="🧬",
    url_path="analyse",
    default=True,
)

independent_tools_page = st.Page(
    "views/independent_tools.py",
    title="Outils indépendants",
    icon="🧪",
    url_path="outils-independants",
)

pg = st.navigation([home_page, independent_tools_page])
pg.run()