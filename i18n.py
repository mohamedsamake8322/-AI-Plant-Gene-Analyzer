import streamlit as st

DEFAULT_LANG = "fr"
SUPPORTED_LANGS = ["fr", "en", "tr"]

TRANSLATIONS = {
    "fr": {
        "ui": {
            "app_title": "Analyseur de gènes végétaux",
            "about": "À propos",
            "settings": "Paramètres",
            "database": "Base de données",
            "language": "Langue",
            "sequence_analysis": "Analyse de séquence",
            "independent_tools": "Outils indépendants",
            "sequence_input": "Entrée de séquence",
            "quick_demo": "Démo rapide",
            "load_demo_sequence": "Charger la séquence de démonstration",
            "top_matches": "Meilleurs résultats de la base à afficher",
            "top_matches_help": "Nombre de gènes les plus proches à afficher.",
            "deep_search": "Activer la recherche de similarité approfondie",
            "deep_search_help": "Désactive le préfiltre de longueur d’alignement et évalue plus de candidats. C’est plus lent, mais augmente la sensibilité pour les séquences courtes ou divergentes.",
            "window_size": "Fenêtre glissante (profil GC)",
            "window_size_help": "Taille de la fenêtre (pb) pour le graphique du profil GC.",
            "reading_frame": "Cadre de lecture pour la traduction",
            "input_type": "Type d’entrée",
            "input_type_help": "Choisissez le type de séquence ou laissez l’application le détecter automatiquement.",
            "input_accepts": "Cette application acceptera **{type}** comme entrée et ajustera l’analyse en conséquence.",
            "auto_detect": "Détection automatique",
            "dna": "ADN",
            "protein": "Protéine",
            "file_support_help": "Les fichiers FASTA commençant par des lignes d’en-tête '>' sont pris en charge.",
            "sequence_placeholder": "Collez une séquence ADN, protéique ou FASTA ici.\nExemple : ATGCGTAGCTAGCGATCGATCGAATTCG...",
            "file_loaded": "Fichier chargé :",
            "selected_sequence": "Sélectionné :",
            "batch_summary_text": "{count} séquences seront analysées en lot.",
            "search_gene": "Rechercher un ID, symbole ou trait de gène",
            "preview_metadata": "Aperçu des métadonnées des gènes",
            "showing_sample": "Affichage d’un échantillon de 20 enregistrements de métadonnées. Utilisez la recherche pour filtrer les gènes.",
            "full_db_load": "La base complète des gènes avec séquences sera chargée au lancement d’une analyse.",
            "build_about": "Analyse les séquences ADN végétales pour :",
            "standalone_tools_message": "Outils autonomes qui fonctionnent sur n’importe quelle séquence que vous collez ici — aucune analyse principale n’est nécessaire.",
            "about_items": [
                "Teneur en GC et statistiques nucléotidiques",
                "Similarité avec la base de données de gènes",
                "Détection des mutations",
                "Interprétation biologique par IA",
                "Recommandations agricoles",
            ],
            "upload_file": "Téléverser un fichier de séquence (.fasta / .fa / .txt)",
            "paste_sequence": "Ou collez votre séquence ici :",
            "select_sequence": "Choisir une séquence à analyser",
            "analyze_all": "Analyser toutes les séquences dans cette entrée FASTA",
            "file_loaded": "Fichier chargé :",
            "choose_demo": "Charger une séquence de démonstration",
            "reading_frame_preview": "Aperçu du cadre de lecture",
            "reading_frame_caption": "Résumé en six cadres pour guider le choix du cadre de lecture avant l’analyse.",
            "analyze_button": "Analyser la séquence",
            "batch_complete": "L’analyse par lot est terminée",
            "batch_summary": "Résumé du lot",
            "select_sequence_inspect": "Sélectionner une séquence à inspecter dans ce lot :",
            "independent_tools_title": "Outils d’analyse indépendante",
            "independent_tools_subtitle": "Alignements · Matrice de distance · Phylogénie · Analyse protéique",
            "tab_alignments": "Alignements",
            "tab_distance_matrix": "Matrice de distance",
            "tab_phylogeny": "Phylogénie",
            "tab_protein_analysis": "Analyse protéique",
            "tab_trait_search": "Recherche par thème",
            "hero_subtitle": "Bioinformatique · Interprétation IA · Perspectives agronomiques",
            "metadata_load_help": "L’application charge d’abord les métadonnées légères des gènes pour la recherche et le filtrage. Les séquences complètes ne sont chargées qu’au lancement d’une analyse.",
            "metadata_filter_help": "Filtrer la base de gènes chargée par gene_id, symbole ou trait.",
            "metadata_preview_note": "Affichage d’un échantillon de 20 enregistrements de métadonnées. Utilisez la recherche pour filtrer les gènes.",
            "metadata_load_error": "Impossible de charger l’aperçu des métadonnées des gènes. Retour à la charge complète de la base.",
            "no_genes_available": "Aucun gène disponible : configurez la connexion PostgreSQL ou ajoutez la base de secours.",
            "metadata_count_summary": "Affichage de {count} enregistrements de métadonnées correspondant à la recherche",
            "welcome_title": "Prêt à analyser",
            "welcome_message": "Collez une séquence ADN végétale ci-dessus ou chargez une démo, puis cliquez sur <b>Analyser la séquence</b>.",
            "what_does_app_analyze": "Que peut analyser cette application ?",
            "supported_input_formats": "Formats d’entrée pris en charge",
            "independent_tools_notice": "Vous cherchez des alignements, des matrices de distance, de la phylogénie ou une analyse protéique autonome ? Ils sont maintenant sur une page dédiée — voir **Outils indépendants** dans la barre latérale.",
            "msa_input_hint": "Collez plusieurs séquences FASTA ou une par ligne :",
            "run_msa": "Lancer l’alignement multiple",
            "sequence_1": "Séquence 1",
            "sequence_2": "Séquence 2",
            "align_pairwise": "Aligner en paire",
            "pairwise_missing": "Fournissez deux séquences pour l’alignement par paires.",
            "distance_matrix_title": "Calculer une matrice de distance par paires",
            "distance_input_hint": "Collez des séquences FASTA ou une séquence par ligne :",
            "method": "Méthode",
            "compute_distance_matrix": "Calculer la matrice de distance",
            "phylogeny_title": "Construire un arbre phylogénétique",
            "phylogeny_input_hint": "Collez des séquences pour la phylogénie (FASTA ou lignes) :",
            "tree_algorithm": "Algorithme d’arbre",
            "build_tree": "Construire l’arbre",
            "protein_biochemical_analysis": "Analyse biochimique des protéines",
            "paste_protein_sequence": "Collez une séquence protéique :",
            "analyze_protein": "Analyser la protéine",
            "translation_label": "Traduction",
            "language_selector_label": "Language / Langue / Dil",
            "menu": "Menu",
            "title_plant_gene_analyzer": "Analyseur de gènes végétaux",
            "sequence_analysis_tab": "Analyse de séquence",
        },
        "results": {
            "gene_records_available": "enregistrements de métadonnées disponibles",
            "genes_loaded": "gènes chargés",
            "sequence_length": "longueur de séquence",
            "gc_content": "teneur en GC",
            "best_match": "meilleur correspondant",
            "similarity_score": "score de similarité",
            "mutations": "mutations",
            "translation": "traduction",
            "protein_analysis": "analyse protéique",
        },
        "errors": {
            "no_sequence": "Veuillez saisir ou coller une séquence avant l’analyse.",
            "db_connection": "Impossible de se connecter à la base de données des gènes.",
            "invalid_file": "Fichier invalide ou encodage non pris en charge.",
            "database_missing": "Aucun gène disponible : configurez la connexion PostgreSQL ou ajoutez la base de secours.",
        },
        "glossary": {
            "gene": "gène",
            "sequence": "séquence",
            "similarity": "similarité",
            "mutation": "mutation",
            "database": "base de données",
            "metadata": "métadonnées",
            "translation": "traduction",
            "reading_frame": "cadre de lecture",
            "gc_content": "teneur en GC",
            "alignment": "alignement",
            "trait": "trait",
            "species": "espèce",
            "annotation": "annotation",
            "expression": "expression",
            "protein": "protéine",
            "organism": "organisme",
            "nucleotide": "nucléotide",
        },
    },
    "en": {
        "ui": {
            "app_title": "Plant Gene Analyzer",
            "about": "About",
            "settings": "Settings",
            "database": "Database",
            "language": "Language",
            "sequence_analysis": "Sequence Analysis",
            "independent_tools": "Independent Tools",
            "sequence_input": "Sequence Input",
            "quick_demo": "Quick Demo",
            "load_demo_sequence": "Load Demo Sequence",
            "top_matches": "Top database matches to show",
            "top_matches_help": "Number of best-matching genes to display.",
            "deep_search": "Enable deep similarity search",
            "deep_search_help": "Disable the alignment length prefilter and evaluate more candidates. This is slower, but increases sensitivity for short or divergent queries.",
            "window_size": "Sliding window (GC profile)",
            "window_size_help": "Window size (bp) for the GC content profile chart.",
            "reading_frame": "Reading frame for translation",
            "input_type": "Input type",
            "input_type_help": "Choose the sequence type or let the app detect it automatically.",
            "input_accepts": "This app will accept **{type}** input and adjust analysis accordingly.",
            "auto_detect": "Auto detect",
            "dna": "DNA",
            "protein": "Protein",
            "file_support_help": "FASTA files starting with '>' header lines are supported.",
            "sequence_placeholder": "Paste raw DNA, protein sequence, or FASTA format here.\nExample: ATGCGTAGCTAGCGATCGATCGAATTCG...",
            "file_loaded": "File loaded:",
            "selected_sequence": "Selected:",
            "batch_summary_text": "{count} sequences will be analyzed as a batch.",
            "search_gene": "Search gene ID, symbol, or trait",
            "preview_metadata": "Preview gene metadata",
            "showing_sample": "Showing a sample of 20 gene metadata records. Use search to filter specific genes.",
            "full_db_load": "Full gene database with sequences will be loaded when you start an analysis.",
            "build_about": "Analyze plant DNA sequences for:",
            "standalone_tools_message": "Standalone tools that work on any sequences you paste here — no need to run the main analysis first.",
            "about_items": [
                "GC content & nucleotide stats",
                "Gene database similarity",
                "Mutation detection",
                "AI biological interpretation",
                "Agricultural recommendations",
            ],
            "upload_file": "Upload a sequence file (.fasta / .fa / .txt)",
            "paste_sequence": "Or paste your sequence here:",
            "select_sequence": "Choose a sequence to analyze",
            "analyze_all": "Analyze all sequences in this FASTA input",
            "file_loaded": "File loaded:",
            "choose_demo": "Load a demo sequence",
            "reading_frame_preview": "Reading-frame preview",
            "reading_frame_caption": "Six-frame summary to guide the reading-frame choice before analysis.",
            "analyze_button": "Analyze Sequence",
            "batch_complete": "Batch analysis complete",
            "batch_summary": "Batch summary",
            "select_sequence_inspect": "Select a sequence to inspect in this batch:",
            "independent_tools_title": "Independent Analysis Tools",
            "independent_tools_subtitle": "Alignments · Distance Matrix · Phylogeny · Protein Analysis",
            "tab_alignments": "Alignments",
            "tab_distance_matrix": "Distance Matrix",
            "tab_phylogeny": "Phylogeny",
            "tab_protein_analysis": "Protein Analysis",
            "tab_trait_search": "Trait Search",
            "hero_subtitle": "Bioinformatics · AI Interpretation · Agricultural Insights",
            "metadata_load_help": "The app loads lightweight gene metadata first for search and filtering. Full sequence data is loaded only when an analysis is run.",
            "metadata_filter_help": "Filter the loaded gene database by gene_id, symbol, or trait.",
            "metadata_preview_note": "Showing a sample of 20 gene metadata records. Use search to filter specific genes.",
            "metadata_load_error": "Could not load gene metadata preview. Falling back to full database load.",
            "no_genes_available": "No genes available: configure the PostgreSQL connection or add the fallback database.",
            "metadata_count_summary": "Showing {count} matching gene metadata records",
            "welcome_title": "Ready to Analyze",
            "welcome_message": "Paste a plant DNA sequence above or load a demo,<br>then click <b>Analyze Sequence</b>.",
            "what_does_app_analyze": "What does this app analyze?",
            "supported_input_formats": "Supported input formats",
            "independent_tools_notice": "Looking for alignments, distance matrices, phylogeny, or standalone protein analysis? They now live on their own page — see **Independent Tools** in the sidebar navigation.",
            "msa_input_hint": "Paste multiple FASTA sequences or one per line:",
            "run_msa": "Run MSA",
            "sequence_1": "Sequence 1",
            "sequence_2": "Sequence 2",
            "align_pairwise": "Align pairwise",
            "pairwise_missing": "Provide two sequences for pairwise alignment.",
            "distance_matrix_title": "Compute Pairwise Distance Matrix",
            "distance_input_hint": "Paste FASTA or one sequence per line:",
            "method": "Method",
            "compute_distance_matrix": "Compute Distance Matrix",
            "phylogeny_title": "Build Phylogenetic Tree",
            "phylogeny_input_hint": "Paste sequences for phylogeny (FASTA or lines):",
            "tree_algorithm": "Tree algorithm",
            "build_tree": "Build Tree",
            "protein_biochemical_analysis": "Protein biochemical analysis",
            "paste_protein_sequence": "Paste protein sequence:",
            "analyze_protein": "Analyze protein",
            "translation_label": "Translation",
            "language_selector_label": "Language / Langue / Dil",
            "menu": "Menu",
            "title_plant_gene_analyzer": "Plant Gene Analyzer",
            "sequence_analysis_tab": "Sequence Analysis",
        },
        "results": {
            "gene_records_available": "gene metadata records available",
            "genes_loaded": "genes loaded",
            "sequence_length": "sequence length",
            "gc_content": "GC content",
            "best_match": "best match",
            "similarity_score": "similarity score",
            "mutations": "mutations",
            "translation": "translation",
            "protein_analysis": "protein analysis",
        },
        "errors": {
            "no_sequence": "Please enter or paste a sequence before analyzing.",
            "db_connection": "Could not connect to the gene database.",
            "invalid_file": "Invalid file or unsupported encoding.",
            "database_missing": "No genes available: configure the PostgreSQL connection or add the fallback database.",
        },
        "glossary": {
            "gene": "gene",
            "sequence": "sequence",
            "similarity": "similarity",
            "mutation": "mutation",
            "database": "database",
            "metadata": "metadata",
            "translation": "translation",
            "reading_frame": "reading frame",
            "gc_content": "GC content",
            "alignment": "alignment",
            "trait": "trait",
            "species": "species",
            "annotation": "annotation",
            "expression": "expression",
            "protein": "protein",
            "organism": "organism",
            "nucleotide": "nucleotide",
        },
    },
    "tr": {
        "ui": {
            "app_title": "Bitki Gen Analizörü",
            "about": "Hakkında",
            "settings": "Ayarlar",
            "database": "Veritabanı",
            "language": "Dil",
            "sequence_analysis": "Dizi Analizi",
            "independent_tools": "Bağımsız Araçlar",
            "sequence_input": "Dizi Girişi",
            "quick_demo": "Hızlı Demo",
            "load_demo_sequence": "Demo Diziyi Yükle",
            "top_matches": "Gösterilecek en iyi veritabanı eşleşmeleri",
            "top_matches_help": "Gösterilecek en iyi eşleşen gen sayısı.",
            "deep_search": "Derin benzerlik aramasını etkinleştir",
            "deep_search_help": "Hizalama uzunluğu ön filtrasyonunu kapatır ve daha fazla aday değerlendirir. Daha yavaştır, ancak kısa veya farklı diziler için hassasiyeti artırır.",
            "window_size": "Kaydırma penceresi (GC profili)",
            "window_size_help": "GC profili grafiği için pencere boyutu (bp).",
            "reading_frame": "Çeviri için okuma çerçevesi",
            "input_type": "Giriş türü",
            "input_type_help": "Dizi türünü seçin veya uygulamanın otomatik algılamasına izin verin.",
            "input_accepts": "Bu uygulama **{type}** girişi kabul edecek ve analizi buna göre uyarlayacaktır.",
            "auto_detect": "Otomatik algıla",
            "dna": "DNA",
            "protein": "Protein",
            "file_support_help": "'>' başlıklı FASTA dosyaları desteklenir.",
            "sequence_placeholder": "Ham DNA, protein dizisini veya FASTA formatını buraya yapıştırın.\nÖrnek: ATGCGTAGCTAGCGATCGATCGAATTCG...",
            "file_loaded": "Yüklenen dosya:",
            "selected_sequence": "Seçili:",
            "batch_summary_text": "{count} dizi toplu olarak analiz edilecektir.",
            "search_gene": "Gen kimliği, sembolü veya özelliği ara",
            "preview_metadata": "Gen meta verilerini ön izleme",
            "showing_sample": "20 gen meta veri kaydı örneği gösteriliyor. Belirli genleri filtrelemek için arama kullanın.",
            "full_db_load": "Tüm gen veritabanı ve diziler, analiz başlatıldığında yüklenecektir.",
            "build_about": "Bitki DNA dizilerini şu amaçlarla analiz eder:",
            "standalone_tools_message": "Buraya yapıştırdığınız herhangi bir dizi üzerinde çalışan bağımsız araçlar — ana analizi çalıştırmaya gerek yok.",
            "about_items": [
                "GC içeriği ve nükleotid istatistikleri",
                "Gen veritabanı benzerliği",
                "Mutasyon tespiti",
                "Yapay zekâ biyolojik yorumlama",
                "Tarımsal öneriler",
            ],
            "upload_file": "Dizi dosyası yükle (.fasta / .fa / .txt)",
            "paste_sequence": "Veya dizinizi buraya yapıştırın:",
            "select_sequence": "Analiz edilecek bir dizi seçin",
            "analyze_all": "Bu FASTA girişindeki tüm dizileri analiz et",
            "file_loaded": "Yüklenen dosya:",
            "choose_demo": "Bir demo dizisi yükle",
            "reading_frame_preview": "Okuma çerçevesi önizlemesi",
            "reading_frame_caption": "Analiz öncesinde okuma çerçevesini seçmeye yardımcı olacak altı çerçeve özeti.",
            "analyze_button": "Diziyi Analiz Et",
            "batch_complete": "Toplu analiz tamamlandı",
            "batch_summary": "Toplu özet",
            "select_sequence_inspect": "Bu topluluktaki bir diziyi incelemek için seçin:",
            "independent_tools_title": "Bağımsız Analiz Araçları",
            "independent_tools_subtitle": "Hizalamalar · Uzaklık Matrisi · Filogenetik · Protein Analizi",
            "tab_alignments": "Hizalamalar",
            "tab_distance_matrix": "Uzaklık Matrisi",
            "tab_phylogeny": "Filogenetik",
            "tab_protein_analysis": "Protein Analizi",
            "tab_trait_search": "Özellik Araması",
            "hero_subtitle": "Biyoinformatik · Yapay Zeka Yorumu · Tarımsal İçgörüler",
            "metadata_load_help": "Uygulama, arama ve filtreleme için önce hafif gen meta verilerini yükler. Tam dizi verileri yalnızca analiz çalıştırıldığında yüklenir.",
            "metadata_filter_help": "Yüklenen gen veritabanını gene_id, sembol veya özellik ile filtreleyin.",
            "metadata_preview_note": "20 adet gen meta veri kaydı örneği gösteriliyor. Belirli genleri filtrelemek için aramayı kullanın.",
            "metadata_load_error": "Gen meta veri önizlemesi yüklenemedi. Tam veritabanı yüklemesine dönülüyor.",
            "no_genes_available": "Kullanılabilir gen yok: PostgreSQL bağlantısını yapılandırın veya yedek veritabanını ekleyin.",
            "metadata_count_summary": "Aramaya uyan {count} gen meta veri kaydı gösteriliyor",
            "welcome_title": "Analize Hazır",
            "welcome_message": "Yukarıdan bir bitki DNA dizisi yapıştırın veya demo yükleyin,<br>ardından <b>Dizi Analizi</b> butonuna tıklayın.",
            "what_does_app_analyze": "Bu uygulama ne analiz eder?",
            "supported_input_formats": "Desteklenen giriş formatları",
            "independent_tools_notice": "Hizalama, uzaklık matrisi, filogenetik veya bağımsız protein analizi mi arıyorsunuz? Bunlar artık kendi sayfasında — kenar çubuğundaki **Bağımsız Araçlar** bölümünden erişilebilir.",
            "msa_input_hint": "Birden çok FASTA dizisini veya satır başına bir dizi yapıştırın:",
            "run_msa": "Çoklu Hizalama Çalıştır",
            "sequence_1": "Dizi 1",
            "sequence_2": "Dizi 2",
            "align_pairwise": "İkili hizala",
            "pairwise_missing": "İkili hizalama için iki dizi sağlayın.",
            "distance_matrix_title": "İkili Uzaklık Matrisi Hesapla",
            "distance_input_hint": "FASTA veya satır başına bir dizi yapıştırın:",
            "method": "Yöntem",
            "compute_distance_matrix": "Uzaklık Matrisi Hesapla",
            "phylogeny_title": "Filogenetik Ağaç Oluştur",
            "phylogeny_input_hint": "Filogenetik için dizileri yapıştırın (FASTA veya satırlar):",
            "tree_algorithm": "Ağaç algoritması",
            "build_tree": "Ağacı Oluştur",
            "protein_biochemical_analysis": "Protein biyokimyasal analizi",
            "paste_protein_sequence": "Protein dizisini yapıştırın:",
            "analyze_protein": "Proteini analiz et",
            "translation_label": "Çeviri",
            "language_selector_label": "Language / Langue / Dil",
            "menu": "Menü",
            "title_plant_gene_analyzer": "Bitki Gen Analizörü",
            "sequence_analysis_tab": "Dizi Analizi",
        },
        "results": {
            "gene_records_available": "gen meta veri kaydı mevcut",
            "genes_loaded": "gen yüklendi",
            "sequence_length": "dizi uzunluğu",
            "gc_content": "GC içeriği",
            "best_match": "en iyi eşleşme",
            "similarity_score": "benzerlik skoru",
            "mutations": "mutasyonlar",
            "translation": "çeviri",
            "protein_analysis": "protein analizi",
        },
        "errors": {
            "no_sequence": "Analiz öncesinde lütfen bir dizi girin veya yapıştırın.",
            "db_connection": "Gen veritabanına bağlanılamadı.",
            "invalid_file": "Geçersiz dosya veya desteklenmeyen kodlama.",
            "database_missing": "Kullanılabilir gen yok: PostgreSQL bağlantısını yapılandırın veya yedek veritabanını ekleyin.",
        },
        "glossary": {
            "gene": "gen",
            "sequence": "dizi",
            "similarity": "benzerlik",
            "mutation": "mutasyon",
            "database": "veritabanı",
            "metadata": "meta veri",
            "translation": "çeviri",
            "reading_frame": "okuma çerçevesi",
            "gc_content": "GC içeriği",
            "alignment": "hizalama",
            "trait": "özellik",
            "species": "tür",
            "annotation": "açıklama",
            "expression": "ifade",
            "protein": "protein",
            "organism": "organizma",
            "nucleotide": "nükleotit",
        },
    },
}


def _resolve_lang(lang=None):
    if lang in SUPPORTED_LANGS:
        return lang

    if "lang" in st.session_state and st.session_state["lang"] in SUPPORTED_LANGS:
        return st.session_state["lang"]

    return DEFAULT_LANG


def _fallback_chain(lang):
    chain = []
    if lang in SUPPORTED_LANGS:
        chain.append(lang)
    for candidate in ("en", "fr", "tr"):
        if candidate not in chain:
            chain.append(candidate)
    return chain


def _lookup(data, path):
    current = data
    for part in path.split("."):
        if not isinstance(current, dict) or part not in current:
            return None
        current = current[part]
    return current


def translate(key, lang=None, default=None, **kwargs):
    selected_lang = _resolve_lang(lang)
    value = None

    for candidate in _fallback_chain(selected_lang):
        locale = TRANSLATIONS.get(candidate, TRANSLATIONS[DEFAULT_LANG])
        value = _lookup(locale, key)
        if value is not None:
            break

    if value is None:
        value = default or key

    if isinstance(value, str):
        try:
            return value.format(**kwargs) if kwargs else value
        except Exception:
            return value

    return str(value)


def get_glossary(term, lang=None):
    selected_lang = _resolve_lang(lang)
    for candidate in _fallback_chain(selected_lang):
        locale = TRANSLATIONS.get(candidate, TRANSLATIONS[DEFAULT_LANG])
        value = locale.get("glossary", {}).get(term)
        if value is not None:
            return value
    return term


def language_selector(*, key="language_selector", label_key="ui.language_selector_label"):
    options = [(lang, lang.upper()) for lang in SUPPORTED_LANGS]
    selected = st.selectbox(
        translate(label_key, default="Language / Langue / Dil"),
        options=[lang for lang in SUPPORTED_LANGS],
        index=SUPPORTED_LANGS.index(_resolve_lang()),
        key=key,
        format_func=lambda lang: {"fr": "FR", "en": "EN", "tr": "TR"}[lang],
    )
    st.session_state["lang"] = selected
    return selected


def translate_input_type(value, lang=None):
    mapping = {
        "Auto detect": {"fr": "Détection automatique", "en": "Auto detect", "tr": "Otomatik algıla"},
        "DNA": {"fr": "ADN", "en": "DNA", "tr": "DNA"},
        "Protein": {"fr": "Protéine", "en": "Protein", "tr": "Protein"},
    }
    selected_lang = _resolve_lang(lang)
    return mapping.get(value, {}).get(selected_lang, value)


def current_lang():
    return _resolve_lang()
