"""
PROFESSIONAL_SCHEMA_APP_PATCH.py
---------------------------------
Example patches for integrating professional schema into app.py.
Copy these changes to your app to unlock new features.

All changes are OPTIONAL - app works without them!
"""

# ════════════════════════════════════════════════════════════════════════════
# PATCH 1: Update imports (optional, for professional schema features)
# ════════════════════════════════════════════════════════════════════════════

# BEFORE (in app.py imports section):
"""
import streamlit as st
import bioinformatics as bio
import similarityengine as sim
"""

# AFTER (add these optional imports):
"""
import streamlit as st
import bioinformatics as bio
import similarityengine as sim

# NEW: Optional imports for professional schema features
try:
    from professional_loader import load_database, ProfessionalGeneDatabase, adapt_for_legacy_app
    PROFESSIONAL_SCHEMA_AVAILABLE = True
except ImportError:
    PROFESSIONAL_SCHEMA_AVAILABLE = False
"""


# ════════════════════════════════════════════════════════════════════════════
# PATCH 2: Update gene database loader (optional, for auto-detection)
# ════════════════════════════════════════════════════════════════════════════

# BEFORE (in app.py):
"""
def load_gene_database_cached(db_path: str = "genes_database.json") -> dict:
    @st.cache_resource
    def _load():
        try:
            if load_gene_database_from_postgres is not None:
                try:
                    db = load_gene_database_from_postgres()
                    if db:
                        return db
                except Exception:
                    pass
        except Exception:
            pass
        db = sim.load_gene_database(db_path)
        return db
    return _load()
"""

# AFTER (with professional schema support):
"""
def load_gene_database_cached(db_path: str = "genes_database.json") -> dict:
    @st.cache_resource
    def _load():
        try:
            if load_gene_database_from_postgres is not None:
                try:
                    db = load_gene_database_from_postgres()
                    if db:
                        return db
                except Exception:
                    pass
        except Exception:
            pass
        
        # NEW: Try professional schema first
        if PROFESSIONAL_SCHEMA_AVAILABLE:
            try:
                db_obj = load_database(db_path)
                # Convert to legacy format for existing code compatibility
                db = adapt_for_legacy_app(db_obj)
                st.sidebar.info("📊 Using professional schema v2.0")
                return db
            except Exception:
                pass
        
        # Fallback to legacy loader
        db = sim.load_gene_database(db_path)
        return db
    return _load()
"""


# ════════════════════════════════════════════════════════════════════════════
# PATCH 3: Add professional schema info sidebar (optional feature showcase)
# ════════════════════════════════════════════════════════════════════════════

# Add this function to app.py:
"""
def show_professional_schema_info(db_path: str = "genes_database.json"):
    '''Display professional schema information in sidebar.'''
    if not PROFESSIONAL_SCHEMA_AVAILABLE:
        return
    
    try:
        from professional_loader import ProfessionalGeneDatabase
        db = ProfessionalGeneDatabase(db_path)
        
        with st.sidebar.expander("📊 Data Quality Metrics"):
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Genes", len(db))
            with col2:
                st.metric("Organisms", len(db.list_organisms()))
            with col3:
                types = ", ".join(db.list_sequence_types())
                st.write(f"**Sequence Types**\\n{types}")
            
            # Show sample quality scores
            pro_genes = list(db.pro_genes.values())
            if pro_genes:
                quality_scores = [g.metadata.quality_score for g in pro_genes if g.metadata.quality_score]
                if quality_scores:
                    avg_quality = sum(quality_scores) / len(quality_scores)
                    st.metric("Avg Quality Score", f"{avg_quality:.2f}")
        
        with st.sidebar.expander("🌾 Organisms"):
            for org in db.list_organisms():
                genes_count = len(db.filter_by_organism(org))
                st.write(f"• {org}: {genes_count} genes")
    
    except Exception as e:
        pass  # Silently fail - professional features optional
"""

# Call this in your main() function:
"""
def main():
    # ... existing code ...
    db = load_gene_database_cached("genes_database.json")
    
    # NEW: Show professional schema info
    if PROFESSIONAL_SCHEMA_AVAILABLE:
        show_professional_schema_info("genes_database.json")
    
    # ... rest of app ...
"""


# ════════════════════════════════════════════════════════════════════════════
# PATCH 4: Advanced search with professional schema
# ════════════════════════════════════════════════════════════════════════════

# Add this function for organism-specific searches:
"""
def search_genes_by_organism(db_path: str, organism: str) -> list:
    '''Search genes by organism using professional schema.'''
    if not PROFESSIONAL_SCHEMA_AVAILABLE:
        return []
    
    try:
        from professional_loader import ProfessionalGeneDatabase
        db = ProfessionalGeneDatabase(db_path)
        results = db.filter_by_organism(organism)
        return list(results.values())
    except Exception:
        return []

def search_genes_by_type(db_path: str, seq_type: str) -> list:
    '''Search genes by sequence type using professional schema.'''
    if not PROFESSIONAL_SCHEMA_AVAILABLE:
        return []
    
    try:
        from professional_loader import ProfessionalGeneDatabase
        db = ProfessionalGeneDatabase(db_path)
        results = db.filter_by_sequence_type(seq_type)
        return list(results.values())
    except Exception:
        return []
"""


# ════════════════════════════════════════════════════════════════════════════
# PATCH 5: Display professional metadata in gene details
# ════════════════════════════════════════════════════════════════════════════

# Update gene detail display to show professional metadata:
"""
def display_gene_details(gene_id: str, db_path: str = "genes_database.json"):
    '''Display gene details with professional schema metadata if available.'''
    
    # ... existing code ...
    
    # NEW: Add professional metadata if available
    if PROFESSIONAL_SCHEMA_AVAILABLE:
        try:
            from professional_loader import ProfessionalGeneDatabase
            db = ProfessionalGeneDatabase(db_path)
            pro_gene = db.get_gene_professional(gene_id)
            
            if pro_gene:
                st.divider()
                st.subheader("📊 Professional Metadata")
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Quality Score", pro_gene.metadata.quality_score or "N/A")
                with col2:
                    completeness = pro_gene.metadata.data_completeness or 0
                    st.metric("Data Completeness", f"{completeness:.0%}")
                with col3:
                    st.metric("Source", pro_gene.metadata.source)
                
                # Show taxonomy
                if pro_gene.core.organism:
                    with st.expander("🌳 Taxonomy Details"):
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**Genus**: {pro_gene.core.organism.genus}")
                        with col2:
                            st.write(f"**Species**: {pro_gene.core.organism.species}")
                        if pro_gene.core.organism.taxonomy_id:
                            st.write(f"**NCBI Taxonomy ID**: {pro_gene.core.organism.taxonomy_id}")
                
                # Show sequence details
                if pro_gene.sequences:
                    with st.expander("🧬 Sequence Details"):
                        for seq in pro_gene.sequences:
                            st.write(f"**Type**: {seq.type}")
                            st.write(f"**Length**: {seq.length:,} bp")
                            if seq.gc_content:
                                st.write(f"**GC Content**: {seq.gc_content:.1%}")
                            if seq.sequence_hash:
                                st.write(f"**Hash**: {seq.sequence_hash}")
        
        except Exception:
            pass  # Silently fail - professional features optional
"""


# ════════════════════════════════════════════════════════════════════════════
# PATCH 6: Export with professional metadata (complete example)
# ════════════════════════════════════════════════════════════════════════════

"""
def export_with_professional_metadata(genes: list, db_path: str) -> dict:
    '''Export genes with professional schema metadata.'''
    if not PROFESSIONAL_SCHEMA_AVAILABLE:
        return {}  # Return legacy format
    
    from professional_loader import ProfessionalGeneDatabase
    db = ProfessionalGeneDatabase(db_path)
    
    enhanced_genes = []
    for gene_id in genes:
        pro_gene = db.get_gene_professional(gene_id)
        if pro_gene:
            # Include professional metadata in export
            export_data = pro_gene.to_dict()
            enhanced_genes.append(export_data)
    
    return {
        "schema": {
            "version": "2.0",
            "description": "Professional schema export with metadata"
        },
        "metadata": {
            "export_date": datetime.utcnow().isoformat() + "Z",
            "gene_count": len(enhanced_genes),
        },
        "genes": enhanced_genes
    }
"""


# ════════════════════════════════════════════════════════════════════════════
# COMPLETE MINIMAL EXAMPLE
# ════════════════════════════════════════════════════════════════════════════

"""
# Minimal changes needed for professional schema support:

# 1. Add to imports
from professional_loader import load_database, adapt_for_legacy_app

# 2. In load_gene_database_cached():
try:
    db_obj = load_database(db_path)
    return adapt_for_legacy_app(db_obj)  # Works with both v1.0 and v2.0
except:
    return sim.load_gene_database(db_path)  # Fallback to legacy

# That's it! Your app now supports professional schema automatically.
"""


if __name__ == "__main__":
    print("📋 PROFESSIONAL SCHEMA APP INTEGRATION PATCHES")
    print("=" * 60)
    print()
    print("This file contains optional patches to integrate")
    print("professional schema features into your app.py.")
    print()
    print("✅ Copy the patches you want to use")
    print("✅ All changes are backward compatible")
    print("✅ Professional features are completely optional")
    print()
    print("Start with PATCH 2 for automatic schema detection.")
    print("Add PATCH 3+ for professional feature showcase.")
