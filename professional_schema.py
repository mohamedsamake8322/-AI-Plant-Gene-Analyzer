"""
professional_schema.py
----------------------
Professional-grade JSON schema for bioinformatics data.
Supports versioning, quality metrics, cloud storage references, and rich metadata.
"""

from __future__ import annotations

from typing import Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime
import hashlib
import re

# ─── Schema Version ────────────────────────────────────────────────────────────
SCHEMA_VERSION = "2.0"
SCHEMA_UPDATED_AT = datetime.utcnow().isoformat() + "Z"


@dataclass
class TaxonomyInfo:
    """Structured taxonomy information."""
    scientific_name: str
    common_name: Optional[str] = None
    taxonomy_id: Optional[int] = None
    genus: Optional[str] = None
    species: Optional[str] = None

    @classmethod
    def from_organism_string(cls, organism_str: str) -> TaxonomyInfo:
        """Parse taxonomy from string like 'Oryza sativa'."""
        parts = organism_str.split()
        genus = parts[0] if len(parts) > 0 else None
        species = parts[1] if len(parts) > 1 else None
        return cls(
            scientific_name=organism_str,
            genus=genus,
            species=species,
        )


@dataclass
class SequenceMetadata:
    """Per-sequence metadata."""
    type: str  # dna, rna, protein, cdna, cds
    length: int
    sequence_hash: Optional[str] = None
    url: Optional[str] = None  # S3 path or cloud reference
    format: str = "fasta"
    gc_content: Optional[float] = None
    
    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass
class ExternalReference:
    """External database reference."""
    database: str  # NCBI-Nucleotide, NCBI-Protein, Ensembl, etc.
    url: str
    id: str
    access_date: Optional[str] = None
    
    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass
class GeneMetadata:
    """Gene-level metadata."""
    created_at: str
    source: str  # NCBI, Ensembl, EBI, etc.
    quality_score: Optional[float] = None  # 0-1
    data_completeness: Optional[float] = None  # 0-1
    accession: Optional[str] = None
    collection_date: Optional[str] = None
    
    def to_dict(self) -> dict:
        return {k: v for k, v in asdict(self).items() if v is not None}


@dataclass
class GeneCore:
    """Core gene identification fields."""
    gene_id: str
    symbol: str
    name: Optional[str] = None
    organism: Optional[TaxonomyInfo] = None
    
    def to_dict(self) -> dict:
        d = asdict(self)
        if self.organism:
            d["organism"] = asdict(self.organism)
        return d


@dataclass
class GeneRecord:
    """Professional gene record with versioning."""
    _id: str  # Unique identifier: source_accession
    metadata: GeneMetadata
    core: GeneCore
    sequences: list[SequenceMetadata] = field(default_factory=list)
    annotations: dict[str, Any] = field(default_factory=dict)
    relationships: dict[str, Any] = field(default_factory=dict)
    data: dict[str, Any] = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        """Convert to dict, excluding None/empty values."""
        return {
            "_id": self._id,
            "metadata": self.metadata.to_dict(),
            "core": self.core.to_dict(),
            "sequences": [s.to_dict() for s in self.sequences],
            "annotations": {k: v for k, v in self.annotations.items() if v},
            "relationships": {k: v for k, v in self.relationships.items() if v},
            "data": {k: v for k, v in self.data.items() if v},
        }


def transform_legacy_gene(legacy_gene: dict) -> GeneRecord:
    """Convert old-style gene record to professional schema."""
    gene_id = legacy_gene.get("gene_id") or legacy_gene.get("symbol") or legacy_gene.get("accession", "unknown")
    source = legacy_gene.get("source", "NCBI")
    
    # Build ID
    record_id = f"gene_{source.lower()}_{gene_id}"
    
    # Metadata
    metadata = GeneMetadata(
        created_at=legacy_gene.get("date_added", datetime.utcnow().isoformat() + "Z"),
        source=source,
        accession=legacy_gene.get("accession") or gene_id,
        quality_score=0.75,  # Default middle score
        data_completeness=_calculate_completeness(legacy_gene),
    )
    
    # Organism
    org_str = legacy_gene.get("organism", "Unknown")
    organism = TaxonomyInfo.from_organism_string(org_str)
    
    # Core
    core = GeneCore(
        gene_id=gene_id,
        symbol=legacy_gene.get("symbol", gene_id),
        name=legacy_gene.get("name") or legacy_gene.get("description"),
        organism=organism,
    )
    
    # Sequences
    sequence_text = legacy_gene.get("sequence")
    sequences = []
    if sequence_text:
        seq_type = legacy_gene.get("sequence_type", "dna")
        seq_len = legacy_gene.get("length", len(sequence_text))
        seq_hash = _compute_sequence_hash(sequence_text)
        gc = _compute_gc_content(sequence_text) if seq_type == "dna" else None
        
        sequences.append(SequenceMetadata(
            type=seq_type,
            length=seq_len,
            sequence_hash=seq_hash,
            format="fasta",
            gc_content=gc,
        ))
    
    # External references
    external_refs = []
    if legacy_gene.get("external_links"):
        for link_key, link_val in legacy_gene.get("external_links", {}).items():
            if link_val:
                if isinstance(link_val, dict):
                    db = link_val.get("database", link_key)
                    url = link_val.get("url", "")
                    id_val = link_val.get("id", gene_id)
                elif isinstance(link_val, str):
                    db = link_key
                    url = link_val
                    id_val = gene_id
                else:
                    continue
                
                if url:
                    external_refs.append(ExternalReference(
                        database=db,
                        url=url,
                        id=id_val,
                    ))
    
    relationships = {
        "external_references": [asdict(ref) for ref in external_refs],
        "homologs": legacy_gene.get("homologs", []),
        "orthologs": legacy_gene.get("orthologs", []),
    }
    
    # Annotations
    annotations = {
        "description": legacy_gene.get("description"),
        "function": legacy_gene.get("function"),
        "keywords": legacy_gene.get("keywords", []),
    }
    
    # Data (expression, pathways, etc)
    data = {
        "expression_profiles": legacy_gene.get("expression_profiles", []),
        "pathways": legacy_gene.get("pathways", []),
        "publications": legacy_gene.get("publications", []),
        "functional_annotations": legacy_gene.get("annotations", {}),
    }
    
    return GeneRecord(
        _id=record_id,
        metadata=metadata,
        core=core,
        sequences=sequences,
        annotations=annotations,
        relationships=relationships,
        data=data,
    )


def get_legacy_compatible_view(pro_gene: GeneRecord) -> dict:
    """
    Convert professional schema back to legacy format for app compatibility.
    This allows the app to work without modification.
    """
    legacy = {
        "gene_id": pro_gene.core.gene_id,
        "symbol": pro_gene.core.symbol,
        "organism": pro_gene.core.organism.scientific_name if pro_gene.core.organism else "Unknown",
        "description": pro_gene.annotations.get("description"),
        "source": pro_gene.metadata.source,
        "date_added": pro_gene.metadata.created_at,
        "accession": pro_gene.metadata.accession,
    }
    
    # Add first sequence
    if pro_gene.sequences:
        seq = pro_gene.sequences[0]
        legacy["sequence"] = ""  # Don't inline for compatibility
        legacy["sequence_type"] = seq.type
        legacy["length"] = seq.length
    
    # External links
    if pro_gene.relationships.get("external_references"):
        legacy["external_links"] = {
            ref["database"]: {
                "url": ref["url"],
                "id": ref["id"],
            }
            for ref in pro_gene.relationships["external_references"]
        }
    
    # Data
    legacy["expression_profiles"] = pro_gene.data.get("expression_profiles", [])
    legacy["pathways"] = pro_gene.data.get("pathways", [])
    legacy["publications"] = pro_gene.data.get("publications", [])
    legacy["annotations"] = pro_gene.data.get("functional_annotations", {})
    
    return legacy


# ─── Helper functions ─────────────────────────────────────────────────────────

def _calculate_completeness(gene: dict) -> float:
    """Calculate data completeness score (0-1)."""
    score = 0.0
    fields = ["gene_id", "symbol", "sequence", "organism", "description"]
    for field in fields:
        if gene.get(field):
            score += 0.2
    return min(score, 1.0)


def _compute_sequence_hash(sequence: str) -> str:
    """Compute SHA256 hash of sequence."""
    return "sha256:" + hashlib.sha256(sequence.upper().encode()).hexdigest()[:16]


def _compute_gc_content(dna_seq: str) -> float:
    """Compute GC content of DNA sequence."""
    seq = dna_seq.upper()
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) if len(seq) > 0 else 0.0


def create_schema_metadata(
    total_records: int,
    unique_organisms: int,
    sequence_types: list[str],
    validated_count: int = 0,
    with_annotations: int = 0,
) -> dict:
    """Create schema-level metadata."""
    return {
        "schema": {
            "version": SCHEMA_VERSION,
            "updated_at": SCHEMA_UPDATED_AT,
            "description": "Multi-type plant sequence database (professional)",
        },
        "metadata": {
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "collection_method": "Automated Pipeline",
            "total_records": total_records,
            "unique_organisms": unique_organisms,
            "sequence_types": sequence_types,
            "data_quality": {
                "validated": validated_count or total_records,
                "with_annotations": with_annotations,
                "with_expression_data": 0,
            },
        },
    }
