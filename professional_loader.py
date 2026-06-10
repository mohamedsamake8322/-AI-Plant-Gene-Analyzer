"""
professional_loader.py
----------------------
Loader for professional schema with backward compatibility for legacy app.
Allows seamless switching between legacy and professional JSON formats.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Union

import professional_schema as pschema


class ProfessionalGeneDatabase:
    """
    Database interface that works with professional schema,
    but provides legacy-compatible views for existing app code.
    """

    def __init__(self, json_path: Union[str, Path]):
        """Load professional schema JSON."""
        self.path = Path(json_path)
        self.pro_genes: dict[str, pschema.GeneRecord] = {}
        self.legacy_view: dict[str, dict] = {}
        self._load()

    def _load(self) -> None:
        """Load and parse professional schema JSON."""
        if not self.path.exists():
            raise FileNotFoundError(f"Database not found: {self.path}")

        with open(self.path, "r", encoding="utf-8") as f:
            data = json.load(f)

        # Extract metadata
        self.schema_version = data.get("schema", {}).get("version", "1.0")
        self.metadata = data.get("metadata", {})

        # Parse genes
        for gene_dict in data.get("genes", []):
            gene_record = self._dict_to_record(gene_dict)
            if gene_record:
                self.pro_genes[gene_record._id] = gene_record
                # Create legacy view
                self.legacy_view[gene_record.core.gene_id] = pschema.get_legacy_compatible_view(
                    gene_record
                )

    def _dict_to_record(self, gene_dict: dict) -> Optional[pschema.GeneRecord]:
        """Convert dict to GeneRecord dataclass."""
        try:
            core_data = gene_dict.get("core", {})
            org_data = core_data.get("organism", {})
            organism = pschema.TaxonomyInfo(**org_data) if org_data else None

            core = pschema.GeneCore(
                gene_id=core_data.get("gene_id", "unknown"),
                symbol=core_data.get("symbol", "unknown"),
                name=core_data.get("name"),
                organism=organism,
            )

            meta_data = gene_dict.get("metadata", {})
            metadata = pschema.GeneMetadata(
                created_at=meta_data.get("created_at", ""),
                source=meta_data.get("source", "Unknown"),
                quality_score=meta_data.get("quality_score"),
                data_completeness=meta_data.get("data_completeness"),
                accession=meta_data.get("accession"),
            )

            sequences = [
                pschema.SequenceMetadata(
                    type=s.get("type", "dna"),
                    length=s.get("length", 0),
                    sequence_hash=s.get("sequence_hash"),
                    url=s.get("url"),
                    format=s.get("format", "fasta"),
                    gc_content=s.get("gc_content"),
                )
                for s in gene_dict.get("sequences", [])
            ]

            return pschema.GeneRecord(
                _id=gene_dict.get("_id", "unknown"),
                metadata=metadata,
                core=core,
                sequences=sequences,
                annotations=gene_dict.get("annotations", {}),
                relationships=gene_dict.get("relationships", {}),
                data=gene_dict.get("data", {}),
            )
        except Exception as e:
            print(f"Warning: Could not parse gene record: {e}")
            return None

    def get_legacy_dict(self) -> dict:
        """
        Return database in legacy format.
        This allows existing code to work without modification.
        """
        return self.legacy_view

    def get_professional_dict(self) -> dict:
        """Return database in professional format."""
        return {
            **pschema.create_schema_metadata(
                len(self.pro_genes),
                len(set(g.core.organism.scientific_name for g in self.pro_genes.values()
                       if g.core.organism)),
                list(set(s.type for g in self.pro_genes.values() for s in g.sequences)),
            ),
            "genes": [g.to_dict() for g in self.pro_genes.values()],
        }

    def get_gene_legacy(self, gene_id: str) -> Optional[dict]:
        """Get gene in legacy format by ID."""
        return self.legacy_view.get(gene_id)

    def get_gene_professional(self, gene_id: str) -> Optional[pschema.GeneRecord]:
        """Get gene in professional format by ID."""
        for pro_gene in self.pro_genes.values():
            if pro_gene.core.gene_id == gene_id:
                return pro_gene
        return None

    def list_organisms(self) -> list[str]:
        """List all organisms in database."""
        orgs = set()
        for gene in self.pro_genes.values():
            if gene.core.organism:
                orgs.add(gene.core.organism.scientific_name)
        return sorted(list(orgs))

    def list_sequence_types(self) -> list[str]:
        """List all sequence types in database."""
        types = set()
        for gene in self.pro_genes.values():
            for seq in gene.sequences:
                types.add(seq.type)
        return sorted(list(types))

    def filter_by_organism(self, organism: str) -> dict[str, dict]:
        """Get all genes for an organism (legacy format)."""
        result = {}
        for gene_id, legacy_gene in self.legacy_view.items():
            if legacy_gene.get("organism", "").lower() == organism.lower():
                result[gene_id] = legacy_gene
        return result

    def filter_by_sequence_type(self, seq_type: str) -> dict[str, dict]:
        """Get all genes with specific sequence type (legacy format)."""
        result = {}
        for pro_gene in self.pro_genes.values():
            if any(s.type == seq_type for s in pro_gene.sequences):
                legacy = pschema.get_legacy_compatible_view(pro_gene)
                result[pro_gene.core.gene_id] = legacy
        return result

    def __len__(self) -> int:
        return len(self.pro_genes)

    def __repr__(self) -> str:
        return (
            f"ProfessionalGeneDatabase({self.path.name}, "
            f"{len(self.pro_genes)} genes, "
            f"v{self.schema_version})"
        )


def load_database(db_path: Union[str, Path]) -> Union[ProfessionalGeneDatabase, dict]:
    """
    Smart loader: detects schema version and returns appropriate object.
    - Professional schema (v2.0): Returns ProfessionalGeneDatabase
    - Legacy schema (v1.x): Returns dict (legacy behavior)
    """
    db_path = Path(db_path)

    with open(db_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Check schema version
    schema_version = data.get("schema", {}).get("version", "1.0")

    if schema_version.startswith("2"):
        # Professional schema
        db = ProfessionalGeneDatabase(db_path)
        # For backward compatibility, also support dict-like access
        db.__getitem__ = lambda k: db.get_gene_legacy(k)
        return db
    else:
        # Legacy schema - just return dict
        return data


# ─── Integration with existing code ────────────────────────────────────────────

def adapt_for_legacy_app(db: Union[ProfessionalGeneDatabase, dict]) -> dict:
    """
    Convert database to legacy format for use in existing app.
    Works with both professional and legacy schemas.
    """
    if isinstance(db, ProfessionalGeneDatabase):
        return db.get_legacy_dict()
    else:
        # Already legacy format
        if "genes" in db:
            return {g.get("gene_id", g.get("symbol", "unknown")): g for g in db.get("genes", [])}
        else:
            # Already a dict keyed by gene_id
            return db
