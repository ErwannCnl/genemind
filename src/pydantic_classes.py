#pydantic classes
from typing import List, Optional
from pydantic import BaseModel, Field

class Standard_gene(BaseModel):
    gene_name: str = Field(description="Official HGNC gene symbol, one per gene")
    start: int = Field(description="start location of the gene")
    end: int = Field(default_factory=list, description="end location of the gene")
    replicon: str = Field(description="name of the a region of an organism's genome that is independently replicated. A genomic element e.g., chromosome, plasmid, scaffold, contig")
    normal_biological_function: Optional[str] = Field(default = None, description="short summary of the biological function of the gene, when working properly. if known")
    implicated_diseases: Optional[List[str]] = Field(default_factory=list, description="diseases in which the gene is implicated. Only the disease name, no explanation. All cancer subtypes fall under the global term 'cancer'")
    interactors: Optional[List[str]] = Field(default_factory=list, description="gene names from direct interactors of the gene, if any")
    cellular_compartments: Optional[List[str]] = Field(default_factory=list, description="Cellular compartments in which the expressed gene is mostly present, if known")
    expressed_tissue: Optional[List[str]] = Field(default_factory=list, description="Human tissues in which the gene is expressed.")

class Cancer_gene(BaseModel):
    gene_name: str = Field(description="Official HGNC gene symbol, one per gene")
    hallmarks: Optional[List[str]] = Field(default_factory=list, description="Cancer hallmarks in which the gene is implicated.")
    cancer_types: Optional[List[str]] = Field(default_factory=list, description="Cancer types in which the gene is implicated, if any.")
    cancer_gene_type: Optional[str] = Field(default=None, description="Category of cancer gene, e.g. oncogene, tumor suppressor gene.")
    cancer_function: Optional[str] = Field(default=None,description="Mechanistic role in tumorigenesis when dysregulated")
     