import requests
from xml.etree.ElementTree import Element, SubElement, tostring
from collections import defaultdict
from Bio import Entrez
import xml.etree.ElementTree as ET

# Let the llm decide what attributes it will querry according to the concept we gave him

def biomart_query(attributes, filters=None, dataset=None,
                  host="https://www.ensembl.org",
                  formatter="TSV", header=True, unique_rows=True):
    
    def build_xml():
        q = Element("Query", {
            "virtualSchemaName": "default",
            "formatter": formatter,
            "header": "1" if header else "0",
            "uniqueRows": "1" if unique_rows else "0",
            "count": "0",
            "datasetConfigVersion": "0.6",
        })
        d = SubElement(q, "Dataset", {"name": dataset, "interface": "default"})
        if filters:
            for k, v in filters.items():
                if isinstance(v, (list, tuple, set)):
                    v = ",".join(map(str, v))
                SubElement(d, "Filter", {"name": k, "value": str(v)})
        for a in attributes:
            SubElement(d, "Attribute", {"name": a})
        return tostring(q, encoding="utf-8", method="xml").decode()

    url = host.rstrip("/") + "/biomart/martservice"
    xml_query = build_xml()
    r = requests.post(url, data={"query": xml_query}, timeout=120)
    r.raise_for_status()
    return r.text  # TSV/CSV


def parse_rows(tsv, headers):
    rows = []
    for line in tsv.strip().splitlines():
        cols = line.split("\t")
        rows.append({h: (cols[i] if i < len(cols) else "") for i, h in enumerate(headers)})
    return rows

def filter_rows(tsv, attributes:list):
    
    rows = parse_rows(tsv, headers=attributes)

    keep_domains = {"biological_process"} 
    exclude_evidence = {"IEA", "RCA", "NAS", "ND"}
    rows = [
        r for r in rows
        if (not keep_domains or r["namespace_1003"] in keep_domains)
        and r["go_linkage_type"] not in exclude_evidence
    ]

    evidence_rank = {
        "IDA": 0, "IMP": 0, "IPI": 0, "IGI": 0, "IEP": 0, "EXP": 0,
        "IC": 1, "TAS": 1,
        "ISS": 2, "ISO": 2, "ISA": 2, "ISM": 2, "IGC": 2, "IBA": 2, "IBD": 2, "IKR": 2, "IRD": 2,
        "IEA": 9, "RCA": 9, "NAS": 9, "ND": 9,
    }
    best = {}
    for r in rows:
        key = (r["ensembl_gene_id"], r["go_id"])
        cur = best.get(key)
        score = evidence_rank.get(r["go_linkage_type"], 5)
        if cur is None or score < evidence_rank.get(cur["go_linkage_type"], 5):
            best[key] = r
    rows = list(best.values())
    
    
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["ensembl_gene_id"]].append(r)

    def rank(r):
        return (evidence_rank.get(r["go_linkage_type"], 5), r["name_1006"])

    K = 5
    final_rows = []
    for g, lst in by_gene.items():
        lst.sort(key=rank)
        final_rows.extend(lst[:K])
    
    return final_rows
        
def filter_attributes(proposed_attributes):
    possible_attributes=["ensembl_gene_id",
            "go_id",
            "external_gene_name",
            "name_1006",
            "ensembl_gene_id",
            "chromosome_name",
            "start_position",
            "end_position",
            "strand",
            "name_1006", 
            "namespace_1003",
            "go_linkage_type"
            ]
    
    final_attr = []
    for attr in proposed_attributes:
        if attr in possible_attributes:
            final_attr.append(attr)
    
    return final_attr


def call_querry(attributes, filters:dict, dataset:str="hsapiens_gene_ensembl"):
    
    attr = filter_attributes(proposed_attributes=attributes)
    

    
    targets = ["go_id", "go_linkage_type", "name_1006"]
    if any(t in attr for t in targets):
        if "go_linkage_type" not in attributes:
            attr.append("go_linkage_type")
        if "name_1006" not in attributes:
            attr.append("name_1006")
        if "namespace_1003" not in attributes:
            attr.append("namespace_1003")
            
        tsv = biomart_query(
            attributes=attr,
            filters=filters,
            dataset=dataset
        )
        final_rows = filter_rows(tsv, attr)
        return final_rows
    
    else:
        tsv = biomart_query(
            attributes=attr,
            filters=filters,
            dataset=dataset
        )
        return parse_rows(tsv, headers=attr)[1:]



# Execution example
if __name__ == "__main__":
    print()
    output = call_querry(attributes=["ensembl_gene_id",
                            "go_id",
                            "external_gene_name",
                            "name_1006",
                            "ensembl_gene_id",
                            "chromosome_name",
                            "start_position",
                            "end_position",
                            "strand",
                            "name_1006", 
                            "namespace_1003",
                            "go_linkage_type"
                            ],
                filters={"external_gene_name": ["TLR4", "NOTCH1"]})
    
    print(type(output))
    print(len(output))





def get_gene_info(gene_symbol: str, organism: str = "Homo sapiens") -> str:
    """Fetch NCBI gene summary and representative expression info as text."""
    Entrez.email = "tom.luijts@ugent.be"
    # Step 1: Search for gene ID
    handle = Entrez.esearch(db="gene", term=f"{gene_symbol}[Gene Name] AND {organism}[Organism]")
    record = Entrez.read(handle)
    handle.close()

    gene_ids = record["IdList"]
    if not gene_ids:
        return f"No gene found for {gene_symbol} in {organism}."

    gene_id = gene_ids[0]

    # Step 2: Fetch full gene record (XML)
    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
    gene_data = handle.read()
    handle.close()

    # Step 3: Parse XML
    root = ET.fromstring(gene_data)

    output_lines = [f"Gene: {gene_symbol} ({organism})"]

    # --- Extract Summary ---
    summary = root.find(".//Entrezgene_summary")
    if summary is not None:
        output_lines.append("\n--- Summary ---")
        output_lines.append(summary.text.strip())
    else:
        output_lines.append("\n--- Summary ---")
        output_lines.append("No summary available.")

    # --- Extract Representative Expression ---
    tissues = []
    expr_category = None
    expr_summary = None

    for gc in root.findall(".//Gene-commentary"):
        heading = gc.find("Gene-commentary_heading")
        if heading is not None and heading.text == "Representative Expression":
            for comment in gc.findall(".//Gene-commentary"):
                label = comment.find("Gene-commentary_label")
                if label is not None:
                    if label.text == "Tissue List":
                        tissue_text = comment.find("Gene-commentary_text").text.strip()
                        tissues = [t.strip() for t in tissue_text.split(";") if t.strip()]
                    elif label.text == "Category":
                        expr_category = comment.find("Gene-commentary_text").text.strip()
                    elif label.text == "Text Summary":
                        expr_summary = comment.find("Gene-commentary_text").text.strip()

    # Add to output
    if expr_summary or expr_category or tissues:
        output_lines.append("\n--- Representative Expression ---")
        if expr_summary:
            output_lines.append(f"Summary: {expr_summary}")
        if expr_category:
            output_lines.append(f"Category: {expr_category}")
        if tissues:
            output_lines.append("Tissues: " + ", ".join(tissues))
    else:
        output_lines.append("\n--- Representative Expression ---")
        output_lines.append("No expression info available.")

    return "\n".join(output_lines)


# # Example usage
# if __name__ == "__main__":
#     text_output = get_gene_info("TP53")
#     print(text_output)
