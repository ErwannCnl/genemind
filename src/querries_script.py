import requests
from xml.etree.ElementTree import Element, SubElement, tostring
from collections import defaultdict

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
    if "go_linkage_type" not in attributes:
        attributes.append("go_linkage_type")
    if "name_1006" not in attributes:
        attributes.append("name_1006")
    if "namespace_1003" not in attributes:
        attributes.append("namespace_1003")
        
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
    
    tsv = biomart_query(
        attributes=attr,
        filters=filters,
        dataset=dataset
    )
    
    targets = ["go_id", "go_linkage_type", "name_1006"]
    if any(t in attributes for t in targets):
        final_rows = filter_rows(tsv, attributes)
        return final_rows
    else:
        return parse_rows(tsv, headers=attributes)


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