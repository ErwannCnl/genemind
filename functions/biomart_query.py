from biomart import BiomartServer

def biomart_query(gene_name, dataset="hsapiens_gene_ensembl"):
    """
    Retrieve gene coordinates from Ensembl BioMart for a given gene.
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., 'BRCA1')
    dataset : str, optional
        BioMart dataset (default: human 'hsapiens_gene_ensembl').
        Examples: 'mmusculus_gene_ensembl' for mouse
    
    Returns
    -------
    dict or None
        Dictionary with chromosome, start, end, strand, and Ensembl ID.
        None if gene not found.
    """
    # Connect to Ensembl BioMart
    server = BiomartServer("http://www.ensembl.org/biomart")
    mart = server.datasets[dataset]
    
    # Query BioMart
    response = mart.search({
        'filters': {'hgnc_symbol': gene_name},
        'attributes': [
            'ensembl_gene_id',
            'chromosome_name',
            'start_position',
            'end_position',
            'strand'
        ]
    })
    
    lines = response.text.strip().split("\n")
    
    if not lines or lines[0] == "":
        return None
    
    ensembl_id, chrom, start, end, strand = lines[0].split("\t")
    
    return {
        "ensembl_gene_id": ensembl_id,
        "chromosome": chrom,
        "start": int(start),
        "end": int(end),
        "strand": "+" if strand == "1" else "-"
    }

# Example usage:
# coords = get_gene_coordinates("NOTCH1")
# print(coords)
