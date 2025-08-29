from gseapy import enrichr

def enrichr_query(gene_list: list):
    """Run enrichment analysis on a list of genes.

    This tool allows to run enrichment analysis on a list of genes using the `gseapy` library.
    Using this tool, a model can get information about the biological processes enriched in a set of genes.

    Args:
        gene_list: list of genes to run enrichment analysis on

    Returns:
        DataFrame: DataFrame containing the enrichment results
    """
    # Run enrichment
    enr = enrichr(
        gene_list=gene_list,
        gene_sets='GO_Biological_Process_2021',
        organism='Human',
        outdir=None,  # no files will be written
        cutoff=0.05
    )

    # Save results as DataFrame
    df_results = enr.results

    return df_results