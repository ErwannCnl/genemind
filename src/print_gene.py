def format_genes(gene_list):
    """
    Pretty-print gene dictionaries with external_gene_name as the title
    and each key: value on its own line.
    """
    result = []
    for gene in gene_list:
        title = f"### {gene.get('external_gene_name', 'UNKNOWN')}"
        result.append(title)
        result.append("```")
        for key, value in gene.items():
            # Turn lists into readable strings
            if isinstance(value, list):
                value_str = "[" + ", ".join(map(str, value)) + "]"
            else:
                value_str = str(value)
            result.append(f"{key}: {value_str}")
        result.append("```")
        result.append("")  # blank line after each gene
    return "\n".join(result)