from pydantic import BaseModel
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

def print_pydantic(model: BaseModel, title_field: str = None) -> str:
    """
    Pretty-print any Pydantic model in a structured format.
    
    Args:
        model: Pydantic model instance
        title_field: (optional) field name to use as the title (e.g. "gene_name")
    """
    result = []
    
    # Use title if available
    if title_field and hasattr(model, title_field):
        result.append(f"### {getattr(model, title_field)}")
    else:
        result.append(f"### {model.__class__.__name__}")
    
    result.append("```")
    for field_name, field_value in model.model_dump().items():
        # Handle lists
        if isinstance(field_value, list):
            value_str = "[" + ", ".join(map(str, field_value)) + "]"
        else:
            value_str = str(field_value)
        result.append(f"{field_name}: {value_str}")
    result.append("```")
    
    return "\n".join(result)