from pydantic import BaseModel
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from markdown import markdown   
from bs4 import BeautifulSoup   

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


## print to PDF
def markdown_to_story(md_text: str):
    """
    Convert a markdown string into a list of ReportLab flowables.
    """
    styles = getSampleStyleSheet()
    story = []

    # Convert markdown → HTML → parse with BeautifulSoup
    html = markdown(md_text)
    soup = BeautifulSoup(html, "html.parser")

    for element in soup.descendants:
        if element.name == "h1" or element.name == "h2" or element.name == "h3":
            story.append(Paragraph(f"<b>{element.get_text()}</b>", styles["Heading2"]))
            story.append(Spacer(1, 12))
        elif element.name == "p":
            story.append(Paragraph(element.get_text(), styles["Normal"]))
            story.append(Spacer(1, 6))
        elif element.name == "li":
            story.append(Paragraph("• " + element.get_text(), styles["Normal"]))
            story.append(Spacer(1, 4))

    return story

def texts_and_markdown_to_pdf(text_outputs: list[str], md_text: str, filename: str):
    """
    Combine structured text outputs and markdown content into a single PDF.
    """
    styles = getSampleStyleSheet()
    story = []

    # Add all structured text outputs
    for idx, text in enumerate(text_outputs):
        for line in text.splitlines():
            if line.startswith("### "):  # Title
                story.append(Paragraph(f"<b>{line[4:]}</b>", styles["Title"]))
                story.append(Spacer(1, 12))
            elif line.strip() == "```":  # skip markdown fences
                continue
            elif line.strip():  # Normal key-value lines
                story.append(Paragraph(line, styles["Normal"]))
                story.append(Spacer(1, 6))
            else:  # blank line
                story.append(Spacer(1, 12))

        # Page break between genes
        story.append(PageBreak())

    # Add markdown summary
    story += markdown_to_story(md_text)

    # Build PDF
    doc = SimpleDocTemplate(filename, pagesize=letter)
    doc.build(story)