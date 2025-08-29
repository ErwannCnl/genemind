from pydantic import BaseModel, Field
from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")

from typing import List, Optional
from langchain_core.output_parsers import PydanticOutputParser
from langchain_core.prompts import ChatPromptTemplate

from langchain.chat_models import init_chat_model

# Initialise LLM (Gemini 2.5)
llm = init_chat_model(
    model="gemini-2.5-flash",
    model_provider="google_genai",
    temperature=0.0
)

# Canonical user input
paragraph = (
    "I identified these five genes to be significantly more mutated than expected by chance in my cohort of human brain cancer patients: CCLX, TLR4, TLR2, IL1B. Do a GSEA"
)

# Define Pydantic class for input genes and context
class StudyExtraction(BaseModel):
    genes: List[str] = Field(
        description="List of gene symbols mentioned in the text, normalized to official HGNC/NCBI-style symbols if possible."
    )
    organism: Optional[str] = Field(
        description="Scientific name (binomial) of the organism (e.g., 'Homo sapiens', 'Mus musculus')."
    )
    field_of_study: Optional[str] = Field(
        description="High-level biomedical domain, e.g., 'oncology', 'cancer genomics', 'neuroscience', 'immunology', 'microbiology'."
    )
    organ: Optional[str] = Field(
        description="Primary organ or tissue referenced (e.g., 'brain', 'liver', 'blood')."
    )
    analysis_type: Optional[str] = Field(
        description="Concise description of the analysis performed, e.g., 'differential expression', 'mutation enrichment', 'GWAS', 'copy-number analysis', 'metagenomic profiling'."
    )
    GSEA: bool = Field(default=False, 
                       description="Whether the user mentions that a GSEA is needed on the gene set. If no mention, keep it False.")

#alternative implementation to parse as pydantic more robustly
parser = PydanticOutputParser(pydantic_object=StudyExtraction)
format_instructions = parser.get_format_instructions()

prompt = ChatPromptTemplate.from_messages([
    ("system", "Extract per schema:\n{format_instructions}"),
    ("human", "{paragraph}"),
]).partial(format_instructions=format_instructions)
parsing_llm = prompt | llm | parser

# pass raw user input "paragraph"
parsed_input = parsing_llm.invoke({"paragraph": paragraph})
#convert to JSON string
json_output = parsed_input.model_dump_json(indent=2)
# print(json_output)
print("parsed input")
#next step - inject the json to the LLM to determine attributes to fetch from BioMart
from src.querries_script import group_by_gene_dynamic, fill_with_ncbi, call_querry_biomart
from src.print_gene import format_genes
import pandas as pd    

attributes = pd.read_csv("data/attributes.csv")["name"].to_list()

output = call_querry_biomart(attributes=attributes[:15],
                            filters={"external_gene_name": parsed_input.genes})

output = group_by_gene_dynamic(output)

output = fill_with_ncbi(output)

# print(format_genes(output))
print("got biomart result")
from src.tools import enrichr_query
if parsed_input.GSEA: 
    tool_results = enrichr_query(parsed_input.genes)

        # Optionally, format or reduce the output for readability (e.g., top rows)
    if not tool_results.empty:
        # Filter rows where Adjusted P-value < 0.05
        filtered_df = tool_results[tool_results["Adjusted P-value"] < 0.05]

        # Drop the 'Gene_set' column
        filtered_df = filtered_df.drop(columns=["Gene_set","Old P-value","Old Adjusted P-value"])
        gsea_string = "the results of gene set enrichment are:\n "+ filtered_df.head(20).to_string(index = False)
        # Show the top rows
        # print(gsea_string)

#retrieve gene popularity
from src.gene_lookup import _format_popularity_block
file = "data/all_gene_counts.tsv"
popularity_block = _format_popularity_block(file, parsed_input.genes)
#print(popularity_block)
print("done GSEA")
from langchain_core.output_parsers import StrOutputParser

text = str(output)

task ='''
<task>
You are a helpful and biological expert specializing in integrating and interpreting gene-related data from the given information report, combined with the knowledge you have obtained during your training.
When given a report on a set of genes, briefly summarize the key insights for each gene. The final and most important job is to concisely contextualise the user's findings with biological or biomedical background knowledge and find commonalities between the genes given. Stay scientifically accurate. Tailor the response to the context given by the client.
Think hard!
</task>
'''

user_prompt = "I want to know more about these genes {input_g}, with respective popularities (frequency of citation) {popularity_block}. Report these popularities. They were obtained after {analyses} in {organism}. All information I know about these genes is the following: {text} \n Do you find any commonalities or interesting findings about these genes? I'm mainly interested in {context}"
if parsed_input.GSEA:
    user_prompt += gsea_string
    
prompt = ChatPromptTemplate.from_messages([
    ("system", task),
    ("user", user_prompt)
])
chain = prompt | llm | StrOutputParser()
response = chain.invoke({"text": text, "input_g": parsed_input.genes, "context": parsed_input.field_of_study, "analyses": parsed_input.analysis_type, "organism": parsed_input.organism, "popularity_block": popularity_block})

print("got summary")
from src.print_gene import texts_and_markdown_to_pdf

texts_and_markdown_to_pdf([format_genes(output)], response, "gene_report_from_text.pdf")
print("done")