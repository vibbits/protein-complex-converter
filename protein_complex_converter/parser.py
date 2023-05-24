import csv, time
from dataclasses import dataclass
from io import StringIO
from typing import Iterable
from bs4 import BeautifulSoup
from Bio import Entrez

# Define classes for the data to use
@dataclass
class Row:
    complex: str
    uniprot_ids: list[str]
    pubmed_ids: list[str]
    first_author_name: list[str]
    source: str
    taxid: str

@dataclass
class MitabRow:
    uida: str
    uidb: str
    author: str
    pmids: str
    sourcedb: str
    
    
# Extract the identifiers of molecules in complex
def extract_id(name: str) -> str:
    return name.split("(")[0]

# Extract pubmed identifiers
def extract_pubmed_id(crossref:str) -> str:
    return crossref.split("(")[0] 

# Parse pubmed xml to get first author and publication year
def get_pubmedArticles(pubmed_ids):
    pubmed_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&rettype=abstract"
    response = req.get(pubmed_base_url+'&id='+str(pubmed_ids))
    xml = response.text
    pubmedArticleSet = et.fromstring(xml)
    return xml

# Extract first authors using pubmed ids
def extract_first_author(xml:str) -> str:
    xml_data = xml
    soup = BeautifulSoup(xml_data, 'xml')
    author_list = soup.find_all('Author')
    first_author = author_list[0]
    last_name = first_author.find('LastName').text
    publication_year = soup.find('PubDate').find('Year').text
    last_name_with_year = last_name + " " + "et al. " + "(" + publication_year + ")"
    return last_name_with_year

def extract_taxid(taxid:int) -> str:
    return "taxid:"+taxid+"(Homo sapiens)"

# Get the complex, uniprot and pubmed identifiers and parse first author from pmid
def parse_complex_tab(ct: str) -> None:
    reader = csv.DictReader(ct.splitlines(), delimiter="\t")
    rows = []
    for row in reader:
        first_author_name = []
        pubmed_ids = [extract_pubmed_id(crossref) for crossref in row["Cross references"
        ].split("|") if 'pubmed:' in crossref]
        for pmid in pubmed_ids:
            first_author_name.append(extract_first_author(Entrez.efetch(db='pubmed', id=pmid.split('pubmed:')[1], retmode='xml')))
            time.sleep(0.3)
        rows.append(Row(
            complex=row["#Complex ac"],
            uniprot_ids=[extract_id(name) for name in row[
                "Identifiers (and stoichiometry) of molecules in complex"
            ].split("|")],
            pubmed_ids= '|'.join(pubmed_ids),
            first_author_name = '|'.join(first_author_name),
            source = row["Source"],
            taxid = extract_taxid(taxid) for taxid in row["Taxonomy identifier"]
        ))
    return rows

# Convert the complex ids and uniprot ids into mitab format
def convert_to_mitab(rows: list[Row]) -> list[MitabRow]:
    for row in rows:
        for id in row.uniprot_ids:
            yield MitabRow(
                uida = row.complex, 
                uidb = id,
                author = row.first_author_name,
                pmids = row.pubmed_ids,
                sourcedb = row.source
            )

# Serialize the mitab results
def serialize_to_mitab(rows: Iterable[MitabRow]) -> str:
    mitab = StringIO() # StringIO: set string as object file
    writer = csv.DictWriter(mitab, fieldnames=["#uidA","uidB","author","pmids","sourcedb"], dialect="excel-tab", quoting=csv.QUOTE_NONE, quotechar='')
    writer.writeheader()
    for row in rows:
        writer.writerow({"#uidA":row.uida,"uidB":row.uidb,"author":row.author,"pmids":row.pmids,"sourcedb":row.sourcedb})
        print(row.sourcedb)
    return mitab.getvalue()

