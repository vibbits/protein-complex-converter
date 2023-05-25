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
	method: str
	first_author_name: list[str]
	pubmed_ids: list[str]
	taxa: str
	taxb: str
	interactiontype: str
	source: str
	interactionidentifiers: list[str]
	expansion: str
	bio_role_A: str
	bio_role_B: str
	exp_role_A: str
	exp_role_B: str 
	interactortype_A: str
	xref: list[str]
	taxid: str

@dataclass
class MitabRow:
	uida: str
	uidb: str
	method: str
	author: str
	pmids: str
	taxa: str 
	taxb: str
	interactiontype: str
	sourcedb: str
	interactionidentifiers: str
	expansion: str
	bio_role_A: str
	bio_role_B: str
	exp_role_A: str
	exp_role_B: str
	interactortype_A: str
	xref: str
	host_org_taxid: str
    
    
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

# Extract the interaction identifiers
def extract_interaction_identifiers(column1:str, column2:list[str]) -> str:
	interaction_identifiers = []
	if column1 != '-':
		interaction_identifiers.append(column1)
	for ref in column2:
		ref = ref.split('(')[0]
		if not ref.startswith('pubmed'):
			interaction_identifiers.append(ref)
	return '|'.join(interaction_identifiers)

# Extract the interactor types 

# Extract the taxonomy identifier of the interaction, host organism
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
		interactionidentifiers = extract_interaction_identifiers(row["Experimental evidence"], list(row["Cross references"].split('|')))
		rows.append(Row(
			complex=row["#Complex ac"],
			uniprot_ids=[extract_id(name) for name in row[
				"Identifiers (and stoichiometry) of molecules in complex"
			].split("|")],
			method = 'psi-mi:"MI:0000"(NA)',
			first_author_name = '|'.join(first_author_name),
			pubmed_ids= '|'.join(pubmed_ids),
			taxa = "NA",
			taxb = "NA",
			interactiontype = "NA",
			source = row["Source"],
			interactionidentifiers = interactionidentifiers,
			expansion = "bipartite",
			bio_role_A = 'psi-mi:"MI:0000"(unspecified)',
			bio_role_B = 'psi-mi:"MI:0000"(unspecified)',
			exp_role_A = 'psi-mi:"MI:0000"(unspecified)',
			exp_role_B = 'psi-mi:"MI:0000"(unspecified)',
			interactortype_A = 'psi-mi:"MI:0315"(protein complex)',
			xref = row["Go Annotations"],
			taxid = extract_taxid(row["Taxonomy identifier"])
		))
		print(interactionidentifiers)
	return rows

# Convert the complex ids and uniprot ids into mitab format
def convert_to_mitab(rows: list[Row]) -> list[MitabRow]:
	for row in rows:
		for id in row.uniprot_ids:
			yield MitabRow(
				uida = row.complex, 
				uidb = id,
				method = row.method,
				author = row.first_author_name,
				pmids = row.pubmed_ids,
				taxa = row.taxa,
				taxb = row.taxb,
				interactiontype = row.interactiontype,
				sourcedb = row.source,
				interactionidentifiers = row.interactionidentifiers,
				expansion = row.expansion,
				bio_role_A = row.bio_role_A,
				bio_role_B = row.bio_role_B,
				exp_role_A = row.exp_role_A,
				exp_role_B = row.exp_role_B,
				interactortype_A = row.interactortype_A,
				xref = row.xref,
				host_org_taxid = row.taxid
			)

# Serialize the mitab results
def serialize_to_mitab(rows: Iterable[MitabRow]) -> str:
	mitab = StringIO() # StringIO: set string as object file
	writer = csv.DictWriter(mitab, fieldnames=["#uidA","uidB","method","author","pmids","taxa","taxb","interactionType","sourcedb",
	"interactionIdentifier","expansion","biological_role_A","biological_role_B","experimental_role_A","experimental_role_B",
	"interactor_type_A","xrefs_Interaction","Host_organism_taxid"], 
	dialect="excel-tab", quoting=csv.QUOTE_NONE, quotechar='')
	writer.writeheader()
	for row in rows:
		writer.writerow({"#uidA":row.uida,"uidB":row.uidb,"method":row.method,"author":row.author,"pmids":row.pmids,"taxa":row.taxa,"taxb":row.taxb,
		"interactionType":row.interactiontype,"sourcedb":row.sourcedb,"interactionIdentifier":row.interactionidentifiers,"expansion":row.expansion,
		"biological_role_A":row.bio_role_A,"biological_role_B":row.bio_role_B,"experimental_role_A":row.bio_role_A,"experimental_role_B":row.exp_role_B,
		"interactor_type_A":row.interactortype_A,"xrefs_Interaction":row.xref,"Host_organism_taxid":row.host_org_taxid})
	return mitab.getvalue()

