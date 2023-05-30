import csv, time, re
from dataclasses import dataclass
from io import StringIO
from typing import Iterable
from bs4 import BeautifulSoup
from Bio import Entrez
from datetime import date
import requests as req
import xml.etree.ElementTree as et

# Define classes for the data to use
@dataclass
class Row:
	complex: str
	uniprot_ids: list[str]
	altA: str
	altB: str
	aliasA: str
	aliasB: str
	method: str
	first_author_name: list[str]
	pubmed_ids: list[str]
	taxa: str
	taxb: str
	interactiontype: str
	source: str
	interactionidentifiers: list[str]
	confidence: list[str]
	expansion: str
	bio_role_A: str
	bio_role_B: str
	exp_role_A: str
	exp_role_B: str 
	interactortype_A: str
	interactortype_B: list[str]
	xref_A: list[str]
	xref_B: list[str]
	xref_interaction: list[str]
	annotation_A: list[str]
	annotation_B: list[str]
	annotation_interaction: list[str]
	taxid: str
	params_interaction: str
	creation_date: str
	update_date: str
	checksum_A: str
	checksum_B: str
	checksum_interaction: str
	negative: str
	feature_A: str
	feature_B: str
	stoichiometry_A: str
	stoichiometry_B: list[int]
	identification_method_A: str
	identification_method_B: str

@dataclass
class MitabRow:
	uida: str
	uidb: str
	altA: str
	altB: str
	aliasA: str
	aliasB: str
	method: str
	author: str
	pmids: str
	taxa: str 
	taxb: str
	interactiontype: str
	sourcedb: str
	interactionidentifiers: str
	confidence: str
	expansion: str
	bio_role_A: str
	bio_role_B: str
	exp_role_A: str
	exp_role_B: str
	interactortype_A: str
	interactortype_B: str
	xref_A: str
	xref_B: str
	xref_interaction: str
	annotation_A: str
	annotation_B: str
	annotation_interaction: str
	host_org_taxid: str
	params_interaction: str
	creation_date: str
	update_date: str
	checksum_A: str
	checksum_B: str
	checksum_interaction: str
	negative: str
	feature_A: str
	feature_B: str
	stoichiometry_A: str
	stoichiometry_B: int
	identification_method_A: str
	identification_method_B: str
    
    
# Extract the identifiers of molecules in complex
def extract_id(name: str) -> str:
	return name.split("(")[0]

# Extract pubmed identifiers
def extract_pubmed_id(crossref:str) -> str:
	return crossref.split("(")[0] 

# Get pubmed xml for all pubmed ids at once 
def get_pubmedArticles(pubmed_ids):
	pubmed_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&rettype=abstract"
	response = req.get(pubmed_base_url+'&id='+','.join([s.removeprefix('pubmed:') for s in pubmed_ids]))
	print(pubmed_base_url+'&id='+','.join([s.removeprefix('pubmed:') for s in pubmed_ids]))
	xml = response.text
	pubmedArticleSet = et.fromstring(xml)
	return xml

# Parse pubmed xml to get first author and publication year
def extract_first_author(xml:str) -> str:
	soup = BeautifulSoup(xml, 'xml')
	articles = soup.find_all('PubmedArticle')
	author_names = []
	for article in articles:
		author_list = article.find_all('Author')
		first_author = author_list[0]
		last_name = first_author.find('LastName').text
		publication_year = article.find('PubDate').find('Year').text
		last_name_with_year = last_name + " " + "et al. " + "(" + publication_year + ")"
		author_names.append(last_name_with_year)
	return '|'.join(author_names)

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
def extract_interactor_type(uniprot_id:str) -> str:
	if re.search("^CHEBI",uniprot_id):
		return 'psi-mi:"MI:0328"(small molecule)'
	elif re.search("^URS",uniprot_id):
		return 'psi-mi:"MI:0320"(RNA)'
	elif re.search("^CPX",uniprot_id):
		return'psi-mi:"MI:0315"(protein complex)'
	else:
		return'psi-mi:"MI:0326"(protein)'

# Get taxonomy xml for all taxids at once 
def get_taxid_xml(taxids):
	taxonomy_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&rettype=abstract"
	response = req.get(taxonomy_base_url+'&id='+ taxids)
	print(taxonomy_base_url+'&id='+ taxids)
	xml = response.text
	return xml

# Parse taxonomy xml to get the scientific name of the organism starting from NCBI taxid
def extract_organism(xml:str) -> str:
	soup = BeautifulSoup(xml, 'xml')
	organisms = soup.find_all('Taxon')
	scientific_name = ""

	for organism in organisms:
		for child in organism.find_all('Taxon'):
			child.decompose()
		
		scientific_name_tag = organism.find('ScientificName')
		if scientific_name_tag:
			scientific_name = scientific_name_tag.text.strip()
	return scientific_name
	

# Get the current date
def get_current_date():
	today = date.today()
	current_date = today.strftime("%Y/%m/%d")
	return current_date

# Extract the stoichiometry of each molecule in the complex
def extract_stoichiometry(name:str) -> int:
	return name.split("(")[1].rstrip(")")

# Get the complex, uniprot and pubmed identifiers and parse first author from pmid
def parse_complex_tab(ct: str) -> None:
	reader = csv.DictReader(ct.splitlines(), delimiter="\t")
	rows = []
	for row in reader:
		first_author_name = []
		pubmed_ids = [extract_pubmed_id(crossref) for crossref in row["Cross references"
		].split("|") if 'pubmed:' in crossref]
		first_author_name.append(extract_first_author(get_pubmedArticles(pubmed_ids)))
		uniprot_ids = [extract_id(name) for name in row["Identifiers (and stoichiometry) of molecules in complex"
		].split("|")]
		interactionidentifiers = extract_interaction_identifiers(row["Experimental evidence"], list(row["Cross references"].split('|')))
		interactortype_B = [extract_interactor_type(uniprot_id) for uniprot_id in uniprot_ids]
		stoichiometry_B = [int(extract_stoichiometry(name)) for name in row["Identifiers (and stoichiometry) of molecules in complex"]
		.split("|")]
		taxids = row["Taxonomy identifier"]
		taxid = "taxid:"+taxids+"("+extract_organism(get_taxid_xml(taxids))+")"

		rows.append(Row(
			complex=row["#Complex ac"],
			uniprot_ids=uniprot_ids,
			altA = "-",
			altB = "-",
			aliasA = "-",
			aliasB = "-",
			method = 'psi-mi:"MI:0686"(unspecified method)',
			first_author_name = '|'.join(first_author_name),
			pubmed_ids= '|'.join(pubmed_ids),
			taxa = "-",
			taxb = taxid,
			interactiontype = "-",
			source = row["Source"],
			interactionidentifiers = interactionidentifiers,
			confidence = '-',
			expansion = "bipartite",
			bio_role_A = 'psi-mi:"MI:0000"(unspecified)',
			bio_role_B = 'psi-mi:"MI:0000"(unspecified)',
			exp_role_A = 'psi-mi:"MI:0000"(unspecified)',
			exp_role_B = 'psi-mi:"MI:0000"(unspecified)',
			interactortype_A = 'psi-mi:"MI:0315"(protein complex)',
			interactortype_B = interactortype_B,
			xref_A = "-",
			xref_B = "-",
			xref_interaction = row["Go Annotations"],
			annotation_A = "-",
			annotation_B = "-",
			annotation_interaction = "-",
			taxid = taxid,
			params_interaction = "-",
			creation_date = get_current_date(),
			update_date = get_current_date(),
			checksum_A = '-',
			checksum_B = '-',
			checksum_interaction = '-',
			negative = "False",
			feature_A = "-",
			feature_B = "-",
			stoichiometry_A = "-",
			stoichiometry_B = stoichiometry_B,
			identification_method_A = "-",
			identification_method_B = "-"
		))
	return rows

# Convert the complex ids and uniprot ids into mitab format
def convert_to_mitab(rows: list[Row]) -> list[MitabRow]:
	for row in rows:
		for id, interactortype_B, stoichiometry_B in zip(row.uniprot_ids, row.interactortype_B, row.stoichiometry_B):
			yield MitabRow(
				uida = row.complex, 
				uidb = id,
				altA = row.altA,
				altB = row.altB,
				aliasA = row.aliasA,
				aliasB = row.aliasB,
				method = row.method,
				author = row.first_author_name,
				pmids = row.pubmed_ids,
				taxa = row.taxa,
				taxb = row.taxb,
				interactiontype = row.interactiontype,
				sourcedb = row.source,
				interactionidentifiers = row.interactionidentifiers,
				confidence = row.confidence,
				expansion = row.expansion,
				bio_role_A = row.bio_role_A,
				bio_role_B = row.bio_role_B,
				exp_role_A = row.exp_role_A,
				exp_role_B = row.exp_role_B,
				interactortype_A = row.interactortype_A,
				interactortype_B = interactortype_B,
				xref_A = row.xref_A,
				xref_B = row.xref_B,
				xref_interaction = row.xref_interaction,
				annotation_A = row.annotation_A,
				annotation_B = row.annotation_B,
				annotation_interaction = row.annotation_interaction,
				host_org_taxid = row.taxid,
				params_interaction = row.params_interaction,
				creation_date = row.creation_date,
				update_date = row.update_date,
				checksum_A = row.checksum_A,
				checksum_B = row.checksum_B,
				checksum_interaction = row.checksum_interaction,
				negative = row.negative,
				feature_A = row.feature_A,
				feature_B = row.feature_B,
				stoichiometry_A = row.stoichiometry_A,
				stoichiometry_B = stoichiometry_B,
				identification_method_A = row.identification_method_A,
				identification_method_B = row.identification_method_B
			)

# Serialize the mitab results
def serialize_to_mitab(rows: Iterable[MitabRow]) -> str:
	mitab = StringIO() # StringIO: set string as object file
	writer = csv.DictWriter(mitab, fieldnames=["#uidA","uidB","altA","altB","aliasA","aliasB","method","author","pmids","taxa","taxb","interactionType","sourcedb",
	"interactionIdentifier","confidence","expansion","biological_role_A","biological_role_B","experimental_role_A","experimental_role_B",
	"interactor_type_A","interactor_type_B","xrefs_A","xrefs_B","xrefs_Interaction","Annotations_A","Annotations_B","Annotations_Interaction",
	"Host_organism_taxid","parameters_Interaction","Creation_date","Update_date", "Checksum_A","Checksum_B","Checksum_Interaction","Negative",
	"Features_A","Features_B","Stoichiometry_A","Stoichiometry_B","Identification_method_A","Identification_method_B"], 
	dialect="excel-tab", quoting=csv.QUOTE_NONE, quotechar='')
	writer.writeheader()
	for row in rows:
		writer.writerow({"#uidA":row.uida,"uidB":row.uidb,"altA":row.altA,"altB":row.altB,"aliasA":row.aliasA,"aliasB":row.aliasB,"method":row.method,
		"author":row.author,"pmids":row.pmids,"taxa":row.taxa,"taxb":row.taxb,"interactionType":row.interactiontype,"sourcedb":row.sourcedb,
		"interactionIdentifier":row.interactionidentifiers,"confidence":row.confidence,"expansion":row.expansion,"biological_role_A":row.bio_role_A,
		"biological_role_B":row.bio_role_B,"experimental_role_A":row.bio_role_A,"experimental_role_B":row.exp_role_B,"interactor_type_A":row.interactortype_A,
		"interactor_type_B":row.interactortype_B,"xrefs_A":row.xref_A,"xrefs_B":row.xref_B,"xrefs_Interaction":row.xref_interaction,"Annotations_A":row.annotation_A,
		"Annotations_B":row.annotation_B,"Annotations_Interaction":row.annotation_interaction,"Host_organism_taxid":row.host_org_taxid,"parameters_Interaction":row.params_interaction,
		"Creation_date":row.creation_date,"Update_date":row.update_date,"Checksum_A":row.checksum_A,"Checksum_B":row.checksum_B,"Checksum_Interaction":row.checksum_interaction,
		"Negative":row.negative,"Features_A":row.feature_A,"Features_B":row.feature_B,"Stoichiometry_A":row.stoichiometry_A,"Stoichiometry_B":row.stoichiometry_B,
		"Identification_method_A":row.identification_method_A,"Identification_method_B":row.identification_method_B})
	return mitab.getvalue()

