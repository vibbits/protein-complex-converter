import csv
from dataclasses import dataclass

@dataclass
class Row:
    complex: str
    uniprot_ids: list[str]

def extract_id(name: str) -> str:
    return name.split("(")[0]

def parse_complex_tab(ct: str) -> None:
    reader = csv.DictReader(ct.splitlines(), delimiter="\t")
    return [
        Row(
            complex=row["#Complex ac"],
            uniprot_ids=[extract_id(name) for name in row[
                "Identifiers (and stoichiometry) of molecules in complex"
            ].split("|")],
        )
        for row in reader
    ]
