import csv
from dataclasses import dataclass
from io import StringIO
from typing import Iterable

@dataclass
class Row:
    complex: str
    uniprot_ids: list[str]

@dataclass
class MitabRow:
    uida: str
    uidb: str
    

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

def convert_to_mitab(rows: list[Row]) -> list[MitabRow]:
    for row in rows:
        for id in row.uniprot_ids:
            yield MitabRow(
                uida = row.complex, 
                uidb = id
            )

def serialize_to_mitab(rows: Iterable[MitabRow]) -> str:
    mitab = StringIO()
    writer = csv.DictWriter(mitab, fieldnames=["#uidA","uidB"], dialect="excel-tab")
    writer.writeheader()
    for row in rows:
        writer.writerow({"#uidA":row.uida,"uidB":row.uidb})
    return mitab.getvalue()
