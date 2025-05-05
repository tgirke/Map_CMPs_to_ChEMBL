# Mapping of Selleck and LATCA Compounds to ChEMBL

## Selleck
1. Map CAS IDs to DrugBank
2. Map CAS IDs to ChEBI (to get those missing in DrugBank)
3. Map CAS IDs to PubChem mappings (to get those missing in DrugBank/ChEBI) are obtained with webchem. Subsequently, fingerprint/atompair similarities are from ChemmineR are used as secondary query to resove ambiguities.
4. Map CAS IDs obtained from steps 1-3 (DrugBank/ChEBI/PubChem) to ChEMBL
5. Retrieve annotations and physicochemical properties of compounds from various sources including: ChEMBL, ChemmineR, CDK and PubChem
6. Add parent molecule hierarchies from ChEMBL to master table (parent molecule refers to the core structure of a compound from which other forms, like salts and hydrates are derived).

## LATCA
1. Map compound names to ChEMBL (additional InChI mappings didn't retrieve additional mappings)
2. Retrieve annotations and physicochemical properties of compounds from various sources including: ChEMBL, ChemmineR, CDK and PubChem

