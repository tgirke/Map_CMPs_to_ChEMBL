# Map_CMPs_to_ChEMBL

## Selleck
1. Get for Selleck CAS IDs annotations from DrugBank
2. Get CAS IDs from ChEBI that are missing in DrugBank
3. Remaining CAS to PubChem mappings are obtained with webchem where fingerprint/atompair similarity from ChemmineR is used as secondary query to resove ambiguities
4. Use UniChem to map CAS IDs from DrugBank/ChEBI/PubChem to ChEMBL
5. Obtain annotations and physicochemical properties of compounds from various sources including: ChEMBL, ChemmineR, CDK and PubChem
6. Add parent molecule hierarchies from ChEMBL

## LATCA
