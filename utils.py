# utils.py

import streamlit as st

def get_query_smiles(input_mode: str, molecule_client) -> str:
    if input_mode == "SMILES":
        return st.text_input(
            "Enter SMILES to find similar compounds",
            help="Paste a valid SMILES string"
        ).strip()

    elif input_mode == "ChEMBL ID":
        chembl_id = st.text_input(
            "Enter ChEMBL ID to fetch compound",
            help="Example: CHEMBL25"
        ).strip().upper()

        if chembl_id:
            res = molecule_client.filter(molecule_chembl_id=chembl_id).only(['molecule_structures'])
            smiles = res[0].get('molecule_structures', {}).get('canonical_smiles', "")
            if smiles:
                st.success(f"✅ SMILES retrieved: `{smiles}`")
            else:
                st.error("❌ No SMILES found for this ChEMBL ID.")
            return smiles

    return ""
