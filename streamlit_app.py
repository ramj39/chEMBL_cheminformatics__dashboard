import streamlit as st
st.sidebar.title("🔗 Dashboard Info")
st.sidebar.markdown("""
**Launch URL**  
[Open in Streamlit](https://chemblcheminformaticsdashboard-cupxnrw6yf56zrsglzc8uk.streamlit.app/)
""")
import requests
from urllib.parse import quote
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
import pandas as pd
import io
from utils import get_query_smiles
# 👋 Sidebar About Section
with st.sidebar:
    st.markdown("## ℹ️ About This App")
    st.write("""
    A cheminformatics dashboard that allows molecular similarity comparisons 
    using ChEMBL data and RDKit fingerprints 🧪🔍
    - Input a compound using SMILES or ChEMBL ID
    - Compare it to user-defined compound names
    - View similarity scores and structural representations
    Built with ❤️ by [subramanianRamajayam] using Streamlit, RDKit, and ChEMBL.
""")  
import streamlit as st
from feedback import collect_feedback, save_feedback

st.subheader("💬 Share your feedback")
user_input = st.text_area("What did you like? What can we improve?", "")

if st.button("Submit Feedback"):
    data = collect_feedback(user_input)
    if data["status"] == "submitted":
        save_feedback(data)
        st.success("Thanks for your feedback!")
    else:
        st.warning("Please enter something before submitting.")
    
# Initialize ChEMBL clients
molecule_client = new_client.molecule
activity_client = new_client.activity

# Streamlit page setup
st.set_page_config(page_title="ChEMBL Dashboard", layout="wide")
st.title("🔬 ChEMBL Cheminformatics Dashboard")
st.caption("Explore compound data, bioactivity, and chemical similarity")

# Helper function to resolve compound names via API
def resolve_to_chembl_id(compound_name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={quote(compound_name)}"
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        data = res.json()
        molecules = data.get("molecules", [])
        if molecules:
            return molecules[0].get("molecule_chembl_id")
    except requests.exceptions.RequestException as e:
        st.error(f"⚠️ Failed to connect to ChEMBL for name resolution: {e}")
    return None

# Layout with three tabs
tab1, tab2, tab3 = st.tabs(["🧪 Compound Explorer", "🎯 Bioactivity Explorer", "🔍 Similar Compounds"])

# -------------------------------
# 🧪 Tab 1: Compound Explorer
# -------------------------------
with tab1:
    compound_list_raw = st.text_area(
        "Enter compound names (one per line)",
        value="aspirin\nibuprofen\ncaffeine",
        help="Use IUPAC or common names. Each name should be on a separate line."
    ).strip()

    compound_list = [name.strip() for name in compound_list_raw.splitlines() if name.strip()]
    compound_list = list(dict.fromkeys(compound_list))  # Remove duplicates

    if compound_list:
        all_data = []
        smiles_dict = {}

        for compound_name in compound_list:
            st.subheader(f"🔍 {compound_name}")
            results = molecule_client.filter(pref_name__iexact=compound_name).only(
                ['molecule_chembl_id', 'molecule_structures', 'molecule_properties']
            )

            if not results:
                resolved_id = resolve_to_chembl_id(compound_name)
                if resolved_id:
                    results = molecule_client.filter(molecule_chembl_id=resolved_id).only(
                        ['molecule_chembl_id', 'molecule_structures', 'molecule_properties']
                    )

            if results:
                mol_data = results[0]
                chembl_id = mol_data['molecule_chembl_id']
                props = mol_data.get('molecule_properties', {})
                smiles = mol_data.get('molecule_structures', {}).get('canonical_smiles', None)

                #st.markdown(f"**ChEMBL ID:** `{chembl_id}`")
                chembl_link = f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"
                st.markdown(f"**ChEMBL ID:** [{chembl_id}]({chembl_link})")

                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    st.image(Draw.MolToImage(mol, size=(250, 250)), caption=compound_name)
                    smiles_dict[compound_name] = smiles

                display_props = {
                    'Compound Name': compound_name,
                    'ChEMBL ID': chembl_id,
                    'SMILES': smiles,
                    'Molecular Weight': props.get('full_mwt', 'N/A'),
                    'AlogP': props.get('alogp', 'N/A'),
                    'HBA': props.get('hba', 'N/A'),
                    'HBD': props.get('hbd', 'N/A'),
                    'RO5 Violations': props.get('num_ro5_violations', 'N/A')
                }

                st.dataframe(pd.DataFrame([display_props]))
                all_data.append(display_props)
            else:
                st.warning(f"No compound found for '{compound_name}'.")

        if all_data:
            df_export = pd.DataFrame(all_data)
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                df_export.to_excel(writer, index=False, sheet_name='Compound Properties')
            buffer.seek(0)

            st.download_button(
                label="📥 Download All Properties as Excel",
                data=buffer,
                file_name="chembl_compound_data.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    else:
        st.warning("⚠️ Please enter at least one valid compound name.")

# -------------------------------
# 🎯 Tab 2: Bioactivity Explorer
# -------------------------------
with tab2:
    chembl_id_input = st.text_input(
        "Enter ChEMBL ID (e.g., CHEMBL25 for aspirin)",
        help="The ID must start with 'CHEMBL' followed by digits."
    ).strip().upper()

    if chembl_id_input:
        if chembl_id_input.startswith("CHEMBL") and chembl_id_input[6:].isdigit():
            st.subheader(f"🔬 Bioactivity Data for {chembl_id_input}")
            activities = activity_client.filter(molecule_chembl_id=chembl_id_input).only(
                ['target_chembl_id', 'target_organism', 'standard_type', 'standard_value', 'standard_units']
            )
            if activities:
                df_bio = pd.DataFrame(activities)
                st.dataframe(df_bio)
                bio_buffer = io.BytesIO()
                with pd.ExcelWriter(bio_buffer, engine='xlsxwriter') as writer:
                    df_bio.to_excel(writer, index=False, sheet_name='Bioactivity')
                bio_buffer.seek(0)

                st.download_button(
                    label="📥 Download Bioactivity Data",
                    data=bio_buffer,
                    file_name=f"{chembl_id_input}_bioactivity.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
            else:
                st.warning("No bioactivity data found for this compound.")
        else:
            st.error("❌ Invalid ChEMBL ID format. It should look like 'CHEMBL25'.")

## 🔍 Tab 3: Similar Compounds
# -------------------------------
with tab3:
    st.subheader("🔬 Compare to Named Compounds")

    # 1️⃣ Choose input mode
    input_mode = st.radio(
        "Select query type",
        ["SMILES", "ChEMBL ID"],
        key="query_mode",
        help="Choose SMILES directly or look up by ChEMBL ID"
    )

    query_mol = None
    query_smiles = ""

    if input_mode == "SMILES":
        query_smiles = st.text_input("Enter SMILES 🧪", key="smiles_input")
        if query_smiles:
            query_mol = Chem.MolFromSmiles(query_smiles)
    elif input_mode == "ChEMBL ID":
        chembl_id = st.text_input("Enter ChEMBL ID 🔍", key="chembl_input")
        if chembl_id:
            try:
                res = molecule_client.filter(molecule_chembl_id=chembl_id).only(['molecule_structures'])
                query_smiles = res[0].get('molecule_structures', {}).get('canonical_smiles', None) if res else None
                query_mol = Chem.MolFromSmiles(query_smiles) if query_smiles else None
            except Exception:
                query_mol = None

    if query_mol:
        try:
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
            st.success("✅ Query molecule parsed successfully.")

            # 2️⃣ Example compound presets
            preset = st.selectbox(
                "Load example compounds (optional)",
                ["None", "Painkillers", "Stimulants", "Antibiotics"],
                index=0
            )

            preset_dict = {
                "Painkillers": "aspirin, ibuprofen, acetaminophen",
                "Stimulants": "caffeine, amphetamine, modafinil",
                "Antibiotics": "amoxicillin, tetracycline, ciprofloxacin"
            }

            default_text = preset_dict.get(preset, "")
            compound_input = st.text_area(
                "Enter compound names (comma or newline-separated)",
                value=default_text,
                help="Examples: aspirin, ibuprofen, caffeine"
            )

            if st.checkbox("🔄 Clear compound input"):
                compound_input = ""

            raw_names = compound_input.replace("\n", ",")
            sample_names = [name.strip() for name in raw_names.split(",") if name.strip()]
            sim_data = []
            failed_names = []

            if sample_names:
                for name in sample_names:
                    try:
                        res = molecule_client.filter(pref_name__iexact=name).only(['molecule_chembl_id', 'molecule_structures'])
                        smiles = res[0].get('molecule_structures', {}).get('canonical_smiles', None) if res else None
                        mol = Chem.MolFromSmiles(smiles) if smiles else None

                        if mol:
                            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
                            sim_data.append({
                                "Name": name,
                                "SMILES": smiles,
                                "Similarity": round(sim, 3)
                            })
                        else:
                            failed_names.append(name)
                    except Exception:
                        failed_names.append(name)

                if sim_data:
                    df_sim = pd.DataFrame(sim_data).sort_values(by="Similarity", ascending=False)
                    st.success(f"✅ Compared {len(sim_data)} compound(s) successfully.")
                    st.dataframe(df_sim)

                    # Download button
                    csv = df_sim.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="📥 Download Results as CSV",
                        data=csv,
                        file_name="similarity_results.csv",
                        mime="text/csv"
                    )

                if failed_names:
                    st.error(f"❌ Could not process the following compounds: {', '.join(failed_names)}")
            else:
                st.warning("⚠️ Please enter at least one compound name to compare.")

        except Exception as e:
            st.error(f"❌ SMILES parsing failed: {str(e)}")
    else:
        st.info("Enter a valid query above to begin comparison.")
