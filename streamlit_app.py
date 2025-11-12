#from typing import Optional
import streamlit as st
import os
import io
import requests
import pandas as pd
import streamlit as st
from urllib.parse import quote
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
from collections import defaultdict

# -------------------------------
# Page Configuration
# -------------------------------
st.set_page_config(
    page_title="ChEMBL Cheminformatics Dashboard",
    layout="wide"
)

# -------------------------------
# Sidebar: Info & Feedback
# -------------------------------
with st.sidebar:
    st.title("🔗 Dashboard Info")
    st.markdown(
        "**Launch URL**  \n"
        "[Open in Streamlit]"
        "(https://chemblcheminformaticsdashboard-cupxnrw6yf56zrsglzc8uk.streamlit.app/)"
    )
    st.markdown("## ℹ️ About This App")
    st.write(
        """
        A cheminformatics dashboard that allows molecular similarity comparisons 
        using ChEMBL data and RDKit fingerprints 🧪🔍  
        - Input a compound using SMILES or ChEMBL ID  
        - Compare it to user-defined compound names  
        - View similarity scores and structural representations  
        Built with ❤️ by Subramanian Ramajayam.
        """
    )
    st.markdown("---")
    st.subheader("💬 Share your feedback")
    feedback_text = st.text_area(
        "Let us know what you liked, what confused you, or new ideas you'd love!"
    )
    if st.button("Submit Feedback"):
        from feedback import collect_feedback, save_feedback
        fb = collect_feedback(feedback_text)
        if fb["status"] == "submitted":
            save_feedback(fb)
            st.success("Thanks for your feedback!")
        else:
            st.warning("Please enter something before submitting.")
    if os.path.exists("feedback_log.csv"):
        df_fb = pd.read_csv(
            "feedback_log.csv",
            header=None,
            names=["Timestamp", "Status", "Feedback"]
        )
        st.markdown("### 📝 Submitted Feedback History")
        st.dataframe(df_fb)
        with open("feedback_log.csv", "rb") as f:
            st.download_button(
                "⬇️ Download Feedback CSV",
                f,
                file_name="feedback_log.csv",
                mime="text/csv"
            )
    else:
        st.info("No feedback submitted yet.")

# -------------------------------
# Initialize ChEMBL Client
# -------------------------------
molecule_client = new_client.molecule
activity_client = new_client.activity

# -------------------------------
# Helper Functions
# -------------------------------
from typing import Optional

def resolve_to_chembl_id(compound_name: str) -> Optional[str]:
    """Query ChEMBL REST API to resolve a name to a ChEMBL ID."""
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
        f"?q={quote(compound_name)}"
    )
    try:
        res = requests.get(url, timeout=10)
        res.raise_for_status()
        mols = res.json().get("molecules", [])
        if mols:
            return mols[0].get("molecule_chembl_id")
    except Exception as e:
        st.error(f"Name resolution error: {e}")
    return None

#def resolve_to_chembl_id(compound_name: str) -> str | None:
#    """Query ChEMBL REST API to resolve a name to a ChEMBL ID."""
 #   url = (
  #      "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
 #       f"?q={quote(compound_name)}"
#    )
#    try:
#        res = requests.get(url, timeout=10)
 #       res.raise_for_status()
 #       mols = res.json().get("molecules", [])
#        if mols:
#            return mols[0].get("molecule_chembl_id")
 #   except Exception as e:
#        st.error(f"Name resolution error: {e}")
 #   return None

@st.cache_data
def load_synonym_map(csv_path: str = "synonyms.csv") -> dict[str, str]:
    """Load brand→generic mappings from CSV and merge with static defaults."""
    # Static fallbacks
    static_map = {
        "restyl": "alprazolam",
        "tylenol": "acetaminophen",
        "ponstan": "mefenamic acid",
    }
    try:
        df = pd.read_csv(csv_path)
        csv_map = {
            row["Brand"].strip().lower(): row["Generic"].strip()
            for _, row in df.iterrows()
        }
    except Exception:
        csv_map = {}
    # CSV entries override static if keys collide
    return {**static_map, **csv_map}

synonym_map = load_synonym_map()

def resolve_synonym(name: str) ->Optional[ str]:
    """Return generic name if known, else the original input."""
    return synonym_map.get(name.strip().lower(), name.strip())

def handle_failed_lookup(name: str) -> Optional[str]:
    """
    Prompt user for generic name when ChEMBL lookup fails.
    Returns the new generic name or None.
    """
    st.warning(f"❌ Couldn't find '{name}'. Provide its generic name:")
    gen = st.text_input("Generic name", key=f"syn-{name}")
    if gen:
        synonym_map[name.lower()] = gen.strip()
        st.success(f"Using '{gen}' for '{name}' now.")
        return gen.strip()
    return None

# -------------------------------
# Tabs Definition
# -------------------------------
tab1, tab2, tab3 = st.tabs([
    "🧪 Compound Explorer",
    "🎯 Bioactivity Explorer",
    "🔍 Similar Compounds"
])

# -------------------------------
# Tab 1: Compound Explorer
# -------------------------------
with tab1:
    st.subheader("Compound Explorer")

    # Input area
    raw_input = st.text_area(
        "Enter compound names (one per line)",
        value="aspirin\nibuprofen\ncaffeine",
        help="Use brand or generic. Mappings will be shown below."
    ).strip()

    # Show original → resolved mapping
    with st.expander("🧾 Generic Name Mappings"):
        mapping = {
            orig: resolve_synonym(orig)
            for orig in raw_input.splitlines() if orig.strip()
        }
        st.table(
            pd.DataFrame.from_dict(
                mapping,
                orient="index",
                columns=["Resolved Generic Name"]
            )
        )

    # Prepare deduplicated compound list
    compounds = [
        resolve_synonym(x) for x in raw_input.splitlines() if x.strip()
    ]
    compounds = list(dict.fromkeys(compounds))

    if not compounds:
        st.warning("⚠️ Please enter at least one valid compound name.")
    else:
        all_data = []
        smiles_dict = {}

        for cmpd in compounds:
            st.markdown(f"---\n#### 🔍 {cmpd}")
            results = molecule_client.filter(
                pref_name__iexact=cmpd
            ).only(
                ["molecule_chembl_id", "molecule_structures", "molecule_properties"]
            )

            if not results:
                cid = resolve_to_chembl_id(cmpd)
                if cid:
                    results = molecule_client.filter(
                        molecule_chembl_id=cid
                    ).only(
                        ["molecule_chembl_id", "molecule_structures", "molecule_properties"]
                    )

            if results:
                mol = results[0]
                chembl_id = mol["molecule_chembl_id"]
                props = mol.get("molecule_properties", {})
                smiles = mol.get("molecule_structures", {}).get("canonical_smiles")
                link = (
                    f"https://www.ebi.ac.uk/chembl/"
                    f"compound_report_card/{chembl_id}/"
                )

                st.markdown(f"**ChEMBL ID:** [{chembl_id}]({link})")

                if smiles:
                    rd_mol = Chem.MolFromSmiles(smiles)
                    st.image(
                        Draw.MolToImage(rd_mol, size=(200, 200)),
                        caption=cmpd
                    )
                    smiles_dict[cmpd] = smiles

                display = {
                    "Compound Name": cmpd,
                    "ChEMBL ID": chembl_id,
                    "SMILES": smiles or "N/A",
                    "Mol Wt": props.get("full_mwt", "N/A"),
                    "AlogP": props.get("alogp", "N/A"),
                    "HBA": props.get("hba", "N/A"),
                    "HBD": props.get("hbd", "N/A"),
                    "RO5 Violations": props.get("num_ro5_violations", "N/A")
                }
                st.dataframe(pd.DataFrame([display]))
                all_data.append(display)
            else:
                fallback = handle_failed_lookup(cmpd)
                if fallback:
                    st.experimental_rerun()

        # Download results
        if all_data:
            df_export = pd.DataFrame(all_data)
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
                df_export.to_excel(
                    writer,
                    index=False,
                    sheet_name="Compound Properties"
                )
            buffer.seek(0)
            st.download_button(
                label="📥 Download All Properties as Excel",
                data=buffer,
                file_name="chembl_compound_data.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

# -------------------------------
# Tab 2: Bioactivity Explorer (Stub)
# -------------------------------
with tab2:
    st.subheader("🎯 Bioactivity Explorer")

    # 1. Accept a compound name or ChEMBL ID
    bio_input = st.text_input(
        "Enter compound name or ChEMBL ID for bioactivity search",
        help="You can type 'aspirin', 'CHEMBL25', etc."
    ).strip()

    if bio_input:
        # 2. Resolve to a ChEMBL ID if needed
        if bio_input.upper().startswith("CHEMBL"):
            chembl_id = bio_input.upper()
        else:
            name = resolve_synonym(bio_input)
            chembl_id = resolve_to_chembl_id(name)

        if chembl_id:
            st.success(f"Fetching activities for **{chembl_id}**")
            # 3. Query the activity client
            acts = activity_client.filter(
                molecule_chembl_id=chembl_id
            ).only([
                "activity_id",
                "assay_chembl_id",
                "standard_type",
                "standard_value",
                "standard_units",
                "pchembl_value",
                "target_chembl_id"
            ])[:50]  # limit to first 50 for speed

            if acts:
                df_acts = pd.DataFrame(acts)
                # Show key columns
                st.dataframe(
                    df_acts[[
                        "activity_id",
                        "assay_chembl_id",
                        "standard_type",
                        "standard_value",
                        "standard_units",
                        "pchembl_value",
                        "target_chembl_id"
                    ]]
                )
                # Optional: plot pChEMBL histogram
                if "pchembl_value" in df_acts.columns:
                    st.markdown("**pChEMBL value distribution**")
                    st.bar_chart(df_acts["pchembl_value"].dropna())
            else:
                st.warning(f"No bioactivity data found for {chembl_id}.")
        else:
            st.error(f"Could not resolve '{bio_input}' to a ChEMBL ID.")

    else:
        st.info("Enter a compound name or CHEMBL ID above to see its bioactivity.")



# -------------------------------
# Tab 3: Similar Compounds (Stub)
# -------------------------------

    from rdkit.Chem import AllChem, DataStructs
    from utils import get_query_smiles  # your helper to fetch SMILES from ID
with tab3:
    st.subheader("🔍 Similar Compounds Explorer")

    # 1. Query compound input
    query_input = st.text_input(
        "Enter query compound (SMILES or ChEMBL ID)",
        help="Type a SMILES string or a ChEMBL ID like CHEMBL25"
    )

    if query_input:
        # convert ID → SMILES if needed
        query_smiles = (
            get_query_smiles(query_input)
            if not query_input.startswith(("C", "c")) or " " in query_input
            else query_input
        )
        try:
            query_mol = Chem.MolFromSmiles(query_smiles)
            query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, 2048)
        except Exception as e:
            st.error(f"Invalid query SMILES: {e}")
            query_fp = None

        if query_fp:
            st.success(f"Query SMILES ready: `{query_smiles}`")
            
            # 2. Upload CSV of target compounds
            upload = st.file_uploader(
                "Upload CSV with a 'SMILES' column to compare",
                  type="csv"
            )

            if upload:
                df_targets = pd.read_csv(upload)
                if "SMILES" not in df_targets.columns:
                    st.error("CSV must contain a 'SMILES' column.")
                else:
                    sims = []
                    for sm in df_targets["SMILES"].dropna().astype(str):
                        try:
                            mol = Chem.MolFromSmiles(sm)
                            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
                            sims.append({"SMILES": sm, "Similarity": sim})
                        except:
                            continue

                    if sims:
                        sim_df = (
                            pd.DataFrame(sims)
                              .sort_values("Similarity", ascending=False)
                              .reset_index(drop=True)
                        )
                        st.markdown("**Similarity Results**")
                        st.dataframe(sim_df)

                        # 3. Download / histogram
                        st.download_button(
                            "⬇️ Download Similarity CSV",
                            sim_df.to_csv(index=False).encode("utf-8"),
                            file_name="similarity_results.csv",
                            mime="text/csv"
                        )
                        st.bar_chart(sim_df["Similarity"])
                    else:
                        st.warning("No valid SMILES found in the uploaded file.")
    else:
        st.info("Enter a query compound to start similarity search.")

