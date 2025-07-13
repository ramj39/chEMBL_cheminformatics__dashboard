🚀 [Launch Dashboard](https://chemblcheminformaticsdashboard-cupxnrw6yf56zrsglzc8uk.streamlit.app/)
# 🔬 ChEMBL Cheminformatics Dashboard

A streamlined, interactive dashboard built with **Streamlit** to explore compound information, bioactivity data, and molecular similarity using the ChEMBL database and cheminformatics tools.

## 🌟 Features

- **🧪 Compound Explorer**  
  Input multiple compound names to retrieve their ChEMBL IDs, properties, and molecular structures.  
  📥 Download all properties as an Excel file.

- **🎯 Bioactivity Explorer**  
  Enter a ChEMBL ID to explore bioactivity data across different targets and organisms.  
  📥 Download full bioactivity data as an Excel file.

- **🔍 Similar Compounds**  
  Input a SMILES string to compute molecular similarity scores against a sample set of known compounds using **RDKit’s Morgan fingerprints**.

## 🚀 Deployment

This app runs seamlessly on **Streamlit Cloud**. To deploy your own version:

1. Fork or clone this repository.
2. Ensure `requirements.txt` includes:
