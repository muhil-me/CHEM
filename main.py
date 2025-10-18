import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st

st.title("Molecular visualiser")
st.write(" This is a 3D Molecule Viewer")


conn = sqlite3.connect("data.db")  
conn.row_factory = sqlite3.Row     
cursor = conn.cursor()

compound_name = st.text_input("Enter compound name: ")

if st.button("Generate 3D Structure"):
    try:
        query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(?)"
        cursor.execute(query, (compound_name,))
        result = cursor.fetchone()
        
        if result:
            smiles = result['smiles']
            st.subheader("Compound Information")
            st.write(f"**Molecular Formula:** {result['formula']}")
            st.write(f"**Molecular Weight:** {result['molecular_weight']}")
            st.write(f"**IUPAC Name:** {result['iupac_name']}")
            st.write(f"**SMILES:** {result['smiles']}")

            if '.' in smiles or '+' in smiles or '-' in smiles:
                st.write("⚠️ 3D structure not available for ionic compounds like salts.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)
                    
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(mol), 'mol')
                view.setStyle({'stick': {}})
                view.zoomTo()
        
                viewer_html = view._make_html()
                st.components.v1.html(viewer_html, height=450)
        

        else:
            st.error("Compound not found in your local database")
            smiles = None
            
    except Exception as e:
        st.error(f"Error: {e}")
        
    conn.close()
