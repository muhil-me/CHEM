from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
import pubchempy as pcp
import time
import psycopg2
from psycopg2.extras import RealDictCursor

st.title("Molecular visualiser")
st.markdown("""
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    [data-testid="stToolbar"] {visibility: hidden;}
    .stDeployButton {visibility: hidden;}
    </style>
    """, unsafe_allow_html=True)


DATABASE_URL = st.secrets["database"]["url"]

def get_db_connection():
    return psycopg2.connect(DATABASE_URL)

def insert_compound(name, formula, weight, iupac, smiles):
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        query = """
        INSERT INTO compounds (name, formula, molecular_weight, iupac_name, smiles)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (name) DO NOTHING
        """
        
        cursor.execute(query, (name, formula, weight, iupac, smiles))
        conn.commit()
        cursor.close()
        conn.close()
        return True
    except Exception as e:
        st.error(f"Error inserting compound: {e}")
        return False

def get_compound_suggestions(search_term, limit=10):
    """Get compound suggestions from database matching the search term"""
    if not search_term:
        # Return initial set when empty
        try:
            conn = get_db_connection()
            cursor = conn.cursor()
            
            query = """
            SELECT DISTINCT name FROM compounds 
            ORDER BY name
            LIMIT %s
            """
            
            cursor.execute(query, (limit,))
            results = cursor.fetchall()
            cursor.close()
            conn.close()
            
            return [row[0] for row in results]
        except Exception as e:
            return []
    
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        
        query = """
        SELECT DISTINCT name FROM compounds 
        WHERE LOWER(name) LIKE LOWER(%s)
        ORDER BY name
        LIMIT %s
        """
        
        cursor.execute(query, (f'%{search_term}%', limit))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        
        return [row[0] for row in results]
    except Exception as e:
        return []

# Initialize session state
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""

# Text input for dynamic searching
search_query = st.text_input(
    "Search compound name:",
    value=st.session_state.search_query,
    key="search_box"
)

# Get dynamic suggestions based on search
suggestions = get_compound_suggestions(search_query, limit=10)

# Dropdown with dynamic suggestions
compound_name = st.selectbox(
    "Select compound:",
    options=[""] + suggestions,
    index=0,
    key="compound_select"
)

if st.button("Generate 3D Structure"):
    if not compound_name:
        st.warning("Please enter a compound name")
    else:
        try:
            # Connect to Neon database
            conn = get_db_connection()
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Search in database
            query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(%s)"
            cursor.execute(query, (compound_name,))
            result = cursor.fetchone()
            
            if result:
                smiles = result['smiles']
                st.subheader("Compound Information")
                st.write(f"**Molecular Formula:** {result['formula']}")
                st.write(f"**Molecular Weight:** {result['molecular_weight']}")
                st.write(f"**IUPAC Name:** {result['iupac_name']}")
                st.write(f"**SMILES:** {result['smiles']}")
                st.info("✅ Data from local database")
                
                if '.' in smiles or '+' in smiles or '-' in smiles:
                    st.write("3D structure not available for ionic compounds like salts.")
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
                # Not found in database, query PubChem
                with st.spinner('🌐 Fetching data from PubChem...'):
                    time.sleep(0.3)  # Rate limiting for PubChem
                    
                    compounds = pcp.get_compounds(compound_name, 'name')
                    
                    if compounds:
                        compound = compounds[0]
                        smiles = compound.isomeric_smiles
                        
                        # Display compound information
                        st.subheader("Compound Information")
                        st.write(f"**Molecular Formula:** {compound.molecular_formula}")
                        st.write(f"**Molecular Weight:** {compound.molecular_weight}")
                        st.write(f"**IUPAC Name:** {compound.iupac_name}")
                        st.write(f"**SMILES:** {smiles}")
                        
                        # Save to database
                        with st.spinner('💾 Adding to database...'):
                            success = insert_compound(
                                name=compound_name,
                                formula=compound.molecular_formula,
                                weight=compound.molecular_weight,
                                iupac=compound.iupac_name,
                                smiles=smiles
                            )
                            
                            if success:
                                st.success("✅ Compound added to database!")
                        
                        st.warning("📡 Data sourced from PubChem. There may be discrepancies with chemical formula.")
                        
                        if '.' in smiles or '+' in smiles or '-' in smiles:
                            st.write("3D structure not available for ionic compounds like salts.")
                        else:
                            with st.spinner('🧬 Generating 3D structure...'):
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
                        st.error("❌ Compound not found in the database or PubChem.")
            
            cursor.close()
            conn.close()
                
        except Exception as e:
            st.error(f"Error: {e}")
