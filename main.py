from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
from streamlit_ketcher import st_ketcher
import pubchempy as pcp
import time
import psycopg2
from psycopg2.extras import RealDictCursor

# Page configuration
st.set_page_config(
    page_title="Molecular Visualizer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Load external CSS file
try:
    with open('style.css') as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
except FileNotFoundError:
    st.warning("CSS file not found. Using default styling.")

# Hide Streamlit branding
st.markdown("""
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    [data-testid="stToolbar"] {visibility: hidden;}
    .stDeployButton {visibility: hidden;}
    .stTitle {border-bottom: none !important;}
    </style>
    """, unsafe_allow_html=True)

# Database configuration
DATABASE_URL = st.secrets["database"]["url"]

# ==================== DATABASE FUNCTIONS ====================
def get_db_connection():
    """Establish connection to PostgreSQL database"""
    return psycopg2.connect(DATABASE_URL)

def insert_compound(name, formula, weight, iupac, smiles):
    """Insert compound data into database"""
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

# ==================== UTILITY FUNCTIONS ====================
def generate_3d_structure(smiles):
    """Generate 3D molecular structure from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        view = py3Dmol.view(width=800, height=500)
        view.addModel(Chem.MolToMolBlock(mol), 'mol')
        view.setStyle({'stick': {}})
        view.setBackgroundColor('white')
        view.zoomTo()
        
        return view._make_html()
    except Exception as e:
        st.error(f"Error generating 3D structure: {e}")
        return None

def display_compound_info(compound_data):
    """Display compound information in a formatted layout"""
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown(f"** Molecular Formula:** `{compound_data.get('formula', 'N/A')}`")
        st.markdown(f"** Molecular Weight:** `{compound_data.get('molecular_weight', 'N/A')} g/mol`")
    
    with col2:
        st.markdown(f"** IUPAC Name:** {compound_data.get('iupac_name', 'N/A')}")
        st.markdown(f"** SMILES:** `{compound_data.get('smiles', 'N/A')}`")

def is_ionic_compound(smiles):
    """Check if compound is ionic (contains salts)"""
    return '.' in smiles or '+' in smiles or '-' in smiles

# ==================== MAIN APPLICATION ====================
st.markdown("<h1 style='text-align: center; color: #1f77b4; font-size: 3rem;'>🧬 Molecular Visualizer</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #666; font-size: 1.2rem; margin-bottom: 2rem;'>Explore molecular structures in stunning 3D</p>", unsafe_allow_html=True)

# Mode selection
st.markdown("---")
option = st.selectbox(" Select Mode", ["Generate 3D Structure", "Draw Molecule"], index=0)
st.markdown("---")

# ==================== MODE 1: GENERATE 3D STRUCTURE ====================
if option == "Generate 3D Structure":
    st.subheader(" Search by Compound Name")
    
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        compound_name = st.text_input(
            "Enter compound name",
            placeholder="e.g., caffeine, aspirin, glucose",
            key="compound_name_input",
            label_visibility="collapsed"
        )
        compound_name = compound_name.strip()
        
        generate_btn = st.button(" Generate 3D Structure", key="generate_3d_btn", use_container_width=True)
    
    if generate_btn:
        if not compound_name:
            st.warning(" Please enter a compound name")
        else:
            try:
                conn = get_db_connection()
                cursor = conn.cursor(cursor_factory=RealDictCursor)
                
                # Search in database
                query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(%s)"
                cursor.execute(query, (compound_name,))
                result = cursor.fetchone()
                
                if result:
                    # Found in database
                    st.success(" Data retrieved from local database")
                    
                    st.markdown("###  Compound Information")
                    display_compound_info(result)
                    
                    smiles = result['smiles']
                    
                    if is_ionic_compound(smiles):
                        st.warning(" 3D structure not available for ionic compounds like salts.")
                    else:
                        st.markdown("###  3D Molecular Structure")
                        with st.spinner('🧬 Generating 3D structure...'):
                            viewer_html = generate_3d_structure(smiles)
                            if viewer_html:
                                st.components.v1.html(viewer_html, height=550)
                
                else:
                    # Not found in database, query PubChem
                    with st.spinner(' Fetching data from PubChem...'):
                        time.sleep(0.3)  # Rate limiting
                        
                        compounds = pcp.get_compounds(compound_name, 'name')
                        
                        if compounds:
                            compound = compounds[0]
                            smiles = compound.isomeric_smiles
                            
                            st.warning(" Data sourced from PubChem")
                            
                            st.markdown("###  Compound Information")
                            compound_data = {
                                'formula': compound.molecular_formula,
                                'molecular_weight': compound.molecular_weight,
                                'iupac_name': compound.iupac_name or 'N/A',
                                'smiles': smiles
                            }
                            display_compound_info(compound_data)
                            
                            # Save to database
                            with st.spinner(' Adding to database...'):
                                success = insert_compound(
                                    name=compound_name,
                                    formula=compound.molecular_formula,
                                    weight=compound.molecular_weight,
                                    iupac=compound.iupac_name,
                                    smiles=smiles
                                )
                                
                                if success:
                                    st.success(" Compound added to database!")
                            
                            if is_ionic_compound(smiles):
                                st.warning(" 3D structure not available for ionic compounds like salts.")
                            else:
                                st.markdown("###  3D Molecular Structure")
                                with st.spinner('🧬 Generating 3D structure...'):
                                    viewer_html = generate_3d_structure(smiles)
                                    if viewer_html:
                                        st.components.v1.html(viewer_html, height=550)
                        
                        else:
                            st.error(" Compound not found in the database or PubChem.")
                
                cursor.close()
                conn.close()
                    
            except Exception as e:
                st.error(f" Error: {e}")

# ==================== MODE 2: DRAW MOLECULE ====================
elif option == "Draw Molecule":
    st.subheader(" Draw Your Molecule")
    st.info(" Use the Ketcher editor below to draw your molecule structure")
    
    smiles_drawn = st_ketcher(key="ketcher_draw")
    
    if smiles_drawn:
        st.success(f" Generated SMILES: `{smiles_drawn}`")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("###  Compound Information")
            try:
                with st.spinner(' Identifying compound...'):
                    compounds = pcp.get_compounds(smiles_drawn, 'smiles')
                    
                    if compounds:
                        compound = compounds[0]
                        compound_name = compound.iupac_name or (compound.synonyms[0] if compound.synonyms else 'Unknown')
                        
                        compound_data = {
                            'formula': compound.molecular_formula,
                            'molecular_weight': compound.molecular_weight,
                            'iupac_name': compound_name,
                            'smiles': smiles_drawn
                        }
                        display_compound_info(compound_data)
                    else:
                        st.warning(" No compound found for this SMILES in PubChem database.")
                        st.markdown(f"** SMILES:** `{smiles_drawn}`")
                        
            except Exception as e:
                st.error(f" Error fetching compound information: {e}")
        
        with col2:
            st.markdown("###  3D Structure")
            if not is_ionic_compound(smiles_drawn):
                with st.spinner('🧬 Generating 3D structure...'):
                    viewer_html = generate_3d_structure(smiles_drawn)
                    if viewer_html:
                        st.components.v1.html(viewer_html, height=400)
            else:
                st.warning(" 3D structure not available for ionic compounds.")

# Footer
st.markdown("<br><br>", unsafe_allow_html=True)
st.markdown("<hr>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #999; font-size: 0.9rem;'>Powered by RDKit, PubChem & py3Dmol | Built with Streamlit</p>", unsafe_allow_html=True)
