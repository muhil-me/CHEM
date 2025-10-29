from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
import pubchempy as pcp
import time
import psycopg2
from psycopg2.extras import RealDictCursor

st.title("Molecular Visualiser")
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
        RETURNING id
        """
        cursor.execute(query, (name, formula, weight, iupac, smiles))
        result = cursor.fetchone()
        conn.commit()
        cursor.close()
        conn.close()
        return result
    except Exception as e:
        st.error(f"Database error: {e}")
        return None

def get_initial_compounds(limit=10):
    """Get initial compounds to display"""
    try:
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT name, formula, molecular_weight, iupac_name, smiles 
            FROM compounds 
            ORDER BY name 
            LIMIT %s
        """, (limit,))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return results
    except Exception as e:
        st.error(f"Error fetching compounds: {e}")
        return []

def search_compounds_in_db(query):
    """Search compounds in database by name, formula, or IUPAC name"""
    try:
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        search_pattern = f"%{query}%"
        cursor.execute("""
            SELECT name, formula, molecular_weight, iupac_name, smiles 
            FROM compounds 
            WHERE LOWER(name) LIKE LOWER(%s) 
               OR LOWER(formula) LIKE LOWER(%s)
               OR LOWER(iupac_name) LIKE LOWER(%s)
            ORDER BY name
            LIMIT 20
        """, (search_pattern, search_pattern, search_pattern))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return results
    except Exception as e:
        st.error(f"Search error: {e}")
        return []

def search_pubchem(compound_name):
    """Search PubChem for compound data"""
    try:
        compounds = pcp.get_compounds(compound_name, 'name')
        if compounds:
            compound = compounds[0]
            return {
                'name': compound_name,
                'formula': compound.molecular_formula,
                'molecular_weight': compound.molecular_weight,
                'iupac_name': compound.iupac_name if hasattr(compound, 'iupac_name') else '',
                'smiles': compound.isomeric_smiles
            }
        return None
    except Exception as e:
        st.error(f"PubChem search error: {e}")
        return None

def visualize_molecule(smiles):
    """Visualize molecule using py3Dmol"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error("Invalid SMILES string")
            return
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mol_block = Chem.MolToMolBlock(mol)
        
        view = py3Dmol.view(width=800, height=600)
        view.addModel(mol_block, 'mol')
        view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
        view.setBackgroundColor('white')
        view.zoomTo()
        
        st.components.v1.html(view._make_html(), height=600, scrolling=False)
    except Exception as e:
        st.error(f"Visualization error: {e}")

# Initialize session state
if 'selected_compound' not in st.session_state:
    st.session_state.selected_compound = None
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""

# Search interface
st.subheader("🔍 Search Compounds")

col1, col2 = st.columns([3, 1])

with col1:
    search_input = st.text_input(
        "Search by compound name, formula, or IUPAC name",
        value=st.session_state.search_query,
        placeholder="e.g., aspirin, C9H8O4, caffeine...",
        key="search_input"
    )

with col2:
    search_pubchem_btn = st.button("🌐 Search PubChem", use_container_width=True)

# Display results
if search_input:
    st.session_state.search_query = search_input
    
    # Search in database first
    db_results = search_compounds_in_db(search_input)
    
    if db_results:
        st.success(f"Found {len(db_results)} compound(s) in database")
        
        # Create a selectbox with compound names
        compound_names = [f"{r['name']} ({r['formula']})" for r in db_results]
        selected = st.selectbox("Select a compound:", compound_names)
        
        if selected:
            idx = compound_names.index(selected)
            st.session_state.selected_compound = db_results[idx]
    else:
        st.warning(f"No compounds found in database for '{search_input}'")
        st.info("💡 Click 'Search PubChem' to find this compound online")

# Handle PubChem search
if search_pubchem_btn and search_input:
    with st.spinner(f"Searching PubChem for '{search_input}'..."):
        pubchem_data = search_pubchem(search_input)
        
        if pubchem_data:
            st.success(f"Found '{pubchem_data['name']}' on PubChem!")
            
            # Display compound info
            st.markdown("#### Compound Information")
            col1, col2 = st.columns(2)
            with col1:
                st.write(f"**Name:** {pubchem_data['name']}")
                st.write(f"**Formula:** {pubchem_data['formula']}")
            with col2:
                st.write(f"**Molecular Weight:** {pubchem_data['molecular_weight']:.2f} g/mol")
                st.write(f"**IUPAC Name:** {pubchem_data['iupac_name'][:50]}...")
            
            # Ask if user wants to add to database
            if st.button("➕ Add to Database"):
                result = insert_compound(
                    pubchem_data['name'],
                    pubchem_data['formula'],
                    pubchem_data['molecular_weight'],
                    pubchem_data['iupac_name'],
                    pubchem_data['smiles']
                )
                if result:
                    st.success("Compound added to database!")
                    st.session_state.selected_compound = pubchem_data
                else:
                    st.info("Compound already exists in database")
                    st.session_state.selected_compound = pubchem_data
        else:
            st.error(f"Could not find '{search_input}' on PubChem")

# Display selected compound
if st.session_state.selected_compound:
    compound = st.session_state.selected_compound
    
    st.markdown("---")
    st.subheader(f"📊 {compound['name']}")
    
    # Display compound details
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Formula", compound['formula'])
    with col2:
        st.metric("Molecular Weight", f"{compound['molecular_weight']:.2f} g/mol")
    with col3:
        if compound.get('iupac_name'):
            st.write("**IUPAC Name:**")
            st.caption(compound['iupac_name'])
    
    # Display SMILES
    with st.expander("View SMILES"):
        st.code(compound['smiles'])
    
    # Visualize molecule
    st.markdown("#### 3D Structure")
    visualize_molecule(compound['smiles'])

# Show initial compounds if no search
elif not search_input:
    st.info("💡 Start typing to search, or browse initial compounds below")
    initial_compounds = get_initial_compounds()
    
    if initial_compounds:
        st.markdown("#### Available Compounds")
        compound_names = [f"{c['name']} ({c['formula']})" for c in initial_compounds]
        selected = st.selectbox("Browse compounds:", ["Select a compound..."] + compound_names)
        
        if selected != "Select a compound...":
            idx = compound_names.index(selected)
            st.session_state.selected_compound = initial_compounds[idx]
            st.rerun()
