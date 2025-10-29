from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
import pubchempy as pcp
import psycopg2
from psycopg2.extras import RealDictCursor

# Page config
st.set_page_config(page_title="Molecular Visualizer", page_icon="🧬", layout="wide")

# Custom CSS
st.markdown("""
    <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    [data-testid="stToolbar"] {visibility: hidden;}
    .stDeployButton {visibility: hidden;}
    
    /* Modern styling */
    .main {
        padding: 2rem;
    }
    
    .compound-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 15px;
        color: white;
        margin-bottom: 2rem;
        box-shadow: 0 10px 30px rgba(0,0,0,0.2);
    }
    
    .metric-card {
        background: white;
        padding: 1.5rem;
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        text-align: center;
        border-left: 4px solid #667eea;
    }
    
    .metric-label {
        font-size: 0.9rem;
        color: #666;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    
    .metric-value {
        font-size: 1.5rem;
        color: #333;
        font-weight: bold;
        margin-top: 0.5rem;
    }
    
    .search-header {
        text-align: center;
        margin-bottom: 2rem;
    }
    
    .search-header h1 {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 3rem;
        margin-bottom: 0.5rem;
    }
    
    .stButton>button {
        width: 100%;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        padding: 0.75rem;
        border-radius: 8px;
        font-weight: 600;
        transition: transform 0.2s;
    }
    
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
    }
    
    .info-box {
        background: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #667eea;
        margin: 1rem 0;
    }
    </style>
""", unsafe_allow_html=True)

# Header
st.markdown('<div class="search-header"><h1>🧬 Molecular Visualizer</h1><p style="color: #666;">Search and explore molecular structures in 3D</p></div>', unsafe_allow_html=True)

DATABASE_URL = st.secrets["database"]["url"]

# Database functions
def get_db_connection():
    return psycopg2.connect(DATABASE_URL)

def insert_compound(name, formula, weight, iupac, smiles):
    try:
        conn = get_db_connection()
        cursor = conn.cursor()
        query = """
        INSERT INTO compounds (name, formula, molecular_weight, iupac_name, smiles)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT (name) DO UPDATE 
        SET formula = EXCLUDED.formula,
            molecular_weight = EXCLUDED.molecular_weight,
            iupac_name = EXCLUDED.iupac_name,
            smiles = EXCLUDED.smiles
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

def search_compounds_in_db(query):
    """Search compounds with autocomplete"""
    if not query or len(query) < 2:
        return []
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
            ORDER BY 
                CASE 
                    WHEN LOWER(name) = LOWER(%s) THEN 1
                    WHEN LOWER(name) LIKE LOWER(%s) THEN 2
                    ELSE 3
                END,
                name
            LIMIT 50
        """, (search_pattern, search_pattern, search_pattern, query, query + '%'))
        results = cursor.fetchall()
        cursor.close()
        conn.close()
        return results
    except Exception as e:
        st.error(f"Search error: {e}")
        return []

def get_compound_by_name(name):
    """Get specific compound by exact name"""
    try:
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute("""
            SELECT name, formula, molecular_weight, iupac_name, smiles 
            FROM compounds 
            WHERE name = %s
        """, (name,))
        result = cursor.fetchone()
        cursor.close()
        conn.close()
        return result
    except Exception as e:
        st.error(f"Error: {e}")
        return None

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
        st.error(f"PubChem error: {e}")
        return None

def visualize_molecule(smiles, compound_name="Molecule"):
    """Visualize molecule using py3Dmol"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.error("❌ Invalid SMILES string")
            return
        
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            st.warning("⚠️ Could not generate 3D coordinates, showing 2D structure")
            AllChem.Compute2DCoords(mol)
        else:
            AllChem.MMFFOptimizeMolecule(mol)
        
        mol_block = Chem.MolToMolBlock(mol)
        
        view = py3Dmol.view(width=800, height=600)
        view.addModel(mol_block, 'mol')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'radius': 0.4}})
        view.setBackgroundColor('#f8f9fa')
        view.zoomTo()
        view.spin(True)
        
        st.components.v1.html(view._make_html(), height=620, scrolling=False)
    except Exception as e:
        st.error(f"❌ Visualization error: {e}")

# Initialize session state
if 'selected_compound' not in st.session_state:
    st.session_state.selected_compound = None
if 'search_results' not in st.session_state:
    st.session_state.search_results = []
if 'show_pubchem' not in st.session_state:
    st.session_state.show_pubchem = False

# Search interface
col1, col2 = st.columns([4, 1])

with col1:
    search_query = st.text_input(
        "🔍 Search compounds",
        placeholder="Type compound name, formula, or IUPAC name (e.g., aspirin, C9H8O4, caffeine)...",
        key="search_box",
        label_visibility="collapsed"
    )

with col2:
    search_pubchem_btn = st.button("🌐 PubChem", use_container_width=True)

# Auto-search as user types
if search_query and len(search_query) >= 2:
    db_results = search_compounds_in_db(search_query)
    st.session_state.search_results = db_results
    
    if db_results:
        # Show dropdown suggestions
        st.markdown("#### 📋 Search Results")
        
        # Create selection list with better formatting
        options = [""] + [f"{r['name']} | {r['formula']}" for r in db_results]
        
        selected = st.selectbox(
            "Select a compound from results:",
            options,
            key="compound_selector",
            label_visibility="collapsed"
        )
        
        if selected and selected != "":
            # Extract compound name from selection
            compound_name = selected.split(" | ")[0]
            compound_data = get_compound_by_name(compound_name)
            if compound_data:
                st.session_state.selected_compound = compound_data
                st.session_state.show_pubchem = False
    else:
        st.info("💡 No results in database. Try searching PubChem using the button above!")
        st.session_state.show_pubchem = True

# Handle PubChem search
if search_pubchem_btn and search_query:
    with st.spinner(f"🔍 Searching PubChem for '{search_query}'..."):
        pubchem_data = search_pubchem(search_query)
        
        if pubchem_data:
            st.success(f"✅ Found '{pubchem_data['name']}' on PubChem!")
            
            col1, col2 = st.columns(2)
            with col1:
                st.markdown(f"**Name:** {pubchem_data['name']}")
                st.markdown(f"**Formula:** {pubchem_data['formula']}")
            with col2:
                st.markdown(f"**Molecular Weight:** {pubchem_data['molecular_weight']:.2f} g/mol")
                iupac_display = pubchem_data['iupac_name'][:80] + "..." if len(pubchem_data['iupac_name']) > 80 else pubchem_data['iupac_name']
                st.markdown(f"**IUPAC:** {iupac_display}")
            
            if st.button("➕ Add to Database & View", use_container_width=True):
                result = insert_compound(
                    pubchem_data['name'],
                    pubchem_data['formula'],
                    pubchem_data['molecular_weight'],
                    pubchem_data['iupac_name'],
                    pubchem_data['smiles']
                )
                if result:
                    st.success("✅ Compound saved to database!")
                st.session_state.selected_compound = pubchem_data
                st.rerun()
        else:
            st.error(f"❌ Could not find '{search_query}' on PubChem")

# Display selected compound
if st.session_state.selected_compound:
    compound = st.session_state.selected_compound
    
    st.markdown("---")
    
    # Compound header card
    st.markdown(f"""
    <div class="compound-card">
        <h2 style="margin: 0; font-size: 2rem;">🧪 {compound['name']}</h2>
        <p style="margin: 0.5rem 0 0 0; opacity: 0.9;">Molecular Structure Visualization</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Metrics in cards
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Chemical Formula</div>
            <div class="metric-value">{compound['formula']}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Molecular Weight</div>
            <div class="metric-value">{compound['molecular_weight']:.2f} g/mol</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        atoms = compound['formula']
        atom_count = sum(int(c) if c.isdigit() else 1 for c in compound['formula'] if c.isalnum())
        st.markdown(f"""
        <div class="metric-card">
            <div class="metric-label">Composition</div>
            <div class="metric-value">{atoms}</div>
        </div>
        """, unsafe_allow_html=True)
    
    # IUPAC name section
    if compound.get('iupac_name'):
        with st.expander("📖 IUPAC Name"):
            st.markdown(f"`{compound['iupac_name']}`")
    
    # SMILES section
    with st.expander("🧬 SMILES Notation"):
        st.code(compound['smiles'], language="text")
    
    # 3D Visualization
    st.markdown("### 🌐 3D Molecular Structure")
    st.markdown('<div class="info-box">💡 The molecule will rotate automatically. Drag to rotate manually, scroll to zoom.</div>', unsafe_allow_html=True)
    
    visualize_molecule(compound['smiles'], compound['name'])
    
    # Additional info
    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### 📊 Properties")
        st.markdown(f"- **Name:** {compound['name']}")
        st.markdown(f"- **Formula:** {compound['formula']}")
        st.markdown(f"- **Molecular Weight:** {compound['molecular_weight']:.2f} g/mol")
    
    with col2:
        if st.button("🔄 Clear Selection", use_container_width=True):
            st.session_state.selected_compound = None
            st.rerun()

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 2rem 0;">
    <p>Powered by RDKit, PubChem, and py3Dmol | Built with Streamlit</p>
</div>
""", unsafe_allow_html=True)
