import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st

# Page configuration
st.set_page_config(page_title="Molecular Visualizer", layout="wide", page_icon="🧬")

# Custom CSS Styling
st.markdown("""
    <style>
    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
    [data-testid="stToolbar"] {visibility: hidden;}
    .stDeployButton {visibility: hidden;}
    
    /* Main app background */
    .stApp {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    }
    
    /* Main content area */
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        background: rgba(255, 255, 255, 0.95);
        border-radius: 20px;
        box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
    }
    
    /* Title styling */
    h1 {
        color: #2c3e50;
        text-align: center;
        font-size: 3rem !important;
        font-weight: 800 !important;
        margin-bottom: 1rem !important;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    /* Subheader styling */
    h2, h3 {
        color: #34495e;
        font-weight: 700 !important;
    }
    
    /* Input field styling */
    .stTextInput > div > div > input {
        border-radius: 15px;
        border: 3px solid #667eea;
        padding: 15px;
        font-size: 18px;
        transition: all 0.3s ease;
        background: #f8f9fa;
    }
    
    .stTextInput > div > div > input:focus {
        border-color: #764ba2;
        box-shadow: 0 0 20px rgba(102, 126, 234, 0.3);
        background: white;
    }
    
    /* Button styling */
    .stButton > button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 25px;
        padding: 15px 40px;
        font-size: 18px;
        font-weight: bold;
        transition: all 0.3s ease;
        width: 100%;
        box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
    }
    
    .stButton > button:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.6);
    }
    
    /* Info box styling */
    .info-container {
        background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
        border-left: 5px solid #2196f3;
        padding: 20px;
        border-radius: 10px;
        margin: 20px 0;
        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
    }
    
    /* Legend box styling */
    .legend-box {
        background: white;
        border-radius: 15px;
        padding: 25px;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.2);
        border: 3px solid #667eea;
    }
    
    .legend-title {
        color: #2c3e50;
        font-size: 24px;
        font-weight: bold;
        margin-bottom: 20px;
        text-align: center;
        border-bottom: 3px solid #667eea;
        padding-bottom: 10px;
    }
    
    .legend-item {
        display: flex;
        align-items: center;
        margin: 15px 0;
        padding: 10px;
        border-radius: 8px;
        background: #f8f9fa;
        transition: all 0.3s ease;
    }
    
    .legend-item:hover {
        background: #e9ecef;
        transform: translateX(5px);
    }
    
    .color-circle {
        width: 30px;
        height: 30px;
        border-radius: 50%;
        margin-right: 15px;
        border: 2px solid #333;
        box-shadow: 0 2px 5px rgba(0, 0, 0, 0.2);
    }
    
    .element-info {
        flex: 1;
    }
    
    .element-name {
        font-weight: bold;
        font-size: 16px;
        color: #2c3e50;
    }
    
    .element-desc {
        font-size: 14px;
        color: #7f8c8d;
    }
    
    /* Error message styling */
    .stAlert {
        border-radius: 10px;
    }
    
    /* Molecule viewer container */
    .molecule-viewer {
        background: white;
        border-radius: 15px;
        padding: 20px;
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.1);
        border: 3px solid #667eea;
    }
    </style>
    """, unsafe_allow_html=True)

# Header with emoji
st.markdown("# 🧬 Molecular Visualizer")

# Create two columns - main content and sidebar legend
col1, col2 = st.columns([3, 1])

with col2:
    # Color Legend Sidebar
    st.markdown("""
        <div class="legend-box">
            <div class="legend-title">🎨 Atom Colors</div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #909090;"></div>
                <div class="element-info">
                    <div class="element-name">Carbon (C)</div>
                    <div class="element-desc">Grey/Dark Grey</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #FFFFFF; border: 2px solid #000;"></div>
                <div class="element-info">
                    <div class="element-name">Hydrogen (H)</div>
                    <div class="element-desc">White</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #FF0000;"></div>
                <div class="element-info">
                    <div class="element-name">Oxygen (O)</div>
                    <div class="element-desc">Red</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #0000FF;"></div>
                <div class="element-info">
                    <div class="element-name">Nitrogen (N)</div>
                    <div class="element-desc">Blue</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #FFFF00; border: 2px solid #999;"></div>
                <div class="element-info">
                    <div class="element-name">Sulfur (S)</div>
                    <div class="element-desc">Yellow</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #FFA500;"></div>
                <div class="element-info">
                    <div class="element-name">Phosphorus (P)</div>
                    <div class="element-desc">Orange</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #00FF00;"></div>
                <div class="element-info">
                    <div class="element-name">Chlorine (Cl)</div>
                    <div class="element-desc">Green</div>
                </div>
            </div>
            
            <div class="legend-item">
                <div class="color-circle" style="background: #8B4513;"></div>
                <div class="element-info">
                    <div class="element-name">Bromine (Br)</div>
                    <div class="element-desc">Brown</div>
                </div>
            </div>
        </div>
        
        <br>
        
        <div class="legend-box" style="margin-top: 20px;">
            <div class="legend-title">💡 How to Use</div>
            <div style="padding: 10px; color: #555; font-size: 14px; line-height: 1.6;">
                <p><strong>🖱️ Rotate:</strong> Click and drag</p>
                <p><strong>🔍 Zoom:</strong> Scroll wheel</p>
                <p><strong>↔️ Pan:</strong> Right-click and drag</p>
            </div>
        </div>
    """, unsafe_allow_html=True)

with col1:
    # Database connection
    conn = sqlite3.connect("data.db")  
    conn.row_factory = sqlite3.Row     
    cursor = conn.cursor()
    
    # Input field
    compound_name = st.text_input("🔬 Enter compound name:", placeholder="e.g., aspirin, glucose, caffeine")
    compound_name = compound_name.rstrip()
    
    # Generate button
    if st.button("🚀 Generate 3D Structure"):
        try:
            query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(?)"
            cursor.execute(query, (compound_name,))
            result = cursor.fetchone()
            
            if result:
                smiles = result['smiles']
                
                # Display compound information in styled box
                st.markdown('<div class="info-container">', unsafe_allow_html=True)
                st.subheader("📊 Compound Information")
                st.write(f"**🧪 Molecular Formula:** {result['formula']}")
                st.write(f"**⚖️ Molecular Weight:** {result['molecular_weight']} g/mol")
                st.write(f"**📝 IUPAC Name:** {result['iupac_name']}")
                st.write(f"**🔗 SMILES:** {result['smiles']}")
                st.markdown('</div>', unsafe_allow_html=True)
                
                # Check for ionic compounds
                if '.' in smiles or '+' in smiles or '-' in smiles:
                    st.warning("⚠️ 3D structure not available for ionic compounds like salts.")
                else:
                    # Generate 3D structure
                    with st.spinner("🔄 Generating 3D structure..."):
                        mol = Chem.MolFromSmiles(smiles)
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol)
                        AllChem.MMFFOptimizeMolecule(mol)
                        
                        # Create styled 3D viewer
                        st.markdown('<div class="molecule-viewer">', unsafe_allow_html=True)
                        st.subheader("🎯 3D Molecular Structure")
                        
                        view = py3Dmol.view(width=800, height=500)
                        view.addModel(Chem.MolToMolBlock(mol), 'mol')
                        view.setStyle({'stick': {'radius': 0.2}, 'sphere': {'scale': 0.3}})
                        view.setBackgroundColor('#f8f9fa')
                        view.zoomTo()
                        
                        viewer_html = view._make_html()
                        st.components.v1.html(viewer_html, height=550)
                        st.markdown('</div>', unsafe_allow_html=True)
                        
                        st.success("✅ 3D structure generated successfully!")
            
            else:
                st.error("❌ Compound not found in your local database")
                smiles = None
                
        except Exception as e:
            st.error(f"⚠️ Error: {e}")
            
        conn.close()
    
    # Instructions when no compound is searched yet
    else:
        st.info("👆 Enter a compound name above and click 'Generate 3D Structure' to visualize the molecule!")
