from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
import pubchempy as pcp
import time
import psycopg2
from psycopg2.extras import RealDictCursor

st.title("Molecular visualiser")

# Load external CSS file
with open('style.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

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

compound_name = st.text_input("Enter compound name: ")
compound_name = compound_name.rstrip()

if st.button("Generate 3D Structure"):
    try:
        # Connect to Neon database
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # Search in database (changed ? to %s for PostgreSQL)
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
                
                # Display in two-column layout
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown('<div class="structure-column">', unsafe_allow_html=True)
                    st.markdown('<h3>3D Molecular Structure</h3>', unsafe_allow_html=True)
                    st.components.v1.html(viewer_html, height=450)
                    st.markdown('</div>', unsafe_allow_html=True)
                
                with col2:
                    st.markdown('<div class="legend-column">', unsafe_allow_html=True)
                    st.markdown('<h3>Bond & Atom Legend</h3>', unsafe_allow_html=True)
                    
                    # Bond Colors Legend
                    st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                    st.markdown('<div class="bond-color-swatch" style="background-color: #FFD700;"></div>', unsafe_allow_html=True)
                    st.markdown('''
                    <div class="bond-description">
                        <div class="bond-name">Single Bond (C-C)</div>
                        <div class="bond-type">Single covalent bond</div>
                        <div class="bond-details">One shared electron pair</div>
                    </div>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                    st.markdown('<div class="bond-color-swatch" style="background-color: #FF6B6B;"></div>', unsafe_allow_html=True)
                    st.markdown('''
                    <div class="bond-description">
                        <div class="bond-name">Double Bond (C=C)</div>
                        <div class="bond-type">Double covalent bond</div>
                        <div class="bond-details">Two shared electron pairs</div>
                    </div>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                    st.markdown('<div class="bond-color-swatch" style="background-color: #4ECDC4;"></div>', unsafe_allow_html=True)
                    st.markdown('''
                    <div class="bond-description">
                        <div class="bond-name">Triple Bond (C≡C)</div>
                        <div class="bond-type">Triple covalent bond</div>
                        <div class="bond-details">Three shared electron pairs</div>
                    </div>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                    st.markdown('<div class="bond-color-swatch" style="background-color: #95E1D3;"></div>', unsafe_allow_html=True)
                    st.markdown('''
                    <div class="bond-description">
                        <div class="bond-name">C-H Bond</div>
                        <div class="bond-type">Hydrogen bond</div>
                        <div class="bond-details">Carbon-hydrogen single bond</div>
                    </div>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                    st.markdown('<div class="bond-color-swatch" style="background-color: #F38181;"></div>', unsafe_allow_html=True)
                    st.markdown('''
                    <div class="bond-description">
                        <div class="bond-name">N/O Bonds</div>
                        <div class="bond-type">Heteroatom bonds</div>
                        <div class="bond-details">Nitrogen or oxygen single bonds</div>
                    </div>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    
                    st.markdown('<div class="legend-summary">', unsafe_allow_html=True)
                    st.markdown('''
                    <h4>💡 Understanding the Structure</h4>
                    <p><strong>Stick representation:</strong> Shows atoms as intersections and bonds as lines</p>
                    <p><strong>Bond strength:</strong> Indicated by line multiplicity (single, double, triple)</p>
                    <p><strong>Interactions:</strong> Colors help identify different bond types and atomic connections</p>
                    ''', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
                    st.markdown('</div>', unsafe_allow_html=True)
        
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
                            
                            # Display in two-column layout
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.markdown('<div class="structure-column">', unsafe_allow_html=True)
                                st.markdown('<h3>3D Molecular Structure</h3>', unsafe_allow_html=True)
                                st.components.v1.html(viewer_html, height=450)
                                st.markdown('</div>', unsafe_allow_html=True)
                            
                            with col2:
                                st.markdown('<div class="legend-column">', unsafe_allow_html=True)
                                st.markdown('<h3>Bond & Atom Legend</h3>', unsafe_allow_html=True)
                                
                                # Bond Colors Legend
                                st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                                st.markdown('<div class="bond-color-swatch" style="background-color: #FFD700;"></div>', unsafe_allow_html=True)
                                st.markdown('''
                                <div class="bond-description">
                                    <div class="bond-name">Single Bond (C-C)</div>
                                    <div class="bond-type">Single covalent bond</div>
                                    <div class="bond-details">One shared electron pair</div>
                                </div>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                                st.markdown('<div class="bond-color-swatch" style="background-color: #FF6B6B;"></div>', unsafe_allow_html=True)
                                st.markdown('''
                                <div class="bond-description">
                                    <div class="bond-name">Double Bond (C=C)</div>
                                    <div class="bond-type">Double covalent bond</div>
                                    <div class="bond-details">Two shared electron pairs</div>
                                </div>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                                st.markdown('<div class="bond-color-swatch" style="background-color: #4ECDC4;"></div>', unsafe_allow_html=True)
                                st.markdown('''
                                <div class="bond-description">
                                    <div class="bond-name">Triple Bond (C≡C)</div>
                                    <div class="bond-type">Triple covalent bond</div>
                                    <div class="bond-details">Three shared electron pairs</div>
                                </div>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                                st.markdown('<div class="bond-color-swatch" style="background-color: #95E1D3;"></div>', unsafe_allow_html=True)
                                st.markdown('''
                                <div class="bond-description">
                                    <div class="bond-name">C-H Bond</div>
                                    <div class="bond-type">Hydrogen bond</div>
                                    <div class="bond-details">Carbon-hydrogen single bond</div>
                                </div>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="legend-item">', unsafe_allow_html=True)
                                st.markdown('<div class="bond-color-swatch" style="background-color: #F38181;"></div>', unsafe_allow_html=True)
                                st.markdown('''
                                <div class="bond-description">
                                    <div class="bond-name">N/O Bonds</div>
                                    <div class="bond-type">Heteroatom bonds</div>
                                    <div class="bond-details">Nitrogen or oxygen single bonds</div>
                                </div>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                
                                st.markdown('<div class="legend-summary">', unsafe_allow_html=True)
                                st.markdown('''
                                <h4>💡 Understanding the Structure</h4>
                                <p><strong>Stick representation:</strong> Shows atoms as intersections and bonds as lines</p>
                                <p><strong>Bond strength:</strong> Indicated by line multiplicity (single, double, triple)</p>
                                <p><strong>Interactions:</strong> Colors help identify different bond types and atomic connections</p>
                                ''', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                                st.markdown('</div>', unsafe_allow_html=True)
                
                else:
                    st.error("❌ Compound not found in the database or PubChem.")
        
        cursor.close()
        conn.close()
            
    except Exception as e:
        st.error(f"Error: {e}")
