import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components


st.set_page_config(page_title="Molecular Visualiser", layout="wide")

st.markdown(
    """
    <style>
    :root{--bg:#0b0b0b;--muted:#bfbfbf;--card:#ffffff;--card-text:#000000}
    html, body, .stApp {background: var(--bg); color: var(--muted);}
    #MainMenu {visibility: hidden;} footer {visibility: hidden;} header {visibility: hidden;}
    .stButton>button {background:#000; color:#fff; border-radius:8px;}
    .card {background: var(--card); color: var(--card-text); border-radius:12px; padding:18px; box-shadow: 0 6px 18px rgba(0,0,0,0.6);}
    .card h2{color:var(--card-text); margin-top:0}
    .input-box {background:#111; padding:12px; border-radius:10px;}
    .small {color: var(--muted); font-size:0.9rem}
    /* ensure py3Dmol iframe fits rounded box */
    .mol-frame iframe{border-radius:10px; overflow:hidden}
    </style>
    """,
    unsafe_allow_html=True,
)


def render_card(title: str, inner_html: str):
    """Render a white rounded card with a title and arbitrary inner HTML/markdown."""
    card = f"""
    <div class="card">
      <h2>{title}</h2>
      {inner_html}
    </div>
    """
    st.markdown(card, unsafe_allow_html=True)


conn = sqlite3.connect("data.db")
conn.row_factory = sqlite3.Row
cursor = conn.cursor()


st.title("Molecular visualiser")
st.write("Explore molecules from your local database. Enter a compound name and generate a 3D view.")

left, right = st.columns([1, 2])

with left:
    with st.container():
        st.markdown('<div class="input-box">', unsafe_allow_html=True)
        compound_name = st.text_input("Compound name", value="")
        compound_name = compound_name.rstrip()
        generate = st.button("Generate 3D Structure")
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown("<div class='small'>Tip: try names like aspirin, ethanol, or search your DB.</div>", unsafe_allow_html=True)

with right:
    view_placeholder = st.empty()


if generate:
    try:
        query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(?)"
        cursor.execute(query, (compound_name,))
        result = cursor.fetchone()

        if result:
            smiles = result['smiles']
            info_html = (
                f"<p><strong>Molecular Formula:</strong> {result['formula']}</p>"
                f"<p><strong>Molecular Weight:</strong> {result['molecular_weight']}</p>"
                f"<p><strong>IUPAC Name:</strong> {result['iupac_name']}</p>"
                f"<p><strong>SMILES:</strong> {result['smiles']}</p>"
            )
            # show info in a rounded card
            with left:
                st.markdown(render_card("Compound information", info_html), unsafe_allow_html=True)

            if '.' in smiles or '+' in smiles or '-' in smiles:
                with left:
                    st.warning("3D structure not available for ionic compounds like salts.")
            else:
                mol = Chem.MolFromSmiles(smiles)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.MMFFOptimizeMolecule(mol)

                view = py3Dmol.view(width=600, height=480)
                view.addModel(Chem.MolToMolBlock(mol), 'mol')
                view.setStyle({'stick': {}})
                view.zoomTo()

                viewer_html = view._make_html()

                # Some environments block the default jsDelivr CDN. Prefer the official 3dmol host
                # and render the returned HTML using Streamlit components so script tags execute.
                viewer_html = viewer_html.replace(
                    'https://cdn.jsdelivr.net/npm/3dmol@2.5.3/build/3Dmol-min.js',
                    'https://3dmol.csb.pitt.edu/build/3Dmol-min.js',
                )

                # place viewer inside a rounded white card and use components.html so JS runs
                html_wrapper = f"<div class='card mol-frame'>{viewer_html}</div>"
                components.html(html_wrapper, height=520)

        else:
            st.error("Compound not found in your local database")

    except Exception as e:
        st.error(f"Error: {e}")

conn.close()

