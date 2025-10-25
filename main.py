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

def render_card(title: str, inner_html: str) -> str:
        """Return HTML for a white rounded card with a title and arbitrary inner HTML/markdown."""
        return f"""
        <div class="card">
            <h2>{title}</h2>
            {inner_html}
        </div>
        """

# --- session state defaults -------------------------------------------------
if 'recent_searches' not in st.session_state:
    st.session_state['recent_searches'] = []
if 'compound_name' not in st.session_state:
    st.session_state['compound_name'] = ''


conn = sqlite3.connect("data.db")
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

st.title("Molecular visualiser")
st.write("Explore molecules from your local database. Select a compound and generate a 3D view.")

left, right = st.columns([1, 2])

with left:
    with st.container():
        st.markdown('<div class="input-box">', unsafe_allow_html=True)

        # helper: fetch matching names from DB using prefix match
        def get_matches(prefix: str, limit: int = 5):
            if not prefix:
                return []
            try:
                q = "SELECT name FROM compounds WHERE LOWER(name) LIKE LOWER(?) ORDER BY name LIMIT ?"
                cursor.execute(q, (f"%{prefix}%", limit))
                rows = cursor.fetchall()
                return [r['name'] for r in rows]
            except Exception:
                return []

        # Search box (text input). As the user types we populate a small "Search results" dropdown.
        compound_input = st.text_input("Search compound", value=st.session_state.get('compound_name', ''), key='compound_input')
        compound_input = compound_input.strip()

        # dynamic results dropdown (acts like the Google-like suggestions dropdown)
        suggestions = get_matches(compound_input, limit=5)
        selected_from_results = None
        if suggestions:
            opts = ["-- choose from results --"] + suggestions
            sel = st.selectbox("Search results", opts, index=0, key='search_results')
            if sel and sel != "-- choose from results --":
                selected_from_results = sel

        # Generate button below
        generate = st.button("Generate 3D Structure")
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown("<div class='small'>Tip: Try searching for compounds like aspirin or ethanol from the database.</div>", unsafe_allow_html=True)

        # Recent searches (clickable)
        recent = st.session_state.get('recent_searches', [])
        if recent:
            st.markdown("<div style='margin-top:10px'><strong>Recent searches:</strong></div>", unsafe_allow_html=True)
            rcols = st.columns(min(5, len(recent)))
            for i, r in enumerate(recent):
                with rcols[i % len(rcols)]:
                    if st.button(r, key=f'recent_{r}'):
                        st.session_state['compound_name'] = r
                        st.session_state['generate'] = True
                        st.experimental_rerun()

with right:
    view_placeholder = st.empty()


def generate_for_name(name: str):
    """Query DB and render info + 3D view for the given compound name."""
    if not name:
        st.warning("Please enter a compound name.")
        return

    query = "SELECT * FROM compounds WHERE LOWER(name) = LOWER(?)"
    cursor.execute(query, (name,))
    result = cursor.fetchone()

    if not result:
        st.error("Compound not found in your local database")
        return

    # update recent searches (keep unique, most recent first, limit 8)
    recent = st.session_state.get('recent_searches', [])
    if name in recent:
        recent.remove(name)
    recent.insert(0, name)
    st.session_state['recent_searches'] = recent[:8]

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
        return

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    view = py3Dmol.view(width=600, height=480)
    view.addModel(Chem.MolToMolBlock(mol), 'mol')
    view.setStyle({'stick': {}})
    view.zoomTo()

    viewer_html = view._make_html()
    viewer_html = viewer_html.replace(
        'https://cdn.jsdelivr.net/npm/3dmol@2.5.3/build/3Dmol-min.js',
        'https://3dmol.csb.pitt.edu/build/3Dmol-min.js',
    )

    html_wrapper = f"<div class='card mol-frame'>{viewer_html}</div>"
    components.html(html_wrapper, height=520)


# If the user picked a result from the dynamic dropdown, generate immediately
if 'selected_from_results' in locals() and selected_from_results:
    st.session_state['compound_name'] = selected_from_results
    generate_for_name(selected_from_results)

# If user clicked Generate, use the input value
if generate:
    st.session_state['compound_name'] = compound_input
    generate_for_name(compound_input)

conn.close()
