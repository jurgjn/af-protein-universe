
import ast, random, os, tempfile, time, sqlite3, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, streamlit_ext as ste, st_aggrid, prody, py3Dmol, stmol, Bio

def RowSelectedDataFrame(df_, pre_selected_rows=[0]):
    gb = st_aggrid.GridOptionsBuilder.from_dataframe(df_)
    gb.configure_selection(selection_mode='single', pre_selected_rows=pre_selected_rows)
    gb.configure_grid_options(domLayout='normal')
    gridOptions = gb.build()
    gridResponse = st_aggrid.AgGrid(df_,
        gridOptions=gridOptions,
        columns_auto_size_mode=st_aggrid.ColumnsAutoSizeMode.FIT_ALL_COLUMNS_TO_VIEW,
        height=400,
        width='100%',
        enable_enterprise_modules=False,
    )
    return gridResponse

def RowSelectedDataFrameGet(gr_):
    if not(len(gr_['selected_rows']) > 0): time.sleep(5) # Prevent annoying row-not-selected errors during loading
    return gr_['selected_rows'][0]

def strip_af_cif(s):
    return s.removesuffix('-F1-model_v3.cif').removeprefix('AF-')

def uf(x):
    return '{:,}'.format(x)

@st.cache_resource
def read_af2_v3_(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v3.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

#st.cache_resource.clear()
st.set_page_config(layout='wide')

#entryID = st.selectbox(label='UniProtKB_ac', options=read_clusters_()['repID'])

st.write('# Browse clusters')
entryID = ste.text_input(label='UniProtKB_ac:', value='A0A2R8Y619', key='UniProtKB_ac')

fig3_sqlite = 'pages/Figure_3_Human-Centric_Clusters.sqlite'

def query_entryID(entryID):
    with sqlite3.connect(fig3_sqlite) as con:
        df_ = pd.read_sql(f'select * from members where entryID == "{entryID}"', con)
    return df_

def query_repID(repID):
    with sqlite3.connect(fig3_sqlite) as con:
        df_ = pd.read_sql(f'select * from members where repID == "{repID}"', con)
    return df_

def query_clusters_repID(repID):
    with sqlite3.connect(fig3_sqlite) as con:
        df_ = pd.read_sql(f'select * from clusters where repID == "{repID}"', con)
    return df_

q_ = query_entryID(entryID)
if len(q_) == 1:
    repID = q_['repId'].squeeze()
    st.write(f'{entryID} belongs to cluster {repID}')
else:
    st.write(f'# {entryID} not found..')

col1, col2, col3 = st.columns(3)

with col1:
    st.write(f'### Cluster {repID}')
    st.write('Cluster metadata:')
    st.dataframe(query_clusters_repID(repID).set_index(['repID']), use_container_width=True)
    df_cluster_ = query_repID(repID)
    entryID_index_ = df_cluster_.query('entryID == @entryID').index.values[0]
    #st.write(entryID_index_, len(df_cluster_))
    st.write('Cluster members:')
    gr_cluster_ = RowSelectedDataFrame(df_cluster_, pre_selected_rows=[int(entryID_index_)])

with col2:
    selEntryID = RowSelectedDataFrameGet(gr_cluster_)['entryID']
    st.write(f'### Selected member ({selEntryID})')
    pdb_mobile_ = read_af2_v3_(selEntryID)
    pdb_target_ = read_af2_v3_(repID)

    tmp_mobile_ = tempfile.NamedTemporaryFile(suffix='.pdb')
    tmp_target_ = tempfile.NamedTemporaryFile(suffix='.pdb')

    with open(tmp_mobile_.name, 'w') as fh:
        fh.write(pdb_mobile_)

    with open(tmp_target_.name, 'w') as fh:
        fh.write(pdb_target_)

    mol_mobile_ = prody.parsePDB(tmp_mobile_.name)
    mol_target_ = prody.parsePDB(tmp_target_.name)
    superimposed_mobile_, match_mobile_, match_target_, seqid, overlap = prody.matchAlign(mol_mobile_, mol_target_, seqid=20, overlap=20)

    tmp_super_ = tempfile.NamedTemporaryFile(suffix='.pdb')
    prody.writePDB(tmp_super_.name, superimposed_mobile_)

    with open(tmp_super_.name) as fh:
        pdb_superimposed_ = fh.read()

    # Add structure
    xyzview2 = py3Dmol.view()

    xyzview2.addModel(pdb_superimposed_, format='pdb')
    xyzview2.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})

    #xyzview2.addModel(pdb_mobile_, format='pdb')
    #xyzview2.setStyle({'model': 1}, {'cartoon': {'color': 'red'}})

    # Back matter
    xyzview2.setBackgroundColor('#eeeeee')
    xyzview2.zoomTo()
    stmol.showmol(xyzview2, height = 400, width=400)

with col3:
    st.write(f'### Cluster representative ({repID})')

    # Add structure
    xyzview3 = py3Dmol.view()
    xyzview3.addModel(pdb_target_, format='pdb')
    xyzview3.setStyle({'model': 0}, {'cartoon': {'color': 'grey'}})

    # Back matter
    xyzview3.setBackgroundColor('#eeeeee')
    xyzview3.zoomTo()
    stmol.showmol(xyzview3, height = 400, width=400)
