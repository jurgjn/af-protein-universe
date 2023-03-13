
import ast, random, os, tempfile, time, sqlite3, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, streamlit_ext as ste, st_aggrid, prody, py3Dmol, stmol, Bio

st.set_page_config(
    page_title='Human-Centric Clusters',
    page_icon='🔬',
    layout='wide',
)
st.cache_resource.clear()

def strip_af_cif(s):
    return s.removesuffix('-F1-model_v3.cif').removeprefix('AF-')

def uf(x):
    return '{:,}'.format(x)

def RowSelectedDataFrame(df_, selected_row_index, height=400):
    gb = st_aggrid.GridOptionsBuilder.from_dataframe(df_)
    gb.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[ selected_row_index ])
    gb.configure_grid_options(domLayout='normal')
    #gb.configure_pagination()
    #https://github.com/PablocFonseca/streamlit-aggrid/issues/57
    gb.configure_grid_options(onFirstDataRendered=st_aggrid.JsCode("""
    function(e) { 
        e.api.ensureIndexVisible(%d, 'middle');
    }
    """ % (selected_row_index,)).js_code)
    gridOptions = gb.build()
    gridResponse = st_aggrid.AgGrid(df_,
        gridOptions=gridOptions,
        #update_mode=st_aggrid.GridUpdateMode.SELECTION_CHANGED,
        fit_columns_on_grid_load=True,
        height=height,
        width='100%',
        enable_enterprise_modules=False,
        allow_unsafe_jscode=True,
    )
    return gridResponse

def RowSelectedDataFrameGet(gr_):
    if not(len(gr_['selected_rows']) > 0): time.sleep(5) # Prevent annoying row-not-selected errors during loading
    return gr_['selected_rows'][0]

@st.cache_resource
def read_af2_v3_(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v3.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')
#---

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

def query_clusters():
    with sqlite3.connect(fig3_sqlite) as con:
        df_ = pd.read_sql(f'select * from clusters', con)
    return df_

st.write('# Human-Centric Clusters')

entryID = ste.text_input(label='Search cluster membership for (enter Uniprot ID):', value='A0A2R8Y619', key='entryID')
q_ = query_entryID(entryID)
if len(q_) == 1:
    repID = q_['repID'].squeeze()
    repID_index_ = int(query_clusters().query('repID == @repID').index.values[0])
    st.write(f'{entryID} belongs to cluster {repID}')
else:
    st.write(f'{entryID} not found')
    repID_index_ = 0

gr_clusters_ = RowSelectedDataFrame(query_clusters(), selected_row_index=repID_index_, height=200)
repID = RowSelectedDataFrameGet(gr_clusters_)['repID']
if entryID in set(query_repID(repID)['entryID']):
    st.experimental_set_query_params(repID=repID, entryID=entryID)
else:
    st.experimental_set_query_params(repID=repID, entryID=repID)
    entryID=repID

col1, col2 = st.columns(spec=[1,2])

with col1:
    st.write(f'### Cluster {repID} members:')
    df_cluster_members_ = query_repID(repID)[['entryID', 'taxID']]
    entryID_index_ = df_cluster_members_.query('entryID == @entryID').index.values[0]
    gr_cluster_members_ = RowSelectedDataFrame(df_cluster_members_, selected_row_index=int(entryID_index_))
    entryID = RowSelectedDataFrameGet(gr_cluster_members_)['entryID']
    st.experimental_set_query_params(repID=repID, entryID=entryID)

with col2:
    st.write(f'### Structure ({entryID})')

    entryID_pdb = read_af2_v3_(entryID)

    # Add structure
    xyzview3 = py3Dmol.view()
    xyzview3.addModel(entryID_pdb, format='pdb')
    xyzview3.setStyle({'model': 0}, {'cartoon': {'color': 'blue'}})

    # Back matter
    xyzview3.setBackgroundColor('#eeeeee')
    xyzview3.zoomTo()
    stmol.showmol(xyzview3, height = 400, width=800)
