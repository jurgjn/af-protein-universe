
import ast, random, os, urllib.request
import matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, st_aggrid, py3Dmol, stmol

st.set_page_config(layout="wide")

@st.cache_resource #https://docs.streamlit.io/library/advanced-features/caching
def read_pockets_():
    return pd.read_csv('pockets_score60_pLDDT90.tsv', sep='\t')

@st.cache_resource
def read_deepfri_():
    return pd.read_csv('pockets_score60_pLDDT90_DeepFRI_predictions.tsv', sep='\t')

@st.cache_resource
def read_af2_v3_(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v3.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

#@st.cache_resource
def read_deepfri_summary_():
    df_ = read_deepfri_().sort_values('Score', ascending=False).groupby('Protein').head(1).rename({
        'Protein': 'struct_id',
        'Score': 'DeepFri_max_score',
        'GO_term/EC_number name': 'DeepFri_max_GO/EC_name',
    }, axis=1)
    return df_[['struct_id', 'DeepFri_max_score', 'DeepFri_max_GO/EC_name']]

st.write('# Enzyme activity predictions for dark clusters')
st.write('Click on row to view structure + DeepFRI summary')
df_pockets_ = read_pockets_().drop(['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'cl_file', 'cl_isfile'], axis=1)\
                             .merge(read_deepfri_summary_(), left_on='struct_id', right_on='struct_id', how='left')
#st.dataframe(df_pockets_, height=200, use_container_width=True)

gb = st_aggrid.GridOptionsBuilder.from_dataframe(df_pockets_)
gb.configure_selection('single')
gb.configure_grid_options(domLayout='normal')
gridOptions = gb.build()

grid_response = st_aggrid.AgGrid(df_pockets_,
    gridOptions=gridOptions,
    #https://discuss.streamlit.io/t/is-there-a-way-to-autosize-all-columns-by-default-on-rendering-with-streamlit-aggrid/31841/2
    #fit_columns_on_grid_load=True,
    columns_auto_size_mode=st_aggrid.ColumnsAutoSizeMode.FIT_ALL_COLUMNS_TO_VIEW,
    height=200,
    width='100%',
    enable_enterprise_modules=False,
)
if len(grid_response['selected_rows']) > 0:
    af2_id_ = grid_response['selected_rows'][0]['struct_id']
    resid_ = grid_response['selected_rows'][0]['resid']
else:
    af2_id_ = df_pockets_.head(1).struct_id.squeeze()
    resid_ = df_pockets_.head(1).resid.squeeze()

pocket_resid_ = ast.literal_eval(resid_)

#selected_indices = st.selectbox('Select structure:', df_pockets_['struct_id'])
# https://github.com/deepmind/alphafold/issues/92#issuecomment-1005495687
# https://github.com/Intron7/alpha_viewer
# https://alphafold.ebi.ac.uk/files/AF-A0A1V6PM83-F1-model_v3.pdb

st.write(f'## {af2_id_}')
st.markdown(f'[{af2_id_} in UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry)')
st.markdown(f'[{af2_id_} in AlphaFill](https://alphafill.eu/model?id={af2_id_})')
st.markdown(f'[{af2_id_} in Ensembl Bacteria](https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q={af2_id_};site=ensemblunit)')

col1, col2 = st.columns(2)
with col1:
    st.write('### Structure')
    st.write('Blue = pocket residues')
    pdb_ = read_af2_v3_(af2_id_)
    colors_pocket = {i: '#0072b2' for i in pocket_resid_}

    xyzview = py3Dmol.view(data=pdb_, style={'stick':{}})
    xyzview.setStyle({'cartoon': {
        #'color':'spectrum'
        'colorscheme': {
            'prop': 'resi',
            'map': colors_pocket,
    }}})
    xyzview.setBackgroundColor('#D3D3D3')
    stmol.showmol(xyzview, height = 800, width=800)

with col2:
    st.write('### DeepFRI GO/EC terms')
    st.dataframe(read_deepfri_().query('Protein == @af2_id_').sort_values('Score', ascending=False).reset_index(drop=True), height=600, use_container_width=True)

    #fig, ax = plt.subplots()
    #sns.heatmap([[1,2,3], [2,3,2]], ax=ax)
    #st.write(fig)
