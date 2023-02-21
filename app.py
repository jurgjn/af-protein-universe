
import ast, random, os, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, st_aggrid, py3Dmol, stmol

st.set_page_config(layout='wide')

@st.cache_resource #https://docs.streamlit.io/library/advanced-features/caching
def read_pockets_():
    return pd.read_csv('data/pockets_score60_pLDDT90.tsv', sep='\t')

@st.cache_resource
def read_deepfri_():
    return pd.read_csv('data/pockets_score60_pLDDT90_DeepFRI_predictions.tsv', sep='\t')

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
df_pockets_ = read_pockets_().drop(['pocket_xmin', 'pocket_xmax', 'pocket_ymin', 'pocket_ymax', 'pocket_zmin', 'pocket_zmax'], axis=1)\
                             .merge(read_deepfri_summary_(), left_on='UniProtKB_ac', right_on='struct_id', how='left').sort_values(['DeepFri_max_score'], ascending=False)
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
    resid_ = grid_response['selected_rows'][0]['pocket_resid']
    cl_file_ = grid_response['selected_rows'][0]['pocket_cl_file']
else:
    af2_id_ = df_pockets_.head(1).struct_id.squeeze()
    resid_ = df_pockets_.head(1).pocket_resid.squeeze()
    cl_file_ = df_pockets_.head(1).pocket_cl_file.squeeze()

pocket_resid_ = ast.literal_eval(resid_)

#selected_indices = st.selectbox('Select structure:', df_pockets_['struct_id'])
# https://github.com/deepmind/alphafold/issues/92#issuecomment-1005495687
# https://github.com/Intron7/alpha_viewer
# https://alphafold.ebi.ac.uk/files/AF-A0A1V6PM83-F1-model_v3.pdb

st.write(f'## {af2_id_}')
st.markdown(f'{af2_id_} in [UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry) / [AlphaFill](https://alphafill.eu/model?id={af2_id_}) / [Ensembl Bacteria](https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q={af2_id_};site=ensemblunit)')

col1, col2 = st.columns(2)
with col1:
    st.write('### DeepFRI GO/EC terms')

    #st.dataframe(read_deepfri_().query('Protein == @af2_id_').sort_values('Score', ascending=False).reset_index(drop=True), height=600, use_container_width=True)
    df_ = read_deepfri_().query('Protein == @af2_id_').sort_values('Score', ascending=False).reset_index(drop=True)
    gb = st_aggrid.GridOptionsBuilder.from_dataframe(df_)
    gb.configure_selection('single')
    gb.configure_grid_options(domLayout='normal')
    gridOptions = gb.build()

    grid_response = st_aggrid.AgGrid(df_,
        gridOptions=gridOptions,
        #https://discuss.streamlit.io/t/is-there-a-way-to-autosize-all-columns-by-default-on-rendering-with-streamlit-aggrid/31841/2
        #fit_columns_on_grid_load=True,
        columns_auto_size_mode=st_aggrid.ColumnsAutoSizeMode.FIT_ALL_COLUMNS_TO_VIEW,
        height=400,
        width='100%',
        enable_enterprise_modules=False,
    )

    if len(grid_response['selected_rows']) > 0:
        go_ec_ = grid_response['selected_rows'][0]['GO_term/EC_number']
    else:
        go_ec_ = None

with col2:
    st.write('### Structure')
    #st.write('Blue = pocket residues')
    saliency_ = None
    if not(go_ec_) is None:
        fp_ = f'data/saliency/{af2_id_}_saliency.tsv'
        if os.path.getsize(fp_) > 1:
            df_ = pd.read_csv(fp_, sep='\t')
            if go_ec_ in df_.columns:
                st.write(f'Residues colored by saliency for {go_ec_}')#: {df_[go_ec_].head(3)}')
                saliency_ = df_[go_ec_]
            else:
                st.write('No saliency available (CC/BP?)')

    pdb_ = read_af2_v3_(af2_id_)

    #colors_pocket = {i: '#0072b2' for i in pocket_resid_}
    cmap_ = sns.color_palette("viridis", as_cmap=True)
    if not(saliency_) is None:
        st.write(cmap_)
        colors_pocket = {i + 1: matplotlib.colors.to_hex(cmap_(val_)) for i, val_ in saliency_.items() }
    else:
        colors_pocket = {}
    #st.write(colors_pocket)
    xyzview = py3Dmol.view()

    # Add structure
    xyzview.addModel(pdb_, format='pdb')
    xyzview.setStyle({'model': 0}, {
        'cartoon': {
            #'color':'spectrum'
            'colorscheme': {
                'prop': 'resi',
                'map': colors_pocket,
        }
    }})

    # Add pocket
    with open(os.path.join('data/pockets', os.path.basename(cl_file_))) as fh:
        pocket_ = fh.read()
    xyzview.addModel(pocket_, format='pdb')
    xyzview.setStyle({'model': -1}, {})
    xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'pink'}, {'model': -1})

    # Back matter
    xyzview.setBackgroundColor('#eeeeee')
    xyzview.zoomTo()
    stmol.showmol(xyzview, height = 800, width=800)
