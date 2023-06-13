
import ast, random, os, tempfile, time, sqlite3, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, streamlit_ext as ste, st_aggrid, py3Dmol, stmol

st.set_page_config(
    page_title='Putative Novel Enzymes',
    page_icon='🔬',
    layout='wide',
)
st.cache_resource.clear()

#---
def strip_af_cif(s):
    return s.removesuffix('-F1-model_v3.cif').removeprefix('AF-')

def uf(x):
    return '{:,}'.format(x)

def select_dataframe_row(df_, selected_row_index, height=400):
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
    if not(len(gridResponse['selected_rows']) > 0): time.sleep(5) # Prevent annoying row-not-selected errors during loading
    if len(gridResponse['selected_rows']) == 0: return None
    return gridResponse['selected_rows'][0]

@st.cache_resource
def read_af2_v3_(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v3.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')
#---

@st.cache_resource #https://docs.streamlit.io/library/advanced-features/caching
def read_pockets_():
    return pd.read_csv('pages/Figure_2_Putative_Novel_Enzymes/pockets_score60_pLDDT90.tsv', sep='\t')

@st.cache_resource
def read_deepfri_():
    return pd.read_csv('pages/Figure_2_Putative_Novel_Enzymes/pockets_score60_pLDDT90_DeepFRI_predictions.tsv', sep='\t')

@st.cache_resource
def read_deepfri_summary_():
    df_ = read_deepfri_().sort_values('Score', ascending=False).groupby('Protein').head(1).rename({
        'Protein': 'struct_id',
        'Score': 'DeepFri_max_score',
        'GO_term/EC_number name': 'DeepFri_max_GO/EC_name',
    }, axis=1)
    return df_[['struct_id', 'DeepFri_max_score', 'DeepFri_max_GO/EC_name']]

@st.cache_resource
def read_pockets_with_deepfri_summary_():
    df_ = read_pockets_().merge(read_deepfri_summary_(), left_on='UniProtKB_ac', right_on='struct_id', how='left')\
        .sort_values(['DeepFri_max_score'], ascending=False).reset_index(drop=True)
    df_['struct_resid_in_pockets'] = df_['pocket_nresid'] / df_['struct_nresid']
    return df_

st.write('# Enzyme activity predictions for dark clusters')

tab1, tab2 = st.tabs(['Browse examples', 'Global statistics'])

with tab1:
    st.write('#### All structures/pockets')
    df_pockets_ = read_pockets_with_deepfri_summary_()
    if st.checkbox('Hide structures with a general lack of compactness (struct_resid_in_pockets > 0.4)', value=True):
        df_pockets_ = df_pockets_.query('struct_resid_in_pockets <= 0.4').reset_index(drop=True)

    cols_drop_ = ['pocket_xmin', 'pocket_xmax', 'pocket_ymin', 'pocket_ymax', 'pocket_zmin', 'pocket_zmax',
        'pocket_cl_file', 'pocket_cl_isfile', 'struct_id',
        'pocket_n_points', 'pocket_energy', 'pocket_energy_per_vol', 'pocket_rgyr', 'pocket_buriedness', 'pocket_resid']
    df_pockets_aggrid_ = df_pockets_.drop(cols_drop_, axis=1).round(
        {'pocket_score': 1, 'pocket_mean_pLDDT': 1, 'DeepFri_max_score': 2, 'struct_resid_in_pockets': 2})
    
    try:
        entryID = st.experimental_get_query_params().get('entryID')[0]
        entryID_index_ = int(df_pockets_aggrid_.query('UniProtKB_ac == @entryID').index.values[0])
    except:
        entryID_index_ = 0

    row_ = select_dataframe_row(df_pockets_aggrid_, selected_row_index=entryID_index_)
    st.write(f'{uf(len(df_pockets_aggrid_))} pockets shown')

    af2_id_ = row_['UniProtKB_ac']
    pocket_id_ = row_['pocket_id']

    resid_ = df_pockets_.query('UniProtKB_ac == @af2_id_ & pocket_id == @pocket_id_').squeeze().pocket_resid
    cl_file_ = df_pockets_.query('UniProtKB_ac == @af2_id_ & pocket_id == @pocket_id_').squeeze().pocket_cl_file
    pocket_resid_ = ast.literal_eval(resid_)

    #selected_indices = st.selectbox('Select structure:', df_pockets_['struct_id'])
    # https://github.com/deepmind/alphafold/issues/92#issuecomment-1005495687
    # https://github.com/Intron7/alpha_viewer
    # https://alphafold.ebi.ac.uk/files/AF-A0A1V6PM83-F1-model_v3.pdb

    col1, col2 = st.columns(2)
    with col1:
        st.write(f'#### DeepFRI GO/EC terms for {af2_id_}')
        tab11, tab12 = st.tabs(['Table', 'Heatmap'])
        with tab11:
            df_terms_ = read_deepfri_().query('Protein == @af2_id_').sort_values('Score', ascending=False).round({'Score': 2,}).reset_index(drop=True)
            if st.checkbox('Hide non-significant terms (Score < 0.5)', value=True):
                df_terms_ = df_terms_.query('Score >= 0.5')
            if st.checkbox('Hide terms without saliency data', value=True):
                df_terms_ = df_terms_.query('(DeepFRI_ont == "EC") | DeepFRI_ont == "MF"')

            rs_terms_ = select_dataframe_row(df_terms_, selected_row_index=0)
            st.write(f'{af2_id_} in [UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry) / [AlphaFill](https://alphafill.eu/model?id={af2_id_}) / [Ensembl Bacteria](https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q={af2_id_};site=ensemblunit)')
            try:
                go_ec_ = rs_terms_['GO_term/EC_number']
            except:
                go_ec_ = None
        
        with tab12:
            df_pockets_q_ = df_pockets_.query('UniProtKB_ac == @af2_id_')
            nresid_ = df_pockets_q_['struct_nresid'].drop_duplicates().squeeze()

            df_heatmap_ = pd.DataFrame(index=range(1, nresid_ + 1))
            for i, r in df_pockets_q_.iterrows():
                col_ = f'pocket_id_{r.pocket_id}'
                val_ = ast.literal_eval(r['pocket_resid'])
                df_heatmap_.loc[val_, col_] = 1

            fp_ = f'pages/Figure_2_Putative_Novel_Enzymes/saliency/{af2_id_}_saliency.tsv'
            #st.write(os.path.getsize(fp_))
            if os.path.getsize(fp_) > 1:
                df_ = pd.read_csv(fp_, sep='\t')
                df_.index = range(1, len(df_) + 1)
                df_.index.name='resseq'
                #df_.columns = [col_ + '\nsecond_row' for col_ in df_.columns]

            df_heatmap_ = pd.concat([df_heatmap_.fillna(0), df_], axis=1)
            #st.dataframe(df_heatmap_)

            fig, ax = plt.subplots(figsize=(4, 4))
            sns.heatmap(df_heatmap_, vmin=0, vmax=1, cmap='viridis')
            plt.gca().set_ylabel('resseq')
            #plt.gca().xaxis.tick_top()
            plt.xticks(rotation=70);
            st.pyplot(fig)

    with col2:
        st.write(f'#### Structure/visualisation for {af2_id_}')
        #st.write('Blue = pocket residues')
        saliency_ = None
        if not(go_ec_) is None:
            fp_ = f'pages/Figure_2_Putative_Novel_Enzymes/saliency/{af2_id_}_saliency.tsv'
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
        with open(os.path.join('pages/Figure_2_Putative_Novel_Enzymes/pockets', os.path.basename(cl_file_))) as fh:
            pocket_ = fh.read()
        xyzview.addModel(pocket_, format='pdb')
        xyzview.setStyle({'model': -1}, {})
        xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'pink'}, {'model': -1})

        # Back matter
        xyzview.setBackgroundColor('#eeeeee')
        xyzview.zoomTo()
        stmol.showmol(xyzview, height=600, width=600)
        st.write(cmap_)

with tab2:
    def uf(x):
        return '{:,}'.format(x)

    st.markdown(f'- {uf(len(read_pockets_()))} pockets in {uf(read_pockets_().UniProtKB_ac.nunique())} structures with pocket_score > 60')

    n_deepfri_ = len(read_deepfri_())
    n_deepfri_signif_ = len(read_deepfri_().query("Score > 0.5"))
    st.markdown(f'- {uf(n_deepfri_)} DeepFRI predictions with {uf(n_deepfri_signif_)} significant at 0.5 (as used in the original paper)')

    for ont_ in ['MF', 'EC', 'BP', 'CC']:
        st.markdown(f'#### Term counts for {ont_}')

        col_1_, col_2_ = st.columns(2)

        value_counts_ = read_deepfri_().query('Score > 0.5 & DeepFRI_ont == @ont_')[['GO_term/EC_number', 'GO_term/EC_number name']].value_counts()
        with col_1_:
            fig, ax = plt.subplots(figsize=(4, 4))
            value_counts_.head(10).plot(kind='barh')
            plt.title(f'Top 10 enriched terms ({ont_})')
            plt.gca().invert_yaxis()
            st.pyplot(fig)
            #plt.savefig(f'DeepFRI_top10_{ont_}.svg', format='svg', bbox_inches='tight') # SVG version for manuscript: https://github.com/streamlit/streamlit/issues/796

        with col_2_:
            st.dataframe(value_counts_)

