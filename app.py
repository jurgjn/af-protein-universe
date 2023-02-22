
import ast, random, os, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, st_aggrid, py3Dmol, stmol

st.set_page_config(layout='wide')
st.cache_resource.clear()

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
tab1, tab2 = st.tabs(['Browse examples', 'Global statistics'])

with tab1:
    st.write('#### All structures/pockets')
    cols_drop_ = ['pocket_xmin', 'pocket_xmax', 'pocket_ymin', 'pocket_ymax', 'pocket_zmin', 'pocket_zmax']
    df_pockets_ = read_pockets_().drop(cols_drop_, axis=1)\
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

    col1, col2 = st.columns(2)
    with col1:
        st.write(f'#### DeepFRI GO/EC terms for {af2_id_}')
        tab11, tab12 = st.tabs(['Table', 'Heatmap'])
        with tab11:
            #st.markdown(f'**Selected structure**: {af2_id_} in [UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry) / [AlphaFill](https://alphafill.eu/model?id={af2_id_}) / [Ensembl Bacteria](https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q={af2_id_};site=ensemblunit)')
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
                height=300,
                width='100%',
                enable_enterprise_modules=False,
            )
            st.write(f'{af2_id_} in [UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry) / [AlphaFill](https://alphafill.eu/model?id={af2_id_}) / [Ensembl Bacteria](https://bacteria.ensembl.org/Multi/Search/Results?species=all;idx=;q={af2_id_};site=ensemblunit)')

            if len(grid_response['selected_rows']) > 0:
                go_ec_ = grid_response['selected_rows'][0]['GO_term/EC_number']
            else:
                go_ec_ = None
        
        with tab12:
            df_pockets_q_ = df_pockets_.query('UniProtKB_ac == @af2_id_')
            nresid_ = df_pockets_q_['struct_nresid'].drop_duplicates().squeeze()

            df_heatmap_ = pd.DataFrame(index=range(1, nresid_ + 1))
            for i, r in df_pockets_q_.iterrows():
                col_ = f'pocket_id_{r.pocket_id}'
                val_ = ast.literal_eval(r['pocket_resid'])
                df_heatmap_.loc[val_, col_] = 1

            fp_ = f'data/saliency/{af2_id_}_saliency.tsv'
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

        with col_2_:
            value_counts_ = read_deepfri_().query('Score > 0.5 & DeepFRI_ont == @ont_')[['GO_term/EC_number', 'GO_term/EC_number name']].value_counts()
            st.dataframe(value_counts_)

        with col_1_:
            fig, ax = plt.subplots(figsize=(4, 4))
            value_counts_.head(10).plot(kind='barh')
            plt.gca().invert_yaxis()
            st.pyplot(fig)

