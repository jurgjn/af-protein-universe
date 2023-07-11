
rule DeepFRI:
    input:
        pdb = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdb'),
    output:
        # saliency map files don't seem to be generated in some situations (ontologies?)
        _BP_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_pred_scores.json'),
        _BP_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_predictions.csv'),
        #_BP_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_saliency_maps.json'),
        _CC_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_pred_scores.json'),
        _CC_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_predictions.csv'),
        #_CC_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_saliency_maps.json'),
        _EC_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_pred_scores.json'),
        _EC_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_predictions.csv'),
        _EC_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_saliency_maps.json'),
        _MF_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_pred_scores.json'),
        _MF_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_predictions.csv'),
        _MF_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_saliency_maps.json'),
    params:
        output_fn_prefix = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}'),
    shell: """
        cd software/DeepFRI
        conda run -n DeepFRI-env python predict.py -pdb ../../{input.pdb} -ont bp --saliency -v --output_fn_prefix ../../{params.output_fn_prefix}
        conda run -n DeepFRI-env python predict.py -pdb ../../{input.pdb} -ont cc --saliency -v --output_fn_prefix ../../{params.output_fn_prefix}
        conda run -n DeepFRI-env python predict.py -pdb ../../{input.pdb} -ont ec --saliency -v --output_fn_prefix ../../{params.output_fn_prefix}
        conda run -n DeepFRI-env python predict.py -pdb ../../{input.pdb} -ont mf --saliency -v --output_fn_prefix ../../{params.output_fn_prefix}
        cd ../..
    """

rule DeepFRI_summary:
    input:
        pdb = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdb'),
        _BP_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_pred_scores.json'),
        _BP_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_predictions.csv'),
        _CC_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_pred_scores.json'),
        _CC_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_predictions.csv'),
        _EC_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_pred_scores.json'),
        _EC_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_predictions.csv'),
        _MF_pred_scores = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_pred_scores.json'),
        _MF_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_predictions.csv'),
    output:
        tsv_predictions = pfile(struct_id='{}', step='{prev_steps}.DeepFRI.summary', suffix='/{struct_id}_predictions.tsv'),
        tsv_saliency = pfile(struct_id='{}', step='{prev_steps}.DeepFRI.summary', suffix='/{struct_id}_saliency.tsv'),
    params:
        #_BP_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_BP_saliency_maps.json'),
        #_CC_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_CC_saliency_maps.json'),
        _EC_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_EC_saliency_maps.json'),
        _MF_saliency_maps = pfile(struct_id='{}', step='{prev_steps}.DeepFRI', suffix='/{struct_id}_MF_saliency_maps.json'),
    run:
        df_ec_ = pd.read_csv(input._EC_predictions, comment='#', sep=',').assign(DeepFRI_ont='EC')
        df_ec_['GO_term/EC_number'] = 'EC:' + df_ec_['GO_term/EC_number']
        df_ = pd.concat([
            pd.read_csv(input._BP_predictions, comment='#', sep=',').assign(DeepFRI_ont='BP'),
            pd.read_csv(input._CC_predictions, comment='#', sep=',').assign(DeepFRI_ont='CC'),
            pd.read_csv(input._MF_predictions, comment='#', sep=',').assign(DeepFRI_ont='MF'),
            df_ec_,
        ])
        df_['Protein'] = wildcards.struct_id
        df_.to_csv(output.tsv_predictions, header=True, index=False, sep='\t')        

        #d_ = df_[['GO_term/EC_number', 'GO_term/EC_number name']].set_index('GO_term/EC_number', drop=True).squeeze().to_dict()
        d_ = {r['GO_term/EC_number']: r['GO_term/EC_number name'] for i, r in df_.iterrows()}

        def read_saliency_map(fp_):
            with open(fp_) as fh:
                yaml_ = yaml.load(fh, Loader=yaml.SafeLoader)
                #print(fp_, yaml_, os.path.getsize(fp_))
                return pd.DataFrame.from_dict({go_id if go_id.startswith('GO:') else f'EC:{go_id}': saliency_map for go_id, go_names, saliency_map in zip(yaml_['query_prot']['GO_ids'], yaml_['query_prot']['GO_names'], yaml_['query_prot']['saliency_maps']) })
        l_ = [* map(read_saliency_map, [fp for fp in params if os.path.isfile(fp) and os.path.getsize(fp) > 2]) ] # getsize() > 2 checks for empty output
        if len(l_) > 0:
            df_ = pd.concat(l_, axis=1)
        else:
            df_ = pd.DataFrame()
        df_.index = range(1, len(df_) + 1)
        df_.to_csv(output.tsv_saliency, header=True, index=False, sep='\t')        

        # Write .pdb files with saliency as B-factor
        parser = Bio.PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(input.pdb, input.pdb)

        def write_saliency(structure, d_saliency, fp_out):
            if set(d_saliency.values()) == set(['NaN']): return # Check for all-NaN saliency
            for chains in structure:
                for chain in chains:
                    for residue in chain:
                        resname = residue.get_resname()
                        hetflag, resseq, icode = residue.get_id()
                        for atom in residue:
                            atom.set_bfactor(d_saliency[resseq])

            io = Bio.PDB.PDBIO()
            io.set_structure(structure)
            io.save(fp_out)
            #print(fp_out)

        for col_ in df_.columns:
            fp_ = pfile(
                struct_id=wildcards.struct_id,
                step=f'{wildcards.prev_steps}.DeepFRI.summary',
                suffix=f'/{wildcards.struct_id}_{penc(col_)}_{penc(d_[col_])[:50]}.pdb',
                base=wildcards.base,
            )
            write_saliency(structure, df_[col_].squeeze().to_dict(), fp_)
