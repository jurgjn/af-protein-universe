
def get_resid_pLDDT(fp_):
    resseq_pLDDT = collections.OrderedDict()
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(fp_, fp_)
    for chains in structure:
        for chain in chains:
            for residue in chain:
                resname = residue.get_resname()
                hetflag, resseq, icode = residue.get_id()
                for atom in residue:
                    resseq_pLDDT[resseq] = atom.bfactor
    return resseq_pLDDT

def get_af2_stats(fp_):
    resseq_pLDDT = get_resid_pLDDT(fp_)
    return len(resseq_pLDDT), np.mean(list(resseq_pLDDT.values()))

rule af2:
    """
        gunzip -c resources/afdb/UP000005640_9606_HUMAN_v3/AF-{wildcards.struct_id}-F1-model_v3.pdb.gz > {output.pdb}
        tar -tf UP000005640_9606_HUMAN_v3.tar --wildcards "*.pdb.gz" > UP000005640_9606_HUMAN_v3.txt
        tar -xf UP000005640_9606_HUMAN_v3.tar --wildcards "*.pdb.gz"
        ---
        tar -xOf resources/alphafold/UP000005640_9606_HUMAN.tar AF-{wildcards.struct_id}-F1-model_v1.pdb.gz | gunzip -c > {output.pdb}
        gunzip -c resources/afdb/UP000005640_9606_HUMAN/AF-{wildcards.struct_id}-F1-model_v1.pdb.gz > {output.pdb}
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v1.pdb ||\
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v2.pdb ||\
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v3.pdb
    """
    output:
        pdb = pfile(struct_id='{}', step='af2', suffix='.pdb'),
    shell: """
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-model_v4.pdb
    """

rule obabel_hxr:
    """
    > The  proteins  from  the  dataset  were  converted  to  the  AutoDock’s PDBQT  format  and  gasteiger  charges  were  assigned.
    > It seems that hydrogen atoms are missing in the receptor.
    obabel for protein prep: http://gnina.github.io/gnina/rsc_workshop2021/#/17
    -xr: https://open-babel.readthedocs.io/en/latest/FileFormats/AutoDock_PDBQT_format.html
    --partialcharge: https://openbabel.org/docs/dev/Command-line_tools/babel.html#options
    """
    input:
        pdb = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdb'),
    output:
        pdb = pfile(struct_id='{}', step='{prev_steps}.obabel_hxr', suffix='.pdbqt'),
    shell: """
        obabel {input.pdb} -h -xr --partialcharge gasteiger -O{output.pdb}
    """

rule obabel_d:
    """
    Remove all hydrogens, see also https://open-babel.readthedocs.io/en/latest/Command-line_tools/babel.html#examples
    """
    input:
        pdb = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdb'),
    output:
        pdb = pfile(struct_id='{}', step='{prev_steps}.obabel_d', suffix='.pdb'),
    shell: """
        obabel {input.pdb} -O{output.pdb} -d
    """

rule autosite:
    """
    Receptor prep as recommended by AutoDockFR:
        https://ccsb.scripps.edu/adfr/how-to-create-a-pdbqt-for-my-receptor/
    export PATH=/hps/nobackup/beltrao/jurgen/ADFR/bin:$PATH
    """
    input:
        pdbqt = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdbqt'),
    output:
        dir = directory(pfile(struct_id='{}', step='{prev_steps}.autosite', suffix='/')),
        autosite_summary = pfile(struct_id='{}', step='{prev_steps}.autosite', suffix='/{struct_id}_summary.csv'),
    shell: """
        conda run -n adfr-suite autosite -r {input.pdbqt} -o {output.dir}
    """

def autosite_residues(fp_struct, fp_cl, radius=4.5):
    """
    https://doi.org/10.1038/s41598-020-72906-7
    > AutoSite 1.1 was run with default settings, and the top predicted binding site cluster was analyzed for each protein structure (XXXX_cl_1.pdb). 
    > Actual predicted binding site residues were back-calculated from these point clusters using a 4.5Å distance cutoff. 
    """
    #print(fp_struct, fp_cl)
    struct = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_struct, fp_struct)
    atoms  = Bio.PDB.Selection.unfold_entities(struct, 'A')
    ns = Bio.PDB.NeighborSearch(atoms)

    # PDBParser does not grok AutoSite output due to duplicate atom names; address by reading manually using read_fwf, and constructing Atom objects adhoc
    #clust = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_cl, fp_cl)
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    # https://gist.github.com/tbrittoborges/929ced78855945f3e296
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26), (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78), (78, 80)]
    names = ['ATOM', 'Atom serial number', 'Atom name', 'Alternate location indicator',
             'Residue name', 'Chain identifier', 'Residue sequence number', 'Insertion code',
             'X', 'Y', 'Z', 'Occupancy', 'Temperature factor', 'Segment identifier', 'Element symbol']
    df_clust = pd.read_fwf(fp_cl, names=names, colspecs=colspecs)
    xmin, xmax, ymin, ymax, zmin, zmax = df_clust['X'].min(), df_clust['X'].max(), df_clust['Y'].min(), df_clust['Y'].max(), df_clust['Z'].min(), df_clust['Z'].max()

    residues = set()
    for i, r in df_clust.iterrows():
        cl_ = [r['X'], r['Y'], r['Z']]
        for res_ in ns.search(center=cl_, radius=radius, level='R'):
            residues.add(res_)

    l_resseq = set()
    for residue in residues:
        hetflag, resseq, icode = residue.get_id()
        l_resseq.add(resseq)

    return xmin, xmax, ymin, ymax, zmin, zmax, l_resseq

def get_lbsp_family_id(struct_id):
    """Get family identifier from the LBSP study"""
    # For four-letter identifiers, assume experimental structure with a pdb_id
    if len(struct_id) == 4:
        df_family_ = resources.family_index().query('PDBid == @struct_id')
    # Otherwise, assume an AF2 structure with a UniProt identifier
    else:
        fp_ = 'results/lbsp_af2_check_resid.tsv'
        df_family_ = pd.read_csv(fp_, sep='\t').query('SP_PRIMARY_SET == @struct_id')[['Family', 'SP_PRIMARY_SET']].drop_duplicates()

    print(struct_id, len(df_family_))
    print(df_family_)

    assert(len(df_family_) == 1)
    return df_family_['Family'].values[0]

def get_ubs_resid(family_id):
    return set(resources.ubs_index().query('Fam_number == @family_id')['Res_number'])

def get_pdb_nresid(fp_):
    # Number of residues in the whole protein (for defining FN)
    d_resseq = dict()
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(fp_, fp_)
    n_resid = 0
    for chains in structure:
        for chain in chains:
            for residue in chain:
                #print(residue)
                n_resid += 1
    return n_resid

def get_lbsp_stats(ghecom_resid, n_resid, ubs_resid):
    TP = len(ghecom_resid & ubs_resid)
    FP = len(ghecom_resid - ubs_resid)
    FN = len(ubs_resid - ghecom_resid)
    TN = n_resid - (TP + FP + FN)
    assert TP >= 0 and FP >= 0 and FN >= 0 and TN >= 0

    if (TP + FP) > 0:
        P = TP / (TP + FP)
    else:
        P = 0
    if (TP + FN) > 0:
        R = TP / (TP + FN)
    else:
        R = 0

    if (P + R) > 0:
        F = 2 * (P * R) / (P + R)
    else:
        F = 0

    if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0:
        MCC = ((TP * TN) - (FP * FN)) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    else:
        MCC = float('nan')

    #print('get_lbsp_stats():', F, MCC, TP, FP, FN, TN)
    #print(sorted(ubs_resid))
    #print(sorted(ghecom_resid))
    print(F, MCC)
    return F, MCC

rule autosite_summary:
    input:
        # https://doi.org/10.1371/journal.pcbi.1006705
        # Hydrogen atoms were not considered in the distance calculation for either the protein or the ligand.
        pdb = pfile(struct_id='{}', step='{step1}', suffix='.pdb'),
        autosite_dir = pfile(struct_id='{}', step='{step1}.{step2}.autosite', suffix='/'),
        autosite_summary = pfile(struct_id='{}', step='{step1}.{step2}.autosite', suffix='/{struct_id}_summary.csv'),
        pdb_dist = pfile(struct_id='{}', step='{step1}.obabel_d', suffix='.pdb'), # Structure with hydrogens explicitly removed for distance calculations
    output:
        tsv = pfile(struct_id='{}', step='{step1,af2|af2_v3|af2_fcz}.{step2,obabel_hxr}.autosite.summary', suffix='.tsv'),
    run:
        df_ = pd.read_csv(input.autosite_summary, sep=',', skiprows=1, names=['pocket_id', 'energy', 'n_points', 'rad_gyration', 'energy_per_vol', 'buriedness', 'score',])#.set_index('cluster_id', drop=True)
        df_.insert(0, 'struct_id', str(wildcards.struct_id))
        df_['cl_file'] = [* map(lambda pocket_id: os.path.join(input.autosite_dir, f'{wildcards.struct_id}_cl_{pocket_id:03d}.pdb'), df_['pocket_id']) ]
        df_['cl_isfile'] = [* map(os.path.isfile, df_['cl_file']) ]
        assert all(df_['cl_isfile'])
        cols_ = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'resid']
        l_bbox_resid_ = [* map(lambda fp_cl: autosite_residues(fp_struct=input.pdb_dist, fp_cl=fp_cl), df_['cl_file']) ] # Do not use hydrogens for calculating associated residues)
        if len(l_bbox_resid_) > 0:
            df_[cols_] = l_bbox_resid_
        else:
            df_ = df_.reindex(df_.columns.tolist() + cols_, axis=1)

        df_['n_resid'] = [* map(len, df_['resid']) ]

        if wildcards.step1 in ['lbsp', 'lbsp_af2']:
            family_id = get_lbsp_family_id(wildcards.struct_id)
            ubs_resid = get_ubs_resid(family_id)
            n_resid = get_pdb_nresid(input.pdb)
            if len(df_) > 0:
                df_[['F', 'MCC']] = [* map(lambda resid: get_lbsp_stats(resid, n_resid, ubs_resid), df_['resid']) ]

        if wildcards.step1 in ['af2', 'af2_fcz', 'af2_v3']:
            def get_mean_pLDDT(l_resid):
                resid_pLDDT_ = get_resid_pLDDT(input.pdb)
                return np.mean([resid_pLDDT_[resid] for resid in l_resid])
            df_['mean_pLDDT'] = [* map(get_mean_pLDDT, df_['resid']) ]
        
        df_['resid'] = [* map(sorted, df_['resid']) ]
        df_.to_csv(output.tsv, sep='\t', index=False)

rule autosite_summary_ext_resid:
    input:
        # https://doi.org/10.1371/journal.pcbi.1006705
        # Hydrogen atoms were not considered in the distance calculation for either the protein or the ligand.
        pdb = pfile(struct_id='{}', step='{step1}', suffix='.pdb'),
        autosite_dir = pfile(struct_id='{}', step='{step1}.{step2}.autosite', suffix='/'),
        autosite_summary = pfile(struct_id='{}', step='{step1}.{step2}.autosite', suffix='/{struct_id}_summary.csv'),
        pdb_dist = pfile(struct_id='{}', step='{step1}.obabel_d', suffix='.pdb'), # Structure with hydrogens explicitly removed for distance calculations
    output:
        tsv = pfile(struct_id='{}', step='{step1,af2|af2_v3|af2_fcz}.{step2,obabel_hxr}.autosite.summary_ext_resid', suffix='.tsv'),
    run:
        df_ = pd.read_csv(input.autosite_summary, sep=',', skiprows=1, names=['pocket_id', 'energy', 'n_points', 'rad_gyration', 'energy_per_vol', 'buriedness', 'score',])#.set_index('cluster_id', drop=True)
        df_.insert(0, 'struct_id', str(wildcards.struct_id))
        df_['cl_file'] = [* map(lambda pocket_id: os.path.join(input.autosite_dir, f'{wildcards.struct_id}_cl_{pocket_id:03d}.pdb'), df_['pocket_id']) ]
        df_['cl_isfile'] = [* map(os.path.isfile, df_['cl_file']) ]
        assert all(df_['cl_isfile'])
        #cols_ = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'resid']
        cols_ = ['resid', 'resid_10', 'resid_15']
        l_bbox_resid_ = [* map(lambda fp_cl: autosite_residues(fp_struct=input.pdb_dist, fp_cl=fp_cl)[-1], df_['cl_file']) ] # Do not use hydrogens for calculating associated residues)
        l_bbox_resid_10_ = [* map(lambda fp_cl: autosite_residues(fp_struct=input.pdb_dist, fp_cl=fp_cl, radius=10)[-1], df_['cl_file']) ] # Do not use hydrogens for calculating associated residues)
        l_bbox_resid_15_ = [* map(lambda fp_cl: autosite_residues(fp_struct=input.pdb_dist, fp_cl=fp_cl, radius=15)[-1], df_['cl_file']) ] # Do not use hydrogens for calculating associated residues)
        if len(l_bbox_resid_) > 0:
            df_['resid'] = l_bbox_resid_
            df_['resid_10'] = l_bbox_resid_10_
            df_['resid_15'] = l_bbox_resid_15_
        else:
            df_ = df_.reindex(df_.columns.tolist() + cols_, axis=1)

        df_['n_resid'] = [* map(len, df_['resid']) ]
        df_['n_resid_10'] = [* map(len, df_['resid_10']) ]
        df_['n_resid_15'] = [* map(len, df_['resid_15']) ]

        if wildcards.step1 in ['lbsp', 'lbsp_af2']:
            family_id = get_lbsp_family_id(wildcards.struct_id)
            ubs_resid = get_ubs_resid(family_id)
            n_resid = get_pdb_nresid(input.pdb)
            if len(df_) > 0:
                df_[['F', 'MCC']] = [* map(lambda resid: get_lbsp_stats(resid, n_resid, ubs_resid), df_['resid']) ]

        if wildcards.step1 in ['af2', 'af2_fcz', 'af2_v3']:
            def get_mean_pLDDT(l_resid):
                resid_pLDDT_ = get_resid_pLDDT(input.pdb)
                return np.mean([resid_pLDDT_[resid] for resid in l_resid])
            df_['mean_pLDDT'] = [* map(get_mean_pLDDT, df_['resid']) ]
        
        df_['resid'] = [* map(sorted, df_['resid']) ]
        df_['resid_10'] = [* map(sorted, df_['resid_10']) ]
        df_['resid_15'] = [* map(sorted, df_['resid_15']) ]
        df_.to_csv(output.tsv, sep='\t', index=False)
