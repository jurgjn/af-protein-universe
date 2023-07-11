# Clustering predicted structures at the scale of the known protein universe
Code, intermediate results and an interactive visualisation on prediction of putative novel enzymes and small molecule binding proteins presented in ([Barrio-Hernandez et al. 2023](https://doi.org/10.1101/2023.03.09.531927)).

- [Predicted pockets with score > 60 and mean pLDDT > 90](pipeline/results/af2_v3.obabel_hxr.autosite.summary.score60_pLDDT90.tsv.gzz)
- [Pocket surfaces, one file per pocket](pipeline/results/af2_v3.obabel_hxr.autosite.summary.score60_pLDDT90)
- [Predicted GO/EC terms for all structuresz](pipeline/results/af2_v3.DeepFRI_terms.tsv.gz)
- [Residue-level saliency weights for MF/EC terms, one file per structure](pipeline/results/af2_v3.DeepFRI_saliency)

## Interactive visualisation
A publicly accessible instance is currently accessible [at this address](https://af-protein-universe.streamlit.app) on the Streamlit Community Cloud.

Alternatively, you install the web app locally by forking the [repository](https://github.com/jurgjn/af-protein-universe) and [setting up a conda environment](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) as follows:
```
$ git clone git@github.com:jurgjn/af-protein-universe.git
$ cd af-protein-universe/
$ conda create -p streamlit-env python numpy matplotlib seaborn 'pandas<2.0.0'
$ conda run -p ./streamlit-env pip install -r requirements.txt
```

After that, run the web app locally with:
```
$ conda run -p ./streamlit-env streamlit run app.py
```
