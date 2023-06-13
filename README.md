# Clustering predicted structures at the scale of the known protein universe
## Interactive supplementary
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
$ conda run -p ./streamlit-env streamlit run Main.py
```
