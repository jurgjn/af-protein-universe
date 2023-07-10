
import ast, random, os, tempfile, time, sqlite3, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, streamlit_ext as ste, st_aggrid, py3Dmol, stmol

st.set_page_config(
    page_title='Clustering predicted structures at the scale of the known protein universe',
    page_icon='🔬',
    layout='centered',
)
#st.cache_resource.clear()

st.markdown("""
# Clustering predicted structures at the scale of the known protein universe
## Interactive supplementary
Choose an analysis to browse from the sidebar on the left or go directly to an example from the figures in the manuscript:
- [A0A849TG76](/Figure_2_Putative_Novel_Enzymes?entryID=A0A849TG76)
and
[A0A2D8BRH7](/Figure_2_Putative_Novel_Enzymes?entryID=A0A2D8BRH7)
from Figure 2B
- [A0A849ZK06](/Figure_2_Putative_Novel_Enzymes?entryID=A0A849ZK06)
from Figure 2C
- [S0EUL8](/Figure_2_Putative_Novel_Enzymes?entryID=S0EUL8)
from Figure 2D
- [A0A2R8Y619](/Figure_3_Human-Centric_Clusters?entryID=A0A2R8Y619)
and
[A0A1G5ASE0](/Figure_3_Human-Centric_Clusters?entryID=A0A1G5ASE0)
from Figure 3C, left
- [B4DKH6](/Figure_3_Human-Centric_Clusters?entryID=B4DKH6)
and
[A0A2D5ZNG0](/Figure_3_Human-Centric_Clusters?entryID=A0A2D5ZNG0)
from Figure 3C, middle
- [O14862](/Figure_3_Human-Centric_Clusters?entryID=O14862)
and
[A0A1C5UEQ5](/Figure_3_Human-Centric_Clusters?entryID=A0A1C5UEQ5)
from Figure 3C, right
""")