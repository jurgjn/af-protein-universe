import streamlit as st

st.set_page_config(
    page_title='Clustering predicted structures at the scale of the known protein universe',
    page_icon='🔬',
)

with open('README.md') as fh:
    st.markdown(fh.read())

st.markdown("""
Browse an analysis from the sidebar on the left or go directly to an example from the figures:
- [A0A849TG76](/Figure_2_Putative_Novel_Enzymes?UniProtKB_ac=A0A849TG76)
and
[A0A2D8BRH7](/Figure_2_Putative_Novel_Enzymes?UniProtKB_ac=A0A2D8BRH7)
from Figure 2B
- [A0A849ZK06](/Figure_2_Putative_Novel_Enzymes?UniProtKB_ac=A0A849ZK06)
from Figure 2C
- [S0EUL8](/Figure_2_Putative_Novel_Enzymes?UniProtKB_ac=S0EUL8)
from Figure 2D
- [A0A2R8Y619](/Figure_3_Human-Centric_Clusters?UniProtKB_ac=A0A2R8Y619)
and
[A0A1G5ASE0](/Figure_3_Human-Centric_Clusters?UniProtKB_ac=A0A1G5ASE0)
from Figure 3C, left
- [O14862](/Figure_3_Human-Centric_Clusters?UniProtKB_ac=O14862)
and
[A0A1C5UEQ5](/Figure_3_Human-Centric_Clusters?UniProtKB_ac=A0A1C5UEQ5)
from Figure 3C, right
""")