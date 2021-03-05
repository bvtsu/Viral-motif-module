import sys
import ntpath
import streamlit as st
from upsetplot import from_memberships, plot
from matplotlib import pyplot as plt

@st.cache
def line_counter(filepath):
    with open(filepath) as f:
        line_count = 0
        for line in f:
            line_count += 1
    return(line_count)

def _max_width_():
    max_width_str = f"max-width: 2000px;"
    st.markdown(
        f"""
    <style>
    .reportview-container .main .block-container{{
        {max_width_str}
    }}
    </style>    
    """,
        unsafe_allow_html=True,
    )

list1 = ntpath.basename("{0}".format(sys.argv[1])).split("_")[0]
list2 = ntpath.basename("{0}".format(sys.argv[2])).split("_")[0]
upset_data = from_memberships(
    [[list1],
    [list2],
    [list1,list2],
    ],
    data=[line_counter(sys.argv[1]),line_counter(sys.argv[2]),line_counter(sys.argv[3])]
)
st.set_option('deprecation.showPyplotGlobalUse', False)
plot(upset_data)
plt.suptitle('{0}'.format(sys.argv[4]))
st.pyplot()