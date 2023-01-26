import pandas as pd
from os import path
from constants import ASSETS_FOLDER_NAME

def get_data_to_fit(filename, dirpath=None):
    filepath = filepath_or_default(filename, dirpath)
    
    df = pd.read_excel(
    filepath, 
    skiprows=[0,1], 
    header=None, 
    index_col=0,
    )
    
    df.columns = ["mh", "rhoh", "Cph", "muh", "kh", "Thi", "Tho", "mc", "rhoc", "Cpc", "muc", "kc", "Tci", "Tco"]
    df.index.name = "t"
    
    return df


def filepath_or_default(filename, dirpath=None):
    dirpath = dirpath or path.join(path.dirname(path.dirname(__file__)), ASSETS_FOLDER_NAME)
    filepath = path.abspath(path.join(dirpath, filename))
    return filepath