from src.omicscope.Omicscope import *

def test_omicscope_progenesis():
    df = Omicscope('tests//data//proteins//progenesis.csv',
    Method = 'Progenesis', ControlGroup= 'KO')
    assert df.ctrl == 'KO'