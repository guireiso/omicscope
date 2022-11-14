from src.omicscope.Omicscope import *

def test_omicscope_static():
    df = Omicscope('tests//data//proteins//progenesis.csv',
    Method = 'Progenesis', ControlGroup= 'KO')
    assert df.ctrl == 'KO'