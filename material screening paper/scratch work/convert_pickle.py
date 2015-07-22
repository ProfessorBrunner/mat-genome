import pandas as pd

binary_sys = pd.load('binary_sys')
ternary_sys = pd.load('tern_data')

binary_sys.to_csv('binary_sys.csv', index_label = False, header = True)
ternary_sys.to_csv('tern_data.csv', index_label = False, header = True)

