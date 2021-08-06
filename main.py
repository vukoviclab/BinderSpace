import pandas as pd
from src._sequence import Sequence

df = pd.read_csv('./keat.csv')
Keat = Sequence(df)
print(Keat.df)
