import numpy as np
import pandas as pd
import os
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", type=str, default="data_lupivot_cmd.csv", dest="filename", help="File to be analysed", metavar="FILE")
arg = parser.parse_args()

print("File to be analysed:", arg.filename)

os.chdir(os.path.dirname(os.path.abspath(__file__)))
with open(arg.filename, "r") as f:
    datas = f.readlines()


length = len(datas)
size = np.zeros(length)
data = np.zeros(length)
for i, rows in zip(range(length), datas):
    _, size[i], data[i] = re.findall(r'\d+', rows)

size = np.unique(size)
data.resize((len(size), int(len(data) / len(size))))

dmin = np.min(data, axis=1)
dmax = np.max(data, axis=1)
dmean = np.mean(data, axis=1)
dstd = np.std(data, axis=1, ddof=1)
df = pd.DataFrame([size, dmin, dmax, dmean, dstd], index=['Size', 'Min/ns', 'Max/ns', 'Average/ns', 'Std'])
# print(df.iloc[0])
# print(df)
s = df.to_markdown(mode="str", numalign="center", stralign="center").split('\n')
s.remove(s[0])
s[1], s[0] = s[0], s[1]
s = '\n'.join(s)
print(s)
