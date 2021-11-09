from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--file", type=str, default="data_lupivot_cmd.csv", dest="filename", help="File to be analysed", metavar="FILE")
arg = parser.parse_args()

fig, ax = plt.subplots()
plt.xlabel("Size")
plt.ylabel("Average Time/ns")
ax.grid()

print("File to be analysed:", arg.filename)

os.chdir(os.path.dirname(os.path.abspath(__file__)))
for filename in arg.filename.split(','):
    with open(filename, "r") as f:
        f.readline()
        datas = f.readlines()

    length = len(datas)
    size = np.zeros(length)
    data = np.zeros((length, len(datas[0].split(',')) - 2))
    for i, rows in zip(range(length), datas):
        rows = rows.split(',')
        size[i] = rows[0]
        data[i, ] = rows[1:-1]

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

    plt.plot(size, dmean)


ax.legend(['O2', 'O3'])
ax.set_title("Comparison of Different Optimization")
image_format = 'svg'  # e.g .png, .svg, etc.
image_name = 'op.svg'

fig.savefig(image_name, format=image_format, dpi=1200, transparent=True)
