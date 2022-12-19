import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

# read vertex.csv with system input
filename = str(sys.argv[1])
df = pd.read_csv(filename, header = None)
df.columns = ["x", "y", "z"]

# plot
plt.plot(df["x"], df["z"])
plt.show()
