import numpy as np
import pandas as pd

df = pd.read_csv("meshpoint.csv", header = None)
df = df.iloc[:, :-1]

def to_xyz_string(df):
    xyz_string = ""
    for row in range(0,len(df)):
        
        # number of vertices per iter
        num_vertices = int(len(df.columns)/3)
        xyz_string += str(num_vertices) + "\n"
        
        # "iteration: NUM"
        xyz_string += "iteration: " + str(row) + "\n"
        
        # coordinates "__do + coordinates"
        for i in range(0, int(len(df.columns)/3)):
            xyz_string += "  do"
            
            def format_spaces_decimals(number):
                # 5 decimals, add spaces in front to make 13 char long in total
                if number > 9999:
                    number = 9999.99999
                if number < -9999:
                    number = -9999.99999
                numstr_5deci = "{:.5f}".format(number)
                return " " * (13 - len(numstr_5deci)) + numstr_5deci
            
            for j in range(0,3):
                xyz_string += format_spaces_decimals(df[i*3+j][row])
            
            
            
            xyz_string += "\n"
        print(f"complete:{row} out of {len(df)}")
    return xyz_string
    
f = open("out_trajectory.xyz", "w")
f.write(to_xyz_string(df))
f.close()
