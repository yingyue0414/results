# %%
"""
!python 3.8.10
convert_pdb_to_csv_CoM.py
@author Y Ying

This file converts ALL .pdb file in the directory for gag 
lattices output from NERDSS to a .csv of all centers of
mass (CoMs) of gag.

`convert_and_rotate_all()` further rotates the best
fit plane of the gag lattice to xy-plane and place the'
side with more gag CoMs in positive z-direction

Dependencies:
pandas 1.4.3
numpy 1.23.1
matplotlib 3.5.2
scikit-spatial 6.7.0
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from mpl_toolkits import mplot3d
from matplotlib import cm # colormap
from matplotlib.ticker import LinearLocator
from skspatial.objects import Plane, Points
from skspatial.plotting import plot_3d

# %%
def read_pdb(filename):
    """
    Read .pdb file from NERDSS output for gag lattices
    and return dataframe with header renamed. Note the
    last two columns are only placeholder for pdb file
    format.

    Args:
        filename (str): file name string

    Returns:
        df (pandas.DataFrame): dataframe of imported pdb
    """
    df = pd.read_csv(filename, skiprows = 3, sep = "\s+", header = None) #"\s+" for undefined number of spaces
    df.columns = ["Type", "Index", "Particle", "Molecule", "Molecule Index", "X", "Y", "Z", "Occupancy", "Atom"]
    return df

# %%
def get_COM_coord(df, drop = True):
    """
    Extract coordinate data (XYZ) of COM from the input
    dataframe based on "Particle" property. .

    Args:
        df (pandas.DataFrame): dataframe of imported pdb
        drop (bool, optional): Drop all other properties 
        except for CoM coordinates when set to drop = True.
        Defaults to True.

    Returns:
        df_COM (pandas.DataFrame): dataframe of imported pdb
        with only CoM coordinate data.
    """
    df_COM = df.loc[df['Particle'] == "COM"].reset_index(drop = True)
    if drop:
        df_COM = df_COM.drop(columns =\
                             ["Type", "Particle", "Index", "Molecule", "Molecule Index", "Occupancy", "Atom"])
    return df_COM

# %%
def convert_all_to_csv(readDir = "./", saveDir = "./csv/"):
    """
    Iterate over readDir where this file is in, select
    the files that end with .pdb, and convert them to
    .csv with only COM coord. Save .csv in saveDir.

    Args:
        readDir (str, optional): Input file directory. Defaults to "./".
        saveDir (str, optional): Output file directory. Defaults to "./csv/".

    Returns:
        True
    """
    directory = os.fsencode(readDir) # get directory from str
    
    #create save directory if not exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    
    for file in os.listdir(directory):
        #iteratre over directory
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):

            #read files that end with .pdb
            df = read_pdb(filename)
            df = get_COM_coord(df, drop = True) # convert to COM only

            #convert filename to "x.csv"
            filenameCsv = saveDir + filename[:-3] + "csv"

            #write to .csv
            print(f"Write to {filenameCsv}")
            df.to_csv(filenameCsv, header = False, index = False)

        else:
            #skip file
            print(f"Skip {filename}")
    
    return True

# %%
"""
The .pdb file output from NERDSS cannot be input into continuum
membrane model directly because the model requires the input
gag lattice to be in correct orientation with respect to the
membrane to begin with. Therefore, first we need to rotated the
gag lattice so that the plane of interaction is approximately
parallel to xy-plane. This is achieved by:
(1) Calculate the plane of best fit (with scikit-spatial)
(2) Calculate a rotation matrix (Infinite possible solutions;
    only need 1 particular solution in this case) that rotates
    the unit normal vector of the plane of best fit to unit z-vector:
    k = R . n
(3) Rotate all the gags by R
(4) Translate all the gags so that the center of mass is at (x = 0, y = 0)
    and that the point with lowest z has z = 50.0 - lbond
"""

def get_rotation_matrix(vStart, vEnd):
    """
    Get one solution of rotation matrix that rotates vStart to vEnd.

    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    Get rotation matrix based on this method.
    v = a x b
    c = a dot b -- scalar = 1/(1 + c)
    [vx] = [[0, -v[2], v[1]]
            [v[2], 0, -v[0]]
            [-v[1], v[0], 0]]
    [R] = [I] + [vx] + [vx]^2 1/(1 + c)

    Args:
        vStart (numpy.ndarray(3)): vector before rotation
        vEnd (numpy.ndarray(3)): vector after rotation

    Returns:
        RMat (numpy.ndarray((3,3))): rotation matrix that rotates
        vStart to vEnd
    """
    # v = a x b
    v = np.cross(vStart, vEnd)
    # c = a dot b
    c = np.dot(vStart, vEnd)
    scalar = 1.0 / (1 + c)
    #[vx] = [[0, -v[2], v[1]]
    #   [v[2], 0, -v[0]]
    #   [-v[1], v[0], 0]]
    vxMat = np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])
    # [R] = [I] + [vx] + [vx]^2 1/(1 + c)
    RMat = np.identity(3) + vxMat + np.dot(vxMat, vxMat) * scalar
    return RMat

# %%
def calculate_rotation_min_max(points, RMat):
    """
    Internal method used by convert_and_rotate_all;
    Iterate over points and calculate rotated. Enumerate
    and calculate COM and subtract r_COM from all points
    to zero COM. Calculate minimum and maximum z-coordinate.

    Args:
        points (numpy.ndarray(N, 3)): gag coordinates
        RMat (numpy.ndarray(3, 3)): rotation matrix

    Returns:
        pointsRotated: rotated points with zeroed COM
        minZ: minimum z-coord of rotated points
        maxZ: maximum z-coord of rotated points
    """
    pointsRotated = np.zeros(np.shape(points)) # initialize
    sumCOM = np.zeros([3]) # sum for zeroing CoM
    for (pointStart, pointRotated) in zip(points, pointsRotated):
        pointRotated += np.dot(RMat, pointStart) # MUST USE += here
        sumCOM += pointRotated
            
            # zero CoM
            
    avgCOM = sumCOM / (np.shape(pointsRotated)[0])
    for pointRotated in pointsRotated:
        pointRotated -= avgCOM
            
            # the membrane is default to be at +50.0 (above zero to make sure curvature sign)
            # (is not ambiguous)
            # for lbond = 12.0, the lowest point is set at +36.0
            # find min and max in z-dir
    minZ = 99999.0
    maxZ = -99999.0
    for pointRotated in pointsRotated:
        if pointRotated[2] < minZ:
            minZ = pointRotated[2]
        if pointRotated[2] > maxZ:
            maxZ = pointRotated[2]
    return pointsRotated, minZ, maxZ
    

def convert_and_rotate_all(readDir = "./", saveDir = "./csv/", plot = False):
    """
    Combine conversion to .csv and rotation. Rotate gag
    lattice such that the normal vector of best fitting
    plane of the resulting gag lattice is pointing in
    positive z-dir. Output .csv file.
    
    Dependency: scikit-spatial

    Args:
        readDir (str, optional): Input file directory. Defaults to "./".
        saveDir (str, optional): Output file directory. Defaults to "./csv/".
        plot (bool, optional): Plot scatter in output file directory. Defaults to False.

    Returns:
        True
    """
    directory = os.fsencode(readDir)
    
    #create save directory if not exists
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    
    for file in os.listdir(directory):
        #iteratre over directory
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):

            #read files that end with .pdb
            print(f"Read {filename}") 
            df = read_pdb(filename)
            df = get_COM_coord(df, drop = True)
            
            #get best fit plane
            points = Points(df.to_numpy())
            plane = Plane.best_fit(points)
            
            #calculate rotation matrix RMat
            vStart = plane.normal
            vEnd = np.array([0.0, 0.0, 1.0]) #positive z-dir
            RMat = get_rotation_matrix(vStart, vEnd)

            # iterator over Points and calculate rotated; calculate sum for zeroing CoM
            pointsRotated, minZ, maxZ = calculate_rotation_min_max(points, RMat)

            # rotate by 180 deg if the side with more points is facing downwards
            # this make sure that the "smoother" side of the hemisphere is facing
            # upwards (positive z-direction)
            planeRotated = Plane.best_fit(pointsRotated) #get best fit plane of rotated gag lattice
            zPlaneRotated = planeRotated.point[2] #get z-coord of plane
            nPtWrtPlane = [0, 0] # number of point with respect to plane -> [above, below]
            for pt in pointsRotated:
                if pt[2] >= zPlaneRotated:
                    nPtWrtPlane[0] += 1 # enumerate pt above plane
                else:
                    nPtWrtPlane[1] += 1 # enumerate pt below plane
            if nPtWrtPlane[0] < nPtWrtPlane[1]: # rotate 180 if npt above plane > npt below plane
                vEndRev = np.array([0.0, 0.0, -1.0])
                RMatRev = get_rotation_matrix(vStart, vEndRev)
                pointsRotated, minZ, maxZ = calculate_rotation_min_max(points, RMatRev)

            if (maxZ - minZ > 15.0):
                print("WARNING: Z-DIR DISTANCE OVER 15.0!")
            # set lowest point to +35.0
            moveDistance = 38.0 - minZ
            for pt in pointsRotated:
                pt[2] += moveDistance
            
            #get dataframe
            df_pointsRotated = pd.DataFrame(pointsRotated)
            df_pointsRotated.columns = (["x", "y", "z"])

            # plot scatter3d if plot set to true
            if plot:

                fig = plt.figure(figsize=(24, 12))
                ax = plt.axes(projection ="3d")

                # Plot the surface.
                ax.scatter3D(df_pointsRotated["x"], df_pointsRotated["y"], df_pointsRotated["z"])
                plt.savefig(saveDir + filename[:-3] + "_scatter.png", dpi = 500)


            #convert filename to "x.csv"
            filenameCsv = saveDir + filename[:-3] + "csv"

            #write to .csv
            print(f"Write to {filenameCsv}")
            df_pointsRotated.to_csv(filenameCsv, header = False, index = False)

        else:
            #skip file
            print(f"Skip {filename}")
    
    return True