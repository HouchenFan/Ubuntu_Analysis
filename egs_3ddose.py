#!/userenv/anaconda3/bin/python

import sys, re, os.path
import numpy
import matplotlib.pyplot as plt


# function to find minimum and maximum position in list
def minimum(a, n):
    # inbuilt function to find the position of minimum
    minvalue = min(a)
    minpos = a.index(min(a))

    # inbuilt function to find the position of maximum
    maxvalue = max(a)
    maxpos = a.index(max(a))

    # printing the position
    print("The maximum is ", maxvalue, " at position", maxpos + 1)
    print("The minimum is ", minvalue, " at position", minpos + 1)


if (len(sys.argv) < 2):
    print("This can be used both ways : using command line arguments or interactive mode")
    print('''\ncommand line usage: %s  x|y|z  a,b  filename

    x|y|z       axis of the dose profile
    a,b         coordinates along the other axes
    filename    3ddose file to plot without extension

    example:   egs-doseplot.py z 0,0 filename
        (extract the central z-axis dose profile for the filename.3ddose file)
        ''' % (os.path.basename(sys.argv[0])))
    files = input("Enter the 3ddose File name without extension or type exit to quit :")
    if (files == "exit") or (files == ""):
        sys.exit(1)
    files = files + ".3ddose"
    dosefile = open(files, 'r')
    # get voxel counts on first line
    nx, ny, nz = list(map(int, dosefile.readline().split()))  # number of voxels along x, y, z
    print("The dimesion of dose matrix x,y,z :", nx, ny, nz)
    dosefile.close()
    var = input("Enter the axis for plotting(x/y/z) :")
    if var in ['x', 'y', 'z']:
        axis = ['x', 'y', 'z'].index(var)  # axis number: x=0, y=1, z=3
    else:
        axis = 2  # assume z axis plot
        print("No axis selected. z axis assumed ")
    where = input("Enter the values of other axis(a,b) :")
    if (where != ""):
        where = list(map(float, where.split(',')))  # position of the axis in x,y form
    else:
        where = [0, 0]  # assume 0,0 dd curve
        print("No coordinates for other axis selected. 0,0 assumed")

else:
    # simplistic command-line parsing
    var = sys.argv[1]  # arg1 is the independent variable: x, y, or z
    axis = ['x', 'y', 'z'].index(var)  # axis number: x=0, y=1, z=3
    where = list(map(float, sys.argv[2].split(',')))  # position of the axis in x,y form
    files = sys.argv[3] + ".3ddose"  # array of all files to process

# open 3ddose file
dosefile = open(files, 'r')
# get voxel counts on first line
nx, ny, nz = list(map(int, dosefile.readline().split()))  # number of voxels along x, y, z
print("The dimesion of dose matrix x,y,z :", nx, ny, nz)
Ng = (nx + 1) + (ny + 1) + (nz + 1)  # total number of voxel grid values (voxels+1)
Nd = nx * ny * nz  # total number of data points

# get voxel grid, dose and relative errors
data = list(map(float, dosefile.read().split()))  # read the rest of the file
xgrid = data[:nx + 1]  # voxel boundaries in x (nx+1 values, 0 to nx)
ygrid = data[nx + 1:nx + 1 + ny + 1]  # voxel boundaries in y (ny+1 values, nx+1 to nx+1+ny)
zgrid = data[nx + 1 + ny + 1:Ng]  # voxel boundaries in z (nz+1 values, rest up to Ng-1)
dose = data[Ng:Nd + Ng]  # flat array of Nd = nx*ny*nz dose values
errs = data[Nd + Ng:]  # flat array of Nd = nx*ny*nz relative error values
del data  # make sure we don't refer to data by mistake from now on
# close 3ddose file
dosefile.close()
print('max and min within the whole matrix :\n')
minimum(dose, len(dose))

# setup variables for current axis
grid = [xgrid, ygrid, zgrid]  # voxel boundaries in x, y, z
step = [1, nx, nx * ny]  # step between values along x, y, z
jump = [nx, nx * ny, nx * ny * nz]  # distance between start and end voxels
mygrid = grid[axis]  # grid for plot axis
mystep = step[axis]  # step for plot axis
del grid[axis]  # remove plot axis from grid
del step[axis]  # remove plot axis from step
# print xmgrace header
print('''\
@   title "dose distribution"
@   xaxis label "%s (cm)"
@   yaxis label "dose (Gy)"
@   type xydy''' % var)
# get voxel indices for location (along two remaining axes)
# (say you are plotting along z, then index will hold the indices for the requested x,y positions
index = []
for g, w in zip(grid, where):
    if (w < g[0] or w > g[-1]):
        print("ERROR: location", where, "outside of data grid!\n")
        sys.exit(1)
    i = len(g) - 1
    while (w < g[i]): i -= 1
    index.append(i)

# get the actual dose profiles
start = index[0] * step[0] + index[1] * step[1]  # starting voxel index
end = start + jump[axis]  # "end" voxel index
mydose = dose[start:end:mystep]  # dose slice
myerrs = errs[start:end:mystep]  # relative error slice
print('max and min within the chosen plane :\n')
minimum(dose[start:end:mystep], len(dose[start:end:mystep]))
# print xmgrace format commands for this set
print('''\
@   s%d errorbar length 0
@   s%d legend "%s" ''' % (1, 1, files))
# Numpy tools
te = numpy.empty([len(mydose)])
errors_ = numpy.empty([len(myerrs)])

# print dose profile with absolute errors
print("#\n# %10s %12s %12s" % (var, "dose", "err"))
for i in range(len(mydose)):
    t = (mygrid[i] + mygrid[i + 1]) / 2.0
    te[i] = t
    print("%12.6f %12g %12g" % (t, mydose[i], myerrs[i] * mydose[i]))
    errors_[i] = myerrs[i] * mydose[i]
print

# count += 1
plt.errorbar(te, mydose, yerr=errors_, fmt='b-.', ecolor='g')

plt.title("dose distribution")
plt.xlabel("position (cm)")
plt.ylabel("dose (Gy)")
plt.show()