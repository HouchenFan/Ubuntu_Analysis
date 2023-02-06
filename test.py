import numpy as np
from datetime import datetime

path_phsp = R"D:\Data_fhc\EGSnrc\DOS_FS10_SW15\EGS.IAEAphsp"
with open(path_phsp, 'rb') as fid1:
    startTime = datetime.now()
    file1 = fid1.read()
    endTime = datetime.now()
    print("fileRead Time is {:d} seconds".format((endTime - startTime).seconds))
fid1.close()

with open(path_phsp, 'rb') as fid2:
    startTime = datetime.now()
    file2 = np.fromfile(path_phsp, dtype=float, count=-1, sep="", offset=0)
    endTime = datetime.now()
    print("numpyFromfile Time is {:d} seconds".format((endTime - startTime).seconds))
fid2.close()
