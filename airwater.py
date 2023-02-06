from guiqwt.plot import ImageDialog
from guiqwt.builder import make
import numpy as np


def deriveAreaFluence(self, minX, maxX, minY, maxY, N, histories = None, historyType = 'all'):
    def create_window():
        win = ImageDialog(edit=False, toolbar=True, wintitle="Cross Sections Fluence of %s" % self.phspFileName,
                          options=dict(show_xsection=True, show_ysection=True))
        win.resize(600, 600)
        return win

    fluence = np.zeros([N, N])
    deltaX = (maxX - minX) / float(N)
    deltaY = (maxY - minY) / float(N)

    voxel_x_boundary = np.linspace(minX,maxX,N+1)
    voxel_y_boundary = np.linspace(minY,maxY,N+1)

    if histories is None:
        histories = self.NPPHSP

    select_particle = (self.x >= minX) & (self.x <= maxX) & (self.y >=minY) & (self.y <= maxY)

    for i in range(N):
        for j in range(N):
            select_particle = (self.x >= voxel_x_boundary[i]) & (self.x <= voxel_x_boundary[i+1] ) & (self.y >= voxel_y_boundary[j]) & (self.y <= voxel_y_boundary[j+1] )
            fluence[i, j] = np.sum(self.wt[select_particle] * self.energy[select_particle])



if __name__ == '__main__':

    p1.deriveAreaFluence(-20, 20, -20, 20, 20, p1.NPPHSP, 'photon')