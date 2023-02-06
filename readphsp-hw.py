# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
matplotlib.rcParams['backend'] = 'SVG'
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import struct
# from guiqwt.plot import ImageDialog
# from guiqwt.builder import make
import math
import pandas as pd
import os
# import time

class History(object):
    """docstring for History"""

    def __init__(self, binData):
        super(History, self).__init__()
        self.size = len(binData)
        if self.size == 28:
            self.MODE_RW = 'MODE0'
        elif self.size == 32:
            self.MODE_RW = 'MODE2'
        else:
            assert False, 'The PhaseSpace File is neither mode 2 nor mode 0'
        self.dataType = 'ifffffff' if self.MODE_RW == 'MODE2' else 'iffffff'
        tmp = struct.unpack(self.dataType, binData)
        self.LATCH = tmp[0]
        self.E = tmp[1] if tmp[1] > 0 else -tmp[1]
        assert self.E >= 0, self.E
        self.X = tmp[2]
        self.Y = tmp[3]
        self.U = tmp[4]
        self.V = tmp[5]
        # print(f'self.U={self.U},self.V={self.V}')
        self.W = math.sqrt(1 - self.U ** 2 - self.V ** 2)
        self.WT = tmp[6]
        if self.MODE_RW == 'MODE2':
            self.ZLAST = tmp[7]
        tmp1 = self.LATCH >> 29
        if tmp1 == 0:
            self.IQ = 0
        elif tmp1 == 1:
            self.IQ = 1
        else:
            self.IQ = -1

    def pack(self):
        if self.MODE_RW == 'MODE2':
            return struct.pack(self.dataType, self.LATCH, self.E, self.X, self.Y,
                               self.U, self.V, self.WT, self.ZLAST)
        else:
            return struct.pack(self.dataType, self.LATCH, self.E, self.X, self.Y,
                               self.U, self.V, self.WT)

    def __str__(self):
        if self.MODE_RW == 'MODE2':
            return '%f,%f,%f,%f,%f,%f,%d\n' % (self.E, self.X, self.Y, self.U,
                                               self.V, self.WT, self.LATCH)
        else:
            return '%f,%f,%f,%f,%f,%f,%f,%d\n' % (self.E, self.X, self.Y, self.U,
                                                  self.V, self.WT, self.ZLAST, self.LATCH)


class PhaseSpace(object):
    """docstring for PhaseSpace"""

    def __init__(self, fileName):
        self.phspFileName = fileName
        if fileName.endswith('.IAEAphsp'):
            self.header_file = fileName.replace('.IAEAphsp','.IAEAheader')
        else:
            self.header_file = None


    def readphsp(self):
        headType = 'iifff'
        self.historyTypeSelection = {'all': {-1, 1, 0}, 'electron': {-1}, 'photon': {0}}
        with open(self.phspFileName,'rb') as fid:

            self.MODE_RW = fid.read(5).decode('utf-8')
            assert self.MODE_RW == 'MODE2' or self.MODE_RW == 'MODE0', self.MODE_RW
            self.phspHead = fid.read(4 * 5)
            self.Nothing = fid.read(3) if self.MODE_RW == 'MODE0' else fid.read(7)
            self.NPPHSP, self.NPHOTPHSP, self.EKMAXPHSP, self.EKMINPHSPE, self.NINCPHSP = struct.unpack(headType,self.phspHead)

            type = np.dtype(
                [('LATCH', 'i4'), ('Energy', 'f4'), ('X', 'f4'), ('Y', 'f4'), ('U', 'f4'), ('V', 'f4'), ('WT', 'f4')])
            all_particle_infor = np.fromfile(self.phspFileName, dtype=type, count=-1, sep="", offset=28)
            self.latch = all_particle_infor['LATCH']
            self.energy = all_particle_infor['Energy']
            self.x = all_particle_infor['X']
            self.y = all_particle_infor['Y']
            self.u = all_particle_infor['U']
            self.v = all_particle_infor['V']
            self.w =  np.sqrt(1.0 - np.square(all_particle_infor['U']) - np.square(all_particle_infor['V']))
            self.wt = all_particle_infor['WT']

            self.particle_type = np.round(self.latch/1e9) # 1：electron  0：photon
            self.create_positoin = np.floor(self.latch/16777216)
        print('Finished reading phash space')

    def readIAEAphsp(self):


        with open(self.header_file, 'r', encoding='utf-8') as fid_header:     # 读取IAEAheader文件基本信息   add by FHC 2022/9/7
            self.NINCPHSP = 0
            while (self.NINCPHSP == 0):
                lines = fid_header.readline()  #整行读取数
                if ("$ORIG_HISTORIES:") in lines:
                    self.NINCPHSP = int(next(fid_header).replace("\n"," "))     #History 数量
                # elif ("$PARTICLES:") in lines:
                #     NPPHSP = int(next(fid_header).replace("\n", " "))      #所有的粒子数量
                # elif ("$PHOTONS:") in lines:
                #     NPHOTPHSP = int(next(fid_header).replace("\n", " "))       #所有的光子数量


        type = np.dtype(
            [('particle_type', 'S1'), ('Energy', 'f4'), ('X', 'f4'), ('Y', 'f4'), ('Z', 'f4'), ('U', 'f4'), ('V', 'f4'),
             ('WT', 'f4'), ('LATCH', 'i4')])  # 定义每个粒子的数据类型


        self.TotalNumParticles = int(os.path.getsize(self.phspFileName) / 33)

        Raw = np.fromfile(self.phspFileName, dtype=type, count=-1, sep="", offset=0)  # 按照数据类型读取全部粒子信息
        PtypeIndex = np.where(Raw['particle_type'] == b'\xff')[0]  # 找到Z方向速度为负的粒子的index
        WIndex = np.where((1 - np.square(Raw['U']) - np.square(Raw['V'])) <= 0)[0]  # 找到W^2小于等于0的粒子的index
        Raw_select = np.delete(Raw, list(set(PtypeIndex).union(set(WIndex))), axis=0)  # 删除上述两种粒子的数据

        PhotonIndex = np.where(Raw_select['particle_type'] == b'\x01')[0]  # 找到光子的index
        ElectronIndex = np.where(Raw_select['particle_type'] == b'\x02')[0]  # 找到电子的index
        PositronIndex = np.where(Raw_select['particle_type'] == b'\x03')[0]  # 找到正电子的index
        self.NPPHSP = len(Raw_select)  # 总粒子数
        self.NPHOTPHSP = len(PhotonIndex)  # 总光子数
        self.EKMAXPHSP = max(abs(Raw_select['Energy']))  # 所有粒子最大动能
        self.EKMINPHSPE = min(abs(Raw_select['Energy'][ElectronIndex]))  # 所有电子最小动能
        # 将电子的动能转化为总能量同时保留其正负号
        Raw_select['Energy'][
            list(set(np.where(Raw_select['Energy'] >= 0)[0]).intersection(set(ElectronIndex)))] += 0.511
        Raw_select['Energy'][
            list(set(np.where(Raw_select['Energy'] < 0)[0]).intersection(set(ElectronIndex)))] -= 0.511
        Raw_select['LATCH'][ElectronIndex] += 2 << 29  # 在LATCH中添加Charge信息
        # 将正电子的动能转化为总能量同时保留其正负号
        Raw_select['Energy'][
            list(set(np.where(Raw_select['Energy'] >= 0)[0]).intersection(set(PositronIndex)))] += 0.511
        Raw_select['Energy'][
            list(set(np.where(Raw_select['Energy'] < 0)[0]).intersection(set(PositronIndex)))] -= 0.511
        Raw_select['LATCH'][PositronIndex] += 1 << 29  # 在LATCH中添加Charge信息
        # 将粒子输运到EPID的上平面，即修改XY坐标
        # WList = np.sqrt(1.0 - np.square(Raw_select['U']) - np.square(Raw_select['V']))
        # Raw_select['X'] += (Raw_select['U'] / WList) * (zEpid - Raw_select['Z'])
        # Raw_select['Y'] += (Raw_select['V'] / WList) * (zEpid - Raw_select['Z'])


        self.latch = Raw_select['LATCH']
        self.energy = Raw_select['Energy']
        self.x = Raw_select['X']
        self.y = Raw_select['Y']
        self.u = Raw_select['U']
        self.v = Raw_select['V']
        self.w = np.sqrt(1.0 - np.square(Raw_select['U']) - np.square(Raw_select['V']))
        self.wt = Raw_select['WT']

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
        index_x = np.array([int(i) for i in (self.x[select_particle] - minX) / deltaX])
        index_y = np.array([int(i) for i in (self.y[select_particle] - minY) / deltaY])
        arr = np.vstack((index_x, index_y, self.wt[select_particle], self.energy[select_particle])).transpose()
        df = pd.DataFrame(arr, columns=['x_index', 'y_index', 'wt', 'energy'])

        def getFluence(df):
            # 分组后的df，在相同的voxel
            df['fluence'] = np.sum( df['wt'] * df['energy'] )
            return df.iloc[0][['x_index', 'y_index','fluence']]

        grouped_data = df.groupby(['x_index', 'y_index']).apply(getFluence)
        for i in range(N):
            for j in range(N):
                if (i,j) in grouped_data.index:
                    fluence[i, j] = grouped_data.loc[(i,j)]['fluence']

        fluence = fluence.transpose()
        # for i in range(N):
        #     for j in range(N):
        #         select_particle = (self.x >= voxel_x_boundary[i]) & (self.x <= voxel_x_boundary[i+1] ) & (self.y >= voxel_y_boundary[j]) & (self.y <= voxel_y_boundary[j+1] )
        #         fluence[i, j] = np.sum(self.wt[select_particle] * self.energy[select_particle])
        # X, Y = np.array([minX + deltaX / 2 * (2 * k + 1) for k in range(N)]), np.array([minY + deltaY / 2 * (2 * k + 1) for k in range(N)])
        import guidata
        _app = guidata.qapplication()
        win = create_window()
        image = make.image(data = fluence, xdata = [minX,maxX], ydata = [minY,maxY])
        plot = win.get_plot()
        plot.add_item(image)
        win.exec_()





    def __iter__(self):
        self.reload()
        outLoop = (self.NPPHSP - 1) // (self.cacheSize // self.binSize) + 1
        # self.log()
        count = 0
        for j in range(outLoop):
            tmp = self.phspFile.read(self.cacheSize)
            num = int(len(tmp) / self.binSize)
            assert len(tmp) % self.binSize == 0, len(tmp)
            for i in range(num):
                history = History(tmp[i * self.binSize:(i + 1) * self.binSize])
                count += 1
                yield history
        assert count == self.NPPHSP, 'count不等于self.NPPHSP'
        self.phspFile.close()



    def deriveSpectral(self, minEnergy=None, maxEnergy=None, N=100, historyType='all', histories=None):
        pass








    def deriveFluence(self, minX, maxX, minY, maxY, N, histories=None, historyType='all'):
        self.reload()
        fluence = np.zeros([N, N])
        deltaX = (maxX - minX) / float(N)
        deltaY = (maxY - minY) / float(N)
        index = 0
        if histories is None:
            histories = 1000000 if self.NPPHSP > 1000000 else self.NPPHSP
        for history in self:
            if history.IQ in self.historyTypeSelection[historyType]:
                if minX < history.X < maxX and minY < history.Y < maxY:
                    fluence[
                        int((history.X - minX) / deltaX), int((history.Y - minY) / deltaY)] += history.WT * history.E
                    index += 1
                    if index == histories:
                        break

        X, Y = np.array([minX + deltaX / 2 * (2 * k + 1) for k in range(N)]), np.array(
            [minY + deltaY / 2 * (2 * k + 1) for k in range(N)])
        self.close()
        return X, Y, fluence





    def deriveCircleEnergyFluence(self, N, Rmax, historyType='all'):
        circleArea = Rmax ** 2 * math.pi / N
        index2Radius = lambda i: (math.sqrt(i * circleArea / math.pi) + math.sqrt((i + 1) * circleArea / math.pi)) / 2
        circleEnergyFluence = np.zeros(N, dtype=float)
        circleRadius = np.array([index2Radius(i) for i in range(N)])
        for history in self:
            if history.IQ in self.historyTypeSelection[historyType]:
                index = int((history.X ** 2 + history.Y ** 2) * math.pi / circleArea)
                assert index >= 0, index
                if index < N:
                    if history.IQ == 0:
                        tmp = history.E
                    else:
                        assert history.E >= 0.511, history.E
                        tmp = history.E - 0.511
                    w = math.sqrt(1 - history.U ** 2 - history.V ** 2)
                    if w == 0:
                        w = 0.08716
                    circleEnergyFluence[index] += history.WT * tmp / w
        circleEnergyFluence[:] /= self.NINCPHSP * circleArea
        return circleRadius, circleEnergyFluence

    def deriveDirectionFluence(self, minCord, maxCord, minRef, maxRef, N, histories=None, direction='X'):
        linearInter = lambda x1, y1, x2, y2, y: (y - y1) / (y2 - y1) * (x2 - x1) + x1
        if direction == 'X':
            D, Y, fluence = self.deriveFluence(minCord, maxCord, minRef, maxRef, N, histories)
            D_Fluence = np.add.reduce(fluence, axis=1)
        if direction == 'Y':
            X, D, fluence = self.deriveFluence(minRef, maxRef, minCord, maxCord, N, histories)
            D_Fluence = np.add.reduce(fluence, axis=0)

        averageMax = np.average(D_Fluence[(minRef < D) * (D < maxRef)])
        indexTemp = (np.arange(N))[D_Fluence > averageMax / 2]
        indexTemp1, indexTemp2 = indexTemp[0], indexTemp[-1]
        FWHM_Left = linearInter(D[indexTemp1 - 1], D_Fluence[indexTemp1 - 1], D[indexTemp1], D_Fluence[indexTemp1],
                                averageMax / 2)
        FWHM_Right = linearInter(D[indexTemp2 - 1], D_Fluence[indexTemp2 - 1], D[indexTemp2], D_Fluence[indexTemp2],
                                 averageMax / 2)
        plt.plot(D, D_Fluence)
        plt.plot([FWHM_Left, FWHM_Right], [averageMax / 2, averageMax / 2], '--')
        plt.xlabel(('%s/cm' % direction).lower())
        plt.ylabel('histories')
        plt.title('Fluence distribution on direction %s' % direction)
        plt.text(0, averageMax / 5, 'FWHM of %s = %g' % (self.phspFileName, FWHM_Right - FWHM_Left))
        return D, D_Fluence

    # def deriveCrossSectionFluence(self, minX, maxX, minY, maxY, N, histories=None, historyType='all'):
    #     def create_window():
    #         win = ImageDialog(edit=False, toolbar=True, wintitle="Cross Sections Fluence of %s" % self.phspFileName,
    #                           options=dict(show_xsection=True, show_ysection=True))
    #         win.resize(600, 600)
    #         return win
    #
    #
    #     X, Y, fluence = self.deriveFluence(minX, maxX, minY, maxY, N, histories, historyType)
    #     import guidata
    #     _app = guidata.qapplication()
    #     win = create_window()
    #     image = make.image(data=fluence, xdata=[minX, maxX], ydata=[minY, maxY])
    #     plot = win.get_plot()
    #     plot.add_item(image)
    #     win.exec_()

    def derive3D_Fluence(self, minX, maxX, minY, maxY, N, histories=None):
        X, Y, fluence = self.deriveFluence(minX, maxX, minY, maxY, N, histories)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(X, Y)
        ax.plot_surface(X, Y, fluence, rstride=8, cstride=8, alpha=0.3)
        ax.contour(X, Y, fluence, zdir='z', offset=-100, cmap=cm.coolwarm)
        ax.contour(X, Y, fluence, zdir='x', offset=minX, cmap=cm.coolwarm)
        ax.contour(X, Y, fluence, zdir='y', offset=maxY, cmap=cm.coolwarm)
        plt.show()

    def listHistories(self, num=100):
        self.reload()
        index = 0
        for history in self:
            print(history)
            index += 1
            if index == num:
                break
        self.close()

    def saveAs(self, new_phsp):
        '''
        把相空间当中的粒子经过modify操作后保存为new_phsp
        '''
        # PhotonIdx = np.where((self.latch >> 29) % (2**2) == 0)[0]
        # LowEIdx = np.where(abs(self.energy) <= 0.12)[0]
        # conditions = set(PhotonIdx).intersection(set(LowEIdx))
        conditions = np.abs(self.energy) <= 0.12
        particle_num = self.energy[conditions].shape[0]
        self.NPPHSP = particle_num
        self.NPHOTPHSP = particle_num
        with open(new_phsp, 'wb')  as f:
            f.write(self.MODE_RW.encode('utf-8'))
            phspHead = struct.pack('iifff', self.NPPHSP, self.NPHOTPHSP, self.EKMAXPHSP, self.EKMINPHSPE, self.NINCPHSP)
            f.write(phspHead)
            f.write(self.Nothing)
            index = 0
            for i,c in enumerate(conditions):
                if c:
                    binary = struct.pack('i6f',self.latch[i],self.energy[i],
                                  self.x[i],self.y[i],self.u[i],
                                  self.v[i],self.wt[i])
                    f.write(binary)





        #  ************************

        # self.MODE_RW = fid.read(5).decode('utf-8')
        # assert self.MODE_RW == 'MODE2' or self.MODE_RW == 'MODE0', self.MODE_RW
        # self.phspHead = fid.read(4 * 5)
        # self.Nothing = fid.read(3) if self.MODE_RW == 'MODE0' else fid.read(7)
        # self.NPPHSP, self.NPHOTPHSP, self.EKMAXPHSP, self.EKMINPHSPE, self.NINCPHSP = struct.unpack(headType,
        #                                                                                             self.phspHead)
        #
        # type = np.dtype(
        #     [('LATCH', 'i4'), ('Energy', 'f4'), ('X', 'f4'), ('Y', 'f4'), ('U', 'f4'), ('V', 'f4'), ('WT', 'f4')])
        # all_particle_infor = np.fromfile(self.phspFileName, dtype=type, count=-1, sep="", offset=28)
        # self.latch = all_particle_infor['LATCH']
        # self.energy = all_particle_infor['Energy']
        # self.x = all_particle_infor['X']
        # self.y = all_particle_infor['Y']
        # self.u = all_particle_infor['U']
        # self.v = all_particle_infor['V']
        # self.w = np.sqrt(1.0 - np.square(all_particle_infor['U']) - np.square(all_particle_infor['V']))
        # self.wt = all_particle_infor['WT']
        #
        # self.particle_type = np.round(self.latch / 1e9)  # 1：electron  0：photon
        # self.create_positoin = np.floor(self.latch / 16777216)











def binData2Hisotries(binData):
    num = len(binData) / 32
    assert len(binData) % 32 == 0
    for i in range(num):
        yield binData[i * 32:(i + 1) * 32]


def foldPhspSpace(phsp, N=None, phspFolded=None):
    phspFile = open(phsp, 'rb')
    if not phspFolded:
        phspFolded = phsp[:-9] + '_folded' + phsp[-9:]
    MODE_RW = phspFile.read(5)
    phspHead = phspFile.read(4 * 5)
    headType = 'iifff'
    NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, NINCPHSP = struct.unpack(headType, phspHead)
    tmp = phspFile.read(7)

    if not N:
        N = 512 * 1024 * 1024 / 32  # 512MB
    outLoop = NPPHSP / N + 1 if NPPHSP % N else NPPHSP / N

    phspFoldedFile = open(phspFolded, 'wb')

    phspFoldedFile.write(MODE_RW)
    phspFoldedFile.write(struct.pack(headType, NPPHSP * 4, NPHOTPHSP * 4, EKMAXPHSP, EKMINPHSPE, NINCPHSP * 4))
    phspFoldedFile.write(tmp)

    print('folding phase space...')
    for i in range(outLoop):
        binData = phspFile.read(4 * 8 * N)
        # print i
        for binHistory in binData2Hisotries(binData):
            binDataOut = ''
            history = History(binHistory)
            binDataOut += history.pack()
            history.X, history.U = -history.X, -history.U
            binDataOut += history.pack()
            history.Y, history.V = -history.Y, -history.V
            binDataOut += history.pack()
            history.X, history.U = -history.X, -history.U
            binDataOut += history.pack()
            # print 'history %d is folded.' % i
            phspFoldedFile.write(binDataOut)
    print('the size of folded phsp file is %d' % (5 + 4 * 5 + 7 + 4 * 8 * NPPHSP * 4))
    phspFile.close()
    phspFoldedFile.close()


def combinePhsp(phspFiles, combinePhspFile):
    print(phspFiles)
    combinePhsp = open(combinePhspFile, 'wb')
    MODE_RW = 'MODE2'
    NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, NINCPHSP = 0, 0, 0.0, 1.0e8, 0
    phaseSpaces = [PhaseSpace(phspFile) for phspFile in phspFiles]
    for phaseSpace in phaseSpaces:
        assert phaseSpace.MODE_RW == MODE_RW, phaseSpace.MODE_RW
        NPPHSP += phaseSpace.NPPHSP
        NPHOTPHSP += phaseSpace.NPHOTPHSP
        if phaseSpace.EKMAXPHSP > EKMAXPHSP:
            EKMAXPHSP = phaseSpace.EKMAXPHSP
        if phaseSpace.EKMINPHSPE < EKMINPHSPE:
            EKMINPHSPE = phaseSpace.EKMINPHSPE
        NINCPHSP += phaseSpace.NINCPHSP

    combinePhsp.write(MODE_RW)
    headType = 'iifff'
    combinePhsp.write(struct.pack(headType, NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, NINCPHSP))
    combinePhsp.write(phaseSpaces[0].NOTHING)
    for phaseSpace in phaseSpaces:
        cacheSize = 512 * 1024 * 1024 / 32  # 512MB
        outLoop = phaseSpace.NPPHSP / cacheSize + 1 if phaseSpace.NPPHSP % cacheSize else phaseSpace.NPPHSP / cacheSize

        print('Combining phase space: %s' % phaseSpace.phspFileName)

        for j in range(outLoop):
            tmp = phaseSpace.phspFile.read(cacheSize * 32)
            assert len(tmp) % 32 == 0, len(tmp)
            combinePhsp.write(tmp)
        phaseSpace.close()
    combinePhsp.close()


if __name__ == '__main__':
    # phsp_path = r"M:\Simulation\egsphantom\new\3-分析出模体后的相空间文件\noCT\noCT-Air-FFF-X27Y27-SSD100-PeakEnergy6.6.egsphsp1"
    phsp_path = r"/home/uih/New_4.0/FS27_P2.60_D60/total-ssd60.egsphsp1"
    p1 = PhaseSpace(phsp_path)
    p1.readphsp()
    new_phsp = r'/home/uih/New_4.0/FS27_P2.60_D60/LLLow-ssd60.egsphsp1'
    p1.saveAs(new_phsp)
    # p1.deriveAreaFluence(-20, 20, -20, 20, 100, p1.NPPHSP, 'photon')

    # path = r"M:\Simulation\egsphantom\new\3-分析出模体后的相空间文件\dosxyz里面的模体\withCT\Source21\IBL\1-Water\withCT-Water-Source21-SSD100-IBL-X-15_5Y-7_3.IAEAphsp"
    # p1 = PhaseSpace(path)
    # p1.readIAEAphsp()
    # p1.deriveAreaFluence(-20, 20, -20, 20, 100, p1.NPPHSP, 'photon')