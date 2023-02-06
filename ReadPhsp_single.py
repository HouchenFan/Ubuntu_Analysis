# -*- coding: utf-8 -*-
"""
Created on Thu Jan 05 15:41:52 2020

@author: fuwei.zhao
"""

import struct
import math
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import sys
import os
import scipy.interpolate as spi
import binascii


def Cut(data, start, width):
    if not isinstance(data, int):
        print("Input is not int.")
        return
    cut_data = data >> start
    cut_data = cut_data % 2**width
    return cut_data


def Int2Binary(data, width):
    BASE=2
    quotient = data
    binstr=['0']*width
    for i in range(width):
        remainder = quotient%BASE
        quotient = quotient>>1
        binstr[width-i-1]=str(remainder)
    binstr=''.join(binstr)
    return binstr


class PhspVector:
    def __init__(self, PID, LATCH, Energy, X, Y, U, V, WT, particle_type, ZLast=0,):
        self.LATCH = LATCH
        self.PID = PID
        self.Bremsstrahlung = Cut(self.LATCH, 0, 1)
        self.InteractiveRegion = Cut(self.LATCH, 1, 23)
        self.SecondaryParticle = Cut(self.LATCH, 24, 5)
        self.Charge = Cut(self.LATCH, 29, 2)
        self.MultiScore = Cut(self.LATCH, 31, 1)
        self.ZLast = ZLast
        self.Energy = Energy
        self.EkElectron = 0.5
        global NPHOTPHSP
        if self.Charge == 0:  #判断是否为光子
            self.EkAll = abs(Energy)
            NPHOTPHSP += 1
        elif self.Charge == -1:  #判断是否为电子
            self.EkAll = abs(Energy) - 0.511
            self.EkElectron = self.EkAll
            SumElectron += 1
        else:                       #判断是否为正电子
            self.EkAll = abs(Energy) - 0.511
        global EKMAXPHSP
        global EKMINPHSPE
        if self.EkAll > EKMAXPHSP:  #更新所有粒子的最大动能
            EKMAXPHSP = self.EkAll
        if self.EkElectron < EKMINPHSPE:   #更新所有电子的最小动能
            EKMINPHSPE = self.EkElectron
        self.X = X
        self.Y = Y
        self.U = U
        self.particle_type = particle_type
        if not -1 <= self.U <= 1:
            print('NO.{:8d}:U is wrong:{:8.3f}'.format(self.PID, self.U))
            raise AssertionError
        self.V = V
        if not -1 <= self.V <= 1:
            print('NO.{:8d}:V is wrong:{:8.3f}'.format(self.PID, self.V))
            raise AssertionError
        W_temp = (1 - U ** 2 - V ** 2)
        self.W = math.sqrt(W_temp)
        # if self.particle_type == b'\xff':  #判断W是否要带负号 add by FHC 2022/9/2
        #     self.W = self.W * -1
        self.WT = WT
        if not 0 <= self.WT <= 1:
            print('NO.{:8d}:WT is wrong:{:8.3f}'.format(self.PID, self.WT))
            raise AssertionError
        self.R = math.sqrt(X ** 2 + Y ** 2)

    def ShowTitle(self):
        print(
            "{:^12} |{:^1}|{:^2}|{:^5}|{:^23}|{:^1}"
            "|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^10}|".format("ParticleID", "P", "Q", "SecP",
                                                                  "Interactive Region", "B",
                                                                  "X", "Y", "U", "V", "W", "Energy", "Weight"))

    def Show(self):
        print(
            "{:>12} |{:^1}|{:^2}|{:^5}|{:^23}|{:^1}|{:>8.3f}|{:>8.3f}|{:>8.3f}|{:>8.3f}|{:>8.3f}|{:>8.3f}|{:>10.3e}|".format(self.PID, Int2Binary(self.MultiScore, 1), Int2Binary(self.Charge, 2), Int2Binary(self.SecondaryParticle, 5), Int2Binary(self.InteractiveRegion, 23), Int2Binary(self.Bremsstrahlung, 1), self.X, self.Y, self.U, self.V, self.W, self.Energy, self.WT))


class Phsp:
    def __init__(self, filename,headername):
        self.FileName = filename
        self.HeaderName = headername
        #global NPHOTPHSP
        global NINCPHSP
        global NPPHSP
        #读取相空间文件基本信息  add by FHC 2022/9/7
        with open(filename, 'rb') as fid, \
             open(headername, 'r', encoding='utf-8') as fid_header:
            NINCPHSP = 0
            while (NINCPHSP == 0):
                lines = fid_header.readline()  #整行读取数
                if ("$ORIG_HISTORIES:") in lines:
                    NINCPHSP = int(next(fid_header).replace("\n"," "))     #history 数量
                # elif ("$PARTICLES:") in lines:
                #     NPPHSP = int(next(fid_header).replace("\n", " "))      #所有的粒子数量
                # elif ("$PHOTONS:") in lines:
                #     NPHOTPHSP = int(next(fid_header).replace("\n", " "))       #所有的光子数量
            fid_header.close()

            # self.Mode = str(struct.unpack('5s', fid.read(5))[0])[2:-1]
            # self.TotalNumParticles = struct.unpack('I', fid.read(4))[0]
            # self.PhotonNumParticles = struct.unpack('I', fid.read(4))[0]
            # self.MaxKineticEnergy = struct.unpack('f', fid.read(4))[0]
            # self.MinKineticEnergy = struct.unpack('f', fid.read(4))[0]
            # self.NumIncidentElectron = struct.unpack('f', fid.read(4))[0]
            self.Mode = 'MODE2'
            self.TotalNumParticles = int(os.path.getsize(filename) / 33)
            if self.Mode == 'MODE2':
                # temp = fid.read(7)
                self.offset = 33
            elif self.Mode == 'MODE0':
                # temp = fid.read(3)
                self.offset = 28
        self.MAXBuffer = 100#1000000
        # self.PriEnergyMap = EnergyMap()
        # self.SecEnergyMap = EnergyMap()
        # self.ThiEnergyMap = EnergyMap()
        self.startTime = datetime.now()
        self.stopTime = datetime.now()
        #fid.close()    # add by FHC 2022/9/7
        # 读取IAEAheader文件基本信息   add by FHC 2022/9/7



    def Loop(self, startID=1, stopID=-1, ptype='ope+-'):
        ptypebin = []
        if 'p' in ptype:
            ptypebin.append(0)
        if 'e+-' in ptype:
            ptypebin.append(2)
            ptypebin.append(1)
        elif 'e-' in ptype:
            ptypebin.append(2)
        elif 'e+' in ptype:
            ptypebin.append(1)
        if 'o' in ptype:
            ptypebin.append(3)
        if len(ptypebin) == 0:
            print("No matched particle type.")
            return 0
        In_InteractiveRegion=[int('00000000000000000000001', 2), int('00000000000001100000001', 2), int('00000000000001000000001', 2)]
        In_SecondaryParticle=[int('00001', 2)]
        In_MultiScore=[int('0', 2)]
        In_Bremsstrahlung=[int('1', 2)]
        if stopID == -1:
            stopID = self.MAXBuffer
        phspList = []
        IDoffset = (startID - 1) * self.offset
        with open(self.FileName, 'rb') as fid:
            fid.seek(IDoffset, 0)  #fid.seek(self.offset + IDoffset, 0)
            if startID > self.TotalNumParticles:
                print("Start ID excessed the max particle number.")
                return []
            global EKMAXPHSP  # add by FHC 2022/9/7
            EKMAXPHSP = 0
            global EKMINPHSPE  #add by FHC 2022/9/7
            EKMINPHSPE = 0.5
            global SumElectron
            SumElectron = 0
            for i in range(min(self.TotalNumParticles - startID + 1, self.MAXBuffer, stopID - startID + 1)):
                particle_type = struct.unpack('s', fid.read(1))[0]
                Energy, X, Y, Z, U, V, WT = struct.unpack('7f', fid.read(28))
                LATCH = struct.unpack('i', fid.read(4))[0]
                # if self.Mode =='MODE0':
                #     p = PhspVector(i + startID, LATCH, Energy, X, Y, U, V, WT)
                # elif self.Mode =='MODE2':
                    # ZLast=struct.unpack('f', fid.read(4))[0]
                W_temp = (1 - U ** 2 - V ** 2)
                if particle_type == b'\x01' and abs(W_temp) > 1e-07:
                    p = PhspVector(i + startID, LATCH, Energy, X, Y, U, V, WT, particle_type, Z)   # particle_type add by FHC 2022/9/2
                # if p.Charge in ptypebin and p.SecondaryParticle in In_SecondaryParticle and p.InteractiveRegion in In_InteractiveRegion:
                    print("Particle Type: ", particle_type)
                    p.Show()
                    if (p.Charge in ptypebin) and (p.W > 0):
                        phspList.append(p)
                    # if p.InteractiveRegion == In_InteractiveRegion[0]:
                    #     self.PriEnergyMap.Deposite(p)
                    # else:
                    #     self.SecEnergyMap.Deposite(p)
                    #
                    # self.ThiEnergyMap.Deposite(p)
                    # elif p.InteractiveRegion == int('00000000000001100000001', 2):
                    #     self.SecEnergyMap.Deposite(p)
                if i+startID % 100000 == 1 or i+startID == self.TotalNumParticles:
                    self.stopTime = datetime.now()
                    percent = (i + startID) / self.TotalNumParticles
                    timedelta = (self.stopTime-self.startTime).seconds
                    print("Proceeding {:8.2%}, time used:{:4d} seconds, rest time:{:6.1f} seconds".format(percent, timedelta, (1-percent)*timedelta/percent), end='\r', flush=True)
                    if i+startID == self.TotalNumParticles:
                        print("\nDone.")
        return phspList

    def Show(self, phspList, pNum=0, printout=True):
        if len(phspList) > 0:
            phspList[0].ShowTitle()
        for i, p in enumerate(phspList):
            if printout:
                p.Show()
        return pNum

    def TotalLoop(self, ptype='ope+-'):
        LoopNum = 1 #math.ceil(self.TotalNumParticles / self.MAXBuffer)
        # self.EnergyFluence.reset()
        self.startTime = datetime.now()
        targetevents=0
        global NPHOTPHSP
        NPHOTPHSP = 0
        global phspList
        phspList = []
        for i in range(LoopNum):
            #phspList = self.Loop(i * self.MAXBuffer + 1, (i + 1) * self.MAXBuffer, ptype)
            phspList.extend(self.Loop(i * self.MAXBuffer + 1, (i + 1) * self.MAXBuffer, ptype))
            targetevents = targetevents + len(phspList)
            #pNum = self.show(phspList, ptype, pNum, i, printout=False)
        # self.EnergyFluence.show()
        self.stopTime = datetime.now()
        # self.PriEnergyMap.Show()
        # self.SecEnergyMap.Show()
        # self.ThiEnergyMap.Show()
        print("Total target particles:{:8d}".format(targetevents))
        print("Total time: {:d} seconds".format((self.stopTime-self.startTime).seconds))
        return phspList

    def TranserEPID(self, phspList):
        Z_epid = 45.0
        for i in range(len(phspList)):
            phspList[i].X += (phspList[i].U / phspList[i].W) * (Z_epid - phspList[i].ZLast)
            phspList[i].Y += (phspList[i].V / phspList[i].W) * (Z_epid - phspList[i].ZLast) + 15.0
            phspList[i].ZLast = 0

    def SavePhsp(self,phspList):    #将相空间数据保存为二进制的egsphsp文件
        self.FileName = self.FileName.split('.')[0] + ".egsphsp1"
        self.headtype = 'iifff'   #相空间文件Header 的数据类型
        self.datatype = 'iffffff'  #相空间文件粒子的数据类型
        self.Nothing = b'\x00\x00\x00'
        NPPHSP = len(phspList)
        with open(self.FileName, 'wb') as fout:
            fout.write(b'MODE0')
            phspHead = struct.pack(self.headtype, NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, NINCPHSP)
            fout.write(phspHead)
            fout.write(self.Nothing)
            for i in range(len(phspList)):
                self.binData = struct.pack(self.datatype, phspList[i].LATCH, phspList[i].Energy, phspList[i].X, phspList[i].Y, phspList[i].U, phspList[i].V, phspList[i].WT)
                fout.write(self.binData)
        fout.close()

class EnergyMap:
    def __init__(self, XSize = 512, YSize = 512, XMax = 20.48, XMin = -20.48, YMax = 20.48, YMin = -20.48):
        self.XSize = XSize
        self.YSize = YSize
        self.XMax = XMax
        self.XMin = XMin
        self.YMax = YMax
        self.YMin = YMin
        self.XRes = (self.XMax - self.XMin) / self.XSize
        self.YRes = (self.YMax - self.YMin) / self.YSize
        self.Fluence = np.zeros([self.XSize, self.YSize])
        self.EnergyFluence = np.zeros([self.XSize, self.YSize])
        self.Response = np.zeros([self.XSize, self.YSize])
        self.ipo3 = 0
        self.Statistics = 0
        self.EnergyResponse(r"E:\工作文档\资料\EGS\DetectorResponse.txt")

    def Deposite(self, p):
        if not isinstance(p, PhspVector):
            print("Unexpected input type.")
            return
        if self.XMin < p.X < self.XMax and self.YMin < p.Y < self.YMax:
            xbin = math.ceil((p.X - self.XMin) / self.XRes) - 1
            ybin = math.ceil((p.Y - self.YMin) / self.YRes) - 1
            Response = self.InterPolate(p.Energy)
            self.Fluence[xbin, ybin] = self.Fluence[xbin, ybin] + p.WT
            self.EnergyFluence[xbin, ybin] = self.EnergyFluence[xbin, ybin] + p.WT*p.Energy
            self.Response[xbin, ybin] = self.Response[xbin, ybin] + p.WT*Response
            self.Statistics = self.Statistics+p.WT
            # self.Fluence[xbin, ybin] = self.Fluence[xbin, ybin] + 1
            # self.EnergyFluence[xbin, ybin] = self.EnergyFluence[xbin, ybin] + p.Energy
            # self.Response[xbin, ybin] = self.Response[xbin, ybin] + Response
            # self.Statistics = self.Statistics+1

    def EnergyResponse(self, path):
        data = np.loadtxt(path)
        self.ipo3 = spi.splrep(data[:, 0], data[:, 1], k=3)

    def InterPolate(self, x):
        y = spi.splev(x, self.ipo3)
        return y

    def Show(self):
        # plt.figure(figsize=(12, 9))
        plt.matshow(self.Fluence)
        plt.title(''.join(['Fluence Statictics: ', str(self.Statistics)]))
        plt.matshow(self.EnergyFluence)
        plt.title(''.join(['EnergyFluence Statictics: ', str(self.Statistics)]))
        plt.matshow(self.Response)
        plt.title(''.join(['Response Statictics: ', str(self.Statistics)]))

    def Save(self, path):
        with open(path+"_Fluence.dat", 'wb') as fout:
            for i in range(self.XSize):
                for j in range(self.YSize):
                    bytes=struct.pack('f', self.Fluence[i, j])
                    fout.write(bytes)
        with open(path+"_EnergyFluence.dat", 'wb') as eout:
            for i in range(self.XSize):
                for j in range(self.YSize):
                    bytes=struct.pack('f', self.EnergyFluence[i, j])
                    eout.write(bytes)
        with open(path+"_Response.dat", 'wb') as rout:
            for i in range(self.XSize):
                for j in range(self.YSize):
                    bytes=struct.pack('f', self.Response[i, j])
                    rout.write(bytes)




class EnergyFluence:
    def __init__(self, nbin, Rmax, binshape='circle'):
        self.nbin = nbin
        self.Rmax = Rmax
        self.Rstep = 0
        self.x = []
        self.y = []
        self.yfin = np.zeros(self.nbin)
        self.ynum = []
        self.reset()

    def reset(self):
        self.Rstep = self.Rmax / math.sqrt(nbin)
        self.x = [(math.sqrt(i + 1)) * self.Rstep for i in range(self.nbin)]
        self.y = [0 for i in range(self.nbin)]
        self.yfin = np.zeros(self.nbin)
        self.ynum = [0 for i in range(self.nbin)]

    def accumulate(self, x, y):
        if 0 < x < self.Rmax:
            xbin = math.ceil((x/self.Rstep)**2) - 1
            # xbin = math.ceil(x * self.nbin / self.Rmax) - 1
            self.y[xbin] = self.y[xbin] + y
            self.ynum[xbin] = self.ynum[xbin]+1

    def show(self):
        plt.figure(figsize=(12, 9))
        self.yfin = [self.y[i]/(math.pi * self.Rstep**2) for i in range(self.nbin)]
        plt.plot(self.x, self.yfin)
        plt.xlabel("R/cm")
        plt.ylabel("Energy Fluence")
        plt.title("Energy Fluence vs Position")
        plt.show()


if __name__ == '__main__':
        # print("args:")
        # for n, arg in enumerate(sys.argv):
        #     print("NO.{:2d} arg: {:}".format(n, arg))
        # del n, arg
        # if len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
        #     P = Phsp(sys.argv[1])
        #     P.TotalLoop()
        # elif not os.path.isfile(sys.argv[1]):
        #     print("File doesn't exist.")
        path_phsp = "/home/uih/Head/Source21/IBL_FS27_UpHead_PE2.6/Source21_Phant_FS27.egsphsp1"
        path_header = "/home/uih/Head/Vacuum_Gan0_FS27_Sigma0.07_D85/EGS.IAEAheader"
        ph = Phsp(path_phsp, path_header)
        ph.TotalLoop()
        ph.TranserEPID(phspList)
        ph.SavePhsp(phspList)
