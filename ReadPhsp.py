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
#import sys
import os
import scipy.interpolate as spi
import binascii
import multiprocessing  #add by FHC 2022/9/9
#import threading  #add by FHC 2022/9/9

global EKMAXPHSP  # 所有粒子的最大动能 add by FHC 2022/9/7
EKMAXPHSP = 0
global EKMINPHSPE  # 所有电子的最小动能 add by FHC 2022/9/7
EKMINPHSPE = 0.5
global SumElectron  # 电子数量 add by FHC 2022/9/9
SumElectron = 0
global NPHOTPHSP  # 光子数量  add by FHC 2022/9/9
NPHOTPHSP = 0
global NPPHSP  # 所有粒子的数量
NPPHSP = 0

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
    def __init__(self, PID, LATCH, Energy, X, Y, U, V,W_temp, WT, particle_type, ZLast=0,):
        self.LATCH = LATCH
        self.PID = PID  #粒子在相空间文件中的序号
        self.Bremsstrahlung = Cut(self.LATCH, 0, 1)      #粒子轫致辐射的信息
        self.InteractiveRegion = Cut(self.LATCH, 1, 23)  #粒子在相互作用区域的信息
        self.SecondaryParticle = Cut(self.LATCH, 24, 5)  #次级粒子的信息
        self.Charge = Cut(self.LATCH, 29, 2)             #粒子电荷的信息
        self.MultiScore = Cut(self.LATCH, 31, 1)         #多个score平面的信息
        self.ZLast = ZLast                               #粒子最后相互作用的Z坐标
        self.Energy = Energy                             #粒子的能量
        self.EkElectron = 0.5
        self.particle_type = particle_type
        global NPHOTPHSP
        #global SumElectron
        global EKMAXPHSP
        global EKMINPHSPE
        if self.particle_type == b'\x01':  #判断是否为光子
            self.EkAll = abs(Energy)
            NPHOTPHSP += 1
            self.Charge = 0
        elif self.particle_type == b'\x02':  #判断是否为电子
            self.EkAll = abs(Energy)
            self.EkElectron = self.EkAll
            self.Energy = abs(Energy) + 0.511  #IAEA输出的能量为动能
            self.Charge = 2
            #SumElectron += 1
        else:                       #判断是否为正电子
            self.EkAll = abs(Energy)
            self.Energy = abs(Energy) + 0.511   #IAEA输出的能量为动能
            self.Charge = 1
        self.LATCH += self.Charge << 29

        if self.EkAll > EKMAXPHSP:         #更新所有粒子的最大动能
            EKMAXPHSP = self.EkAll
        if self.EkElectron < EKMINPHSPE:   #更新所有电子的最小动能
            EKMINPHSPE = self.EkElectron

        self.X = X
        self.Y = Y
        self.U = U

        if not -1 <= self.U <= 1:
            print('NO.{:8d}:U is wrong:{:8.3f}'.format(self.PID, self.U))
            raise AssertionError
        self.V = V
        if not -1 <= self.V <= 1:
            print('NO.{:8d}:V is wrong:{:8.3f}'.format(self.PID, self.V))
            raise AssertionError
        self.W = math.sqrt(W_temp)
        # if self.particle_type == b'\xff':  #判断W是否要带负号 add by FHC 2022/9/2
        #     self.W = self.W * -1
        self.WT = WT
        if not 0 <= self.WT <= 1:
            print('NO.{:8d}:WT is wrong:{:8.3f}'.format(self.PID, self.WT))
            raise AssertionError
        #self.R = math.sqrt(X ** 2 + Y ** 2)

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
        self.FileName = filename  #IAEA相空间文件名
        self.HeaderName = headername  #IAEAheader文件名
        self.NINCPHSP = 0  #History数量
        self.lock = multiprocessing.RLock()
        #读取相空间文件基本信息  add by FHC 2022/9/7
        with open(filename, 'rb') as fid, \
             open(headername, 'r', encoding='utf-8') as fid_header:
            while (self.NINCPHSP == 0):
                lines = fid_header.readline()  #整行读取数
                if ("$ORIG_HISTORIES:") in lines:
                    self.NINCPHSP = int(next(fid_header).replace("\n"," "))     #从IAEAheader读取History 数量
            fid_header.close()

            self.Mode = 'MODE2'   #IAEAphsp文件为默认MODE2
            self.TotalNumParticles = int(os.path.getsize(filename) / 33)
            if self.Mode == 'MODE2':
                self.offset = 33
            elif self.Mode == 'MODE0':
                self.offset = 28
        self.LoopNum = multiprocessing.cpu_count()   #获取进程数量 = 逻辑CPU数量
        self.MAXBuffer = math.ceil(self.TotalNumParticles / self.LoopNum)  #1000000  #单个进程内读取的粒子数量
        self.startTime = datetime.now()
        self.stopTime = datetime.now()
        self.NPHOTPHSP = 0



    def Loop(self, startID=1, stopID=-1, ptype='ope+-'):  #单个进程
        ptypebin = []   #可能存在的粒子种类
        if 'p' in ptype:  #光子
            ptypebin.append(0)
        if 'e+-' in ptype:  #正负电子
            ptypebin.append(2)
            ptypebin.append(1)
        elif 'e-' in ptype:  #电子
            ptypebin.append(2)
        elif 'e+' in ptype:  #正电子
            ptypebin.append(1)
        if 'o' in ptype:  #未知
            ptypebin.append(3)
        if len(ptypebin) == 0:
            print("No matched particle type.")
            return 0

        if stopID == -1:
            stopID = self.MAXBuffer
        phspList = []  #创建单个进程中读取的粒子相空间信息
        IDoffset = (startID - 1) * self.offset  #获取当前相空间文件的初始读取位置
        with open(self.FileName, 'rb') as fid:   #以二进制只读的方式，打开相空间文件
            fid.seek(IDoffset, 0)  #fid.seek(self.offset + IDoffset, 0)    将指针移动到当前进程应该在的初始读取位置
            if startID > self.TotalNumParticles:
                print("Start ID excessed the max particle number.")
                return []
            #typelist = []
            for i in range(min(self.TotalNumParticles - startID + 1, self.MAXBuffer, stopID - startID + 1)):  #对进程内部的粒子进行循环读取
                particle_type = struct.unpack('s', fid.read(1))[0]
                #typelist.append(particle_type)#粒子的type，代表Z方向速度的正负号
                Energy, X, Y, Z, U, V, WT = struct.unpack('7f', fid.read(28))  #粒子的能量、坐标（X, Y, Z）、水平速度（U, V）、权重（WT）
                LATCH = struct.unpack('i', fid.read(4))[0]  #粒子的LATCH信息

                W_temp = (1 - U ** 2 - V ** 2)  #获取Z方向速度的平方值，由于下一步的筛选，以免出现对负值开方的情况
                if particle_type != b'\xff' and abs(W_temp) > 1e-07:   #只选择Z方向速度为正 并且 Z速度平方为正的粒子进行后处理
                    p = PhspVector(i + startID, LATCH, Energy, X, Y, U, V, W_temp, WT, particle_type, Z)   # 对粒子的相空间信息进行后处理
                    # if p.Charge in ptypebin and p.SecondaryParticle in In_SecondaryParticle and p.InteractiveRegion in In_InteractiveRegion:
                    print("Particle Type: ", particle_type)
                    p.Show()
                    #if p.Charge in ptypebin:
                    phspList.append(p)  #添加粒子的相空间信息

                # if i+startID % 100000 == 1 or i+startID == self.TotalNumParticles:
                # self.stopTime = datetime.now()
                #     percent = (i + startID) % self.MAXBuffer / self.MAXBuffer
                #     timedelta = (self.stopTime-self.startTime).seconds
                #     print("Proceeding {:8.2%}, time used:{:4d} seconds, rest time:{:6.1f} seconds".format(percent, timedelta, (1-percent)*timedelta/percent), end='\r', flush=True)
                #     if i+startID == self.TotalNumParticles:
                #         print("\nDone.")
        return phspList

    def Show(self, phspList, pNum=0, printout=True):
        if len(phspList) > 0:
            phspList[0].ShowTitle()
        for i, p in enumerate(phspList):
            if printout:
                p.Show()
        return pNum

    def ProcessLoop(self,MAXBuffer, ptype, i):
        # 并行模块 开始
        #self.lock.acquire()  # 请求并行锁
        self.phspList = self.Loop(i * MAXBuffer + 1, (i + 1) * MAXBuffer, ptype)
        self.phspList = self.TranserEPID(self.phspList)
        self.SavePhspProcess(self.phspList, i)  #将当前进程的粒子相空间信息保存为二进制文件
        #self.lock.release()  # 释放并行锁
        # 并行模块 结束

    def TotalLoop(self, ptype='ope+-'):
        #LoopNum = math.ceil(self.TotalNumParticles / self.MAXBuffer)
        # self.EnergyFluence.reset()
        self.startTime = datetime.now()
        #targetevents=0

        ProcessList = []
        # #并行模块 开始
        # threading.Lock().acquire()  #请求并行锁

        for i in range(self.LoopNum):
            #phspList = self.Loop(i * self.MAXBuffer + 1, (i + 1) * self.MAXBuffer, ptype)
            exec('Process_'+ str(i)+ ' = multiprocessing.Process(target = self.ProcessLoop, args=(self.MAXBuffer, ptype, i))')
            eval('ProcessList.append(Process_'+str(i)+')')

        for t in ProcessList:
            #t.daemon(True)
            t.start()
        t.join()

        self.SavePhspAll()

            #targetevents = targetevents + len(phspList)
            #pNum = self.show(phspList, ptype, pNum, i, printout=False)
        #
        # threading.Lock().release()  #释放并行锁
        # #并行模块 结束
        self.stopTime = datetime.now()
        print("Total time: {:d} seconds".format((self.stopTime-self.startTime).seconds))


    def TranserEPID(self, phspList):
        Z_epid = 45.0
        for i in range(len(phspList)):
            phspList[i].X += (phspList[i].U / phspList[i].W) * (Z_epid - phspList[i].ZLast)
            phspList[i].Y += (phspList[i].V / phspList[i].W) * (Z_epid - phspList[i].ZLast) + 15.0
            phspList[i].ZLast = 0
        return phspList

    def SavePhspAll(self):    #将所有进程的粒子相空间信息进行整合
        self.FileName = self.FileName.split('.')[0] + ".egsphsp1"
        self.headtype = 'iifff'   #相空间文件Header 的数据类型
        #self.datatype = 'iffffff'  #相空间文件粒子的数据类型
        self.Nothing = b'\x00\x00\x00'  #Head信息与粒子信息之间的间隔
        #NPPHSP = len(phspList)  #获取粒子总数量
        global NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE
        with open(self.FileName, 'wb') as fout:   #打开总相空间文件
            fout.write(b'MODE0')   #写入MODE信息
            phspHead = struct.pack(self.headtype, NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, self.NINCPHSP)
            fout.write(phspHead)  #写入Head信息
            fout.write(self.Nothing)  #写入Head信息与粒子信息之间的间隔
            for i in range(self.LoopNum):
                self.FinName = self.FileName.split('.')[0] + "_" + str(i) + ".egsphsp1"
                with open(self.FinName, 'rb') as fin:  #打开每一个进程的粒子相空间信息
                    singleHead = fin.read(4*4)
                    npphsp, nphotphsp, ekmaxphsp, ekminphspe = struct.unpack('iiff', singleHead)
                    NPPHSP += npphsp
                    NPHOTPHSP += nphotphsp
                    EKMAXPHSP = max(EKMAXPHSP, ekmaxphsp)
                    EKMINPHSPE = min(EKMINPHSPE, ekminphspe)
                    data = fin.read()
                    fout.write(data)
                fin.close()
        fout.close()
        with open(self.FileName, 'rb+') as fout:
            fout.seek(5)
            phspHead = struct.pack(self.headtype, NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE, self.NINCPHSP)
            fout.write(phspHead)  # 写入Head信息
        fout.close()

    def SavePhspProcess(self,phspList, i):   #将当前进程的粒子相空间信息保存为二进制文件
        self.FileName = self.FileName.split('.')[0] + "_" + str(i) + ".egsphsp1"
        self.datatype = 'iffffff'  # 相空间文件粒子的数据类型
        self.headtype = 'iiff'  # 相空间文件Header 的数据类型
        global NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE
        NPPHSP = len(phspList)  # 获取粒子总数量
        with open(self.FileName, 'wb') as fout:     #打开每一个进程的粒子相空间信息
            phspHead = struct.pack(self.headtype, NPPHSP, NPHOTPHSP, EKMAXPHSP, EKMINPHSPE)
            fout.write(phspHead)  # 写入Head信息
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
        path_header = "/home/uih/Head/Source21/IBL_FS27_UpHead_PE2.6/Source21_Phant_FS27.IAEAheader"
        ph = Phsp(path_phsp, path_header)
        ph.TotalLoop()
        # thread_1 = threading.Thread(target = ph.TotalLoop())
        # thread_2 = threading.Thread(target = ph.TotalLoop())
        #ph.TranserEPID(phspList)
        #ph.SavePhsp(phspList)
