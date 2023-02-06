import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('TkAgg')
from  commonFunction import *

class Cal_Cure:
    pass


def findSortedPosition(theList, target):

    list  = np.array(theList)
    return np.argmin(abs(list-target))

def read3ddose(dir_or_file,point=[0,0,3.3]): #3.3 #0.919   #1.476
    # dir_or_file: 3ddose文件路径 或者文件夹路径
    # point  ： 提取指定深度的profile,按照指定深度查找，如果没指定，默认1.5cm，假如3ddose只有5cm深度的profile，依然能够定位到5cm的深度
    #z_phsp = 0   # 相空间平面的Z坐标  by houchen.fan 2022/8/25
    cal_file_list = []
    if os.path.isdir(dir_or_file):
        cal_file_list = findEveryFolder(dir_or_file, '.3ddose')

    else:
        cal_file_list.append(dir_or_file)

    all_cure_list  = []

    for cal_file in cal_file_list:
        data = Cal_Cure()
        with open(cal_file,'r') as f:

            data.file_name = os.path.split(cal_file)[-1]
            temp = f.readline().split()
            nx = int(temp[0])
            ny = int(temp[1])
            nz = int(temp[2])

            bx = list(map(float, f.readline().split()))
            by = list(map(float, f.readline().split()))
            bz = list(map(float, f.readline().split()))

            dose1D = list(map(float, f.readline().split()))
            # dose1D = list(map(float, f.readline().split()))
        data.dose1D = np.array(dose1D)
        data.dose2D = np.reshape(dose1D, (nx * ny, nz))
        data.dose3D = np.reshape(dose1D, (nz, ny, nx))

        data.cx =  [(bx[i] + bx[i + 1]) / 2 for i in range(len(bx) - 1)]
        data.cy =  [(by[i] + by[i + 1]) / 2 for i in range(len(by) - 1)]
        data.cz =  [(bz[i] + bz[i + 1]) / 2 for i in range(len(bz) - 1)] 

        # # 提取pdd数据，默认提取中心轴的数据
        raw = findSortedPosition(data.cy, 0 + point[1])  # 在cy中找位置为 0 的索引
        column = findSortedPosition(data.cx, 0 + point[0])
        page = findSortedPosition(data.cz, bz[0] + point[2])  # 若用户传入点信息，pdd按照该点进行归一
        if nx > 1 and ny > 1 and nz > 1:
            data.detectorplane = data.dose3D[page]
            save_file = os.path.join(os.path.split(cal_file)[0], 'EPID.npy')     #保存EPID模体中闪烁体层的剂量数据
            np.save(save_file, data.detectorplane)

        pdd = [data.dose3D[page][raw][column] for page in range(len(data.cz))]
        data.DD = pdd
        data.PDD = pdd / max(pdd) if point[2] == 0 else pdd / pdd[findSortedPosition(data.cz, point[2])]

        # 提取指定深度的profile,按照指定深度查找，如果没指定，默认1.5cm，假如3ddose只有5cm深度的profile，依然能够定位到5cm的深度
        index = findSortedPosition(data.cz, 1.5) if point[2] == 0 else page
        oar = data.dose3D[index, raw, :]
        axisDose = oar[findSortedPosition(data.cx, 0)]
        data.profileX = oar
        data.profileXnorm = oar / axisDose


        oar = data.dose3D[index, :, column]
        data.profileY = oar
        data.profileYnorm = oar / axisDose

        #data.cz = (np.array(data.cz) - z_phsp).tolist()    # 按照相空间平面Z坐标归一  by houchen.fan 2022/8/25

        data.cx = np.round(np.array(data.cx) * 10, 2)
        data.cy = np.round(np.array(data.cy) * 10, 2)
        data.cz = np.round(np.array(data.cz) * 10, 2)
        if nx == 1 and ny == 1 and nz > 1:
            save_file_DD = os.path.join(os.path.split(cal_file)[0], 'DD.npy')
            DD = []
            DD.append(data.cz)
            DD.append(data.DD)
            np.save(save_file_DD, DD)     #保存水箱DD数据
            save_file_PDD = os.path.join(os.path.split(cal_file)[0], 'PDD.npy')
            PDD = []
            PDD.append(data.cz)
            PDD.append(data.PDD)
            np.save(save_file_PDD, PDD)
        elif nx > 1 and ny == 1 and nz == 1:
            save_file_X = os.path.join(os.path.split(cal_file)[0], 'profileX.npy')
            PFX = []
            PFX.append(data.cx)
            PFX.append(data.profileX)
            np.save(save_file_X, PFX)     #保存水箱profileX数据
            save_file_Xnorm = os.path.join(os.path.split(cal_file)[0], 'profileXnorm.npy')
            PFXnorm = []
            PFXnorm.append(data.cx)
            PFXnorm.append(data.profileXnorm)
            np.save(save_file_Xnorm, PFXnorm)
        data.Type = 'Crossline'
        data.Depth = 100
        data.FS = [100,100]
        data.Data = np.ones((len(data.cx),2))
        data.Data[:,0]=data.cx
        data.Data[:,1]=data.profileXnorm
        all_cure_list.append(data)

    return all_cure_list




def plotData(data_list):
    
    plt.figure()
    for data in data_list:
        if len(data.PDD)>2:
            plt.plot(data.cz, data.PDD,label=data.file_name)
            plt.title('PDD')
            plt.xlabel('off-axis (mm)')
            plt.ylabel('pdd')
            plt.legend(loc='best')


    plt.figure()
    for data in data_list:
        if len(data.profileXnorm) > 2 :
            linestyle = '-'
            if 'Double' in data.file_name:
                linestyle = '--'
                # data.cx = data.cx - 0.6
            plt.plot(data.cx,data.profileXnorm,linestyle=linestyle, label=data.file_name)

            plt.title('Crossline')
            plt.xlabel('off-axis (mm)')
            plt.ylabel('oar')
            plt.legend(loc='best')

    plt.figure()
    for data in data_list:
        if len(data.profileYnorm) > 2:
            plt.plot(data.cy, data.profileY, label=data.file_name)
            plt.title('Inline')
            plt.xlabel('off-axis (mm)')
            plt.ylabel('oar')
            plt.legend(loc='best')
    #
    # plt.ioff()
    plt.show()




if __name__ == '__main__':
    # path =r"M:\10.8.61.47\3ddose\X5Y5"
    path = r"/home/uih/Water/EPID_Air_FS27_P2.60/EGS.3ddose"
    # wt_list = [0, 1, 5, 10, 15, 20, 25, 30]
    # for wt in wt_list:
    #     path_wt = path + str(wt)
    cure_list = read3ddose(path, [0, 0, 0.919]) #0.919   #1.476  #3.376
    plotData(cure_list)

    # SW = [0,1,5,10,15,20,25,30]
    # DD = np.zeros((3, 8), dtype=float)
    # PDD = DD
    # filename = ["EGS","EGS_FBCTPhantom","EGS_IdeaPhantom"]
    # for i in range(len(SW)):
    #     path =r"D:\Data_fhc\EGSnrc\EPID_DOS_FS10_SW"+str(SW[i])
    #     cure_list =  read3ddose(path)
    #
    # #plotData(cure_list)
    #
    #     for j in range(len(cure_list)):
    #         DD[j, i] = cure_list[j].DD[9]
    # for i in range(len(DD)):
    #     PDD[i] = DD[i] / DD[i,3]
    # plt.figure()
    # for i in range(len(PDD)):
    #     plt.plot(SW, PDD[i,:],label=str(filename[i]))
    #     plt.title('PDD')
    #     plt.xlabel('SolidWater (cm)')
    #     plt.ylabel('pdd')
    #     plt.legend(loc='best')
    # plt.show()




    #print(np.interp(3.3,cure_list[0].cz,cure_list[0].PDD))
    # method = 'varian'    #  inflection  or  varian
    # beam_mode = 'FFF'
    #
    # if method == 'varian':
    #     cure_list =  renormalized(cure_list)
    #
    # for cure in cure_list:
    #     # plt.plot(cure.Data[:, 0], cure.Data[:, 1])
    #     lr_20 = calFieldSize(beam_mode, cure.Data[:,0], cure.Data[:,1], method=method, field_size_percent=0.2)
    #     lr_80 = calFieldSize(beam_mode, cure.Data[:,0], cure.Data[:,1], method=method, field_size_percent=0.8)
    #
    #     pen_left = np.round( lr_80[1] - lr_20[1] ,2)
    #     pen_right = np.round( lr_20[2] - lr_80[2] ,2 )
    #     print(f'Fieldsize {cure.file_name} left penubra is {pen_left} mm ------ right penubra is {pen_right} mm  ----- average:{np.round((pen_right+pen_left)/2,2)}mm')

    print('Done')


