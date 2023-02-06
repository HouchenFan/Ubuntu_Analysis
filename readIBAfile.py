import os
import numpy as  np
import matplotlib.pyplot as plt
from commonFunction import *

class IBA_Cure(object):  # 32102.703
    pass



def  toDic(datalist):
    dic = {}
    for ic in datalist:
        ic = ic.split()
        dic[ic[0].replace('%','')] = list(map(lambda x :x.strip(),ic[1:]))
    return  dic



def readIBAfile(dir_or_file):
    
    iba_file_list = []
    if os.path.isdir(dir_or_file):
        iba_file_list = findEveryFolder(dir_or_file,'.asc')

    else:
        iba_file_list.append(dir_or_file)


    data_list = []
    for file_path in iba_file_list:

        with open(file_path, 'rb') as f:
            content = f.readline()
            if b'\0' in content:
                print('文件格式是asc，但是内容不是文本文件，请检查.')
                return

        with open(file_path, 'r') as f:
    
            temp = f.readline().split()
            if temp:
                cureNum = int(temp[1])
                temp = [f.readline() for i in range(3)]
    
                for i in range(cureNum):
                    curve_obj = IBA_Cure()
                    curve_obj.file_name = os.path.split(file_path)[-1]
                    temp = [f.readline() for i in range(25)]
                    curve_obj.heard_infor = toDic(temp)
    
                    curve_obj.Type = 'PDD' if curve_obj.heard_infor['SCN'][0] == 'DPT' else 'Profile'
                    curve_obj.Date = curve_obj.heard_infor['DAT']
                    curve_obj.Time = curve_obj.heard_infor['TIM']
                    curve_obj.FS = [np.round(float(i) , 1) for i in curve_obj.heard_infor['FSZ']]
                    curve_obj.Mode = curve_obj.heard_infor['BMT'][0]
                    curve_obj.Energy = curve_obj.heard_infor['BMT'][1]
                    curve_obj.SSD = np.round(float(curve_obj.heard_infor['SSD'][0]) , 1)
    
                    temp = [f.readline() for i in range(3)]  # DUMP
                    cureData = ''
                    while True:
                        temp = f.readline()
                        if temp[0] != '=':
                            break
                        else:
                            cureData += '[' + ','.join(temp.replace('=', '').split()) + '],'
    
                    temp = [f.readline() for i in range(2)]  # end
    
                    originalData = eval('[' + cureData[:-1] + ']')
                    originalData = np.array(originalData).reshape(-1, 4)
                    curve_obj.Data = np.zeros((len(originalData[:, 0]), 2))
                    if curve_obj.Type == 'PDD':
                        curve_obj.Depth = 0
                        curve_obj.Data[:, 0] = np.round(originalData[:, 2] , 2)
                        curve_obj.Data[:, 1] = np.round(originalData[:, 3] / np.max(originalData[:, 3]) * 100, 1)
                        curve_obj.Data = curve_obj.Data[np.argsort(curve_obj.Data[:, 0])]  # 对z坐标按照 从小到大排序
                    elif curve_obj.Type == 'Profile':
                        # 防止出现0.1 没有检测出来的问题，以及离轴扫描的profile
                        if np.abs(originalData[0, 0] - originalData[1, 0]) <= 0.1 and np.abs(
                                originalData[0, 1] - originalData[1, 1]) >= 0.09:  # crossline  x那一列为0，为crossline。软件bug
                            #                        第一列坐标都相同                                  第二列坐标逐渐变化
                            curve_obj.Type = 'Crossline'
                            curve_obj.Depth = np.round(originalData[0, 2], 2)  # 单位为cm
                            curve_obj.Data[:, 0] = np.round(originalData[:, 1] , 2)
                            curve_obj.Data[:, 1] = originalData[:, 3]
                        elif np.abs(originalData[0, 0] - originalData[1, 0]) >= 0.09 and np.abs(originalData[0, 1] - originalData[1, 1]) <= 0.1:
                            curve_obj.Type = 'Inline'
                            curve_obj.Depth = np.round(originalData[0, 2], 2)  # 单位为mm
                            curve_obj.Data[:, 0] = np.round(originalData[:, 0] , 2)
                            curve_obj.Data[:, 1] = originalData[:, 3]
                        elif np.abs(originalData[0, 0] - originalData[1, 0]) > 0.1 and np.abs(
                                originalData[0, 1] - originalData[1, 1]) > 0.1:
                            curve_obj.Type = 'Dia--/++' if originalData[0, 0] * originalData[0, 1] > 0 else 'Dia-+/+-'
                            curve_obj.Depth = np.round(originalData[0, 2] , 2)  # 单位为mm
                            x = originalData[:, 0]
                            y = originalData[:, 1]
                            offaxis = np.sqrt(x * x + y * y)
    
                            id = [i for i in range(len(x) - 1) if x[i] * x[i + 1] <= 0]  # 找到对角线x轴的 正负坐标的分界点
    
                            offaxis[0: id[0]] = - offaxis[0:id[0]]
                            curve_obj.Data[:, 0] = np.round(offaxis , 2)
                            curve_obj.Data[:, 1] = originalData[:, 3]
    
                    data_list.append(curve_obj)

    return data_list


if __name__ == '__main__':
    path =r"\\dataserver03\rt\06_PH\107-装机医院\110033-北京协和医院\01-水箱数据\协和asc\FF\2-Profile&PDD测量-MLC&Jaw开野-CC13探头"

    all_cures = readIBAfile(path)

    conditions = {'FS':[[100,100],[200,200]],
                 'Type': ['PDD'],
                 'Depth': [0],
                  }

    select_cures = selectCure(all_cures, conditions)

    plt.ion()
    fig,ax = plt.subplots()
    for cure in select_cures:
        ax.plot(cure.Data[:, 0], cure.Data[:, 1],label=cure.file_name)
        ax.set_xlabel('Off-axis (mm)')
        ax.set_ylabel('Relative Dose ')
        ax.set_title('IBA curves')

    ax.legend(loc=0)
    plt.ioff()
    plt.show()




