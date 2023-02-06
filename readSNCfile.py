
# read ptw file(*.mcc)
# input  dir or full file name
# output: list  , each item  is obj.  the obj contain cures informationin .like ,depth,ssd,type,ect..

import os
import numpy as np
import matplotlib.pyplot as plt
from commonFunction import *

class SNC_Cure(object):  # 32102.703
    pass



def readSNCfile(dir_or_file):
    snc_file_list = []
    if os.path.isdir(dir_or_file):
        snc_file_list = findEveryFolder(dir_or_file,'.snctxt')

    else:
        snc_file_list.append(dir_or_file)

    all_cure_list  = []

    for snc_file in snc_file_list:

        with open(snc_file,'r') as f:
            all_str = f.read()
            all_cure_infor_list = strfindlr(all_str,'BEGIN SCAN','BEGIN RAW DATA')  # 得到每一条cure的 分段信息
            for cure_infor in all_cure_infor_list:
                temp_obj = SNC_Cure()

                # temp_obj.Type =  getCureType(cure_infor)
                X_FS = int(float(strfindlr(cure_infor, 'Summary FieldSize X (cm)\t', '\t\n')[0]))*10  #unit :mm
                Y_FS = int(float(strfindlr(cure_infor, 'Summary FieldSize Y (cm)\t', '\t\n')[0]))*10
                temp_obj.FS = [X_FS,Y_FS]
                temp_obj.Type = strfindlr(cure_infor, 'Scan Type	', '\t\n')[0]



                temp = strfindlr(cure_infor, 'Relative Dose (%)', 'END DOSE TABLE')[0]
                temp = re.sub(r'\s+', ' ', temp)
                temp = temp.split(' ')
                temp = [float(i)  for i in temp if i]
                Data  = np.array(temp).reshape(-1,4)

                if temp_obj.Type == 'Depth Scan':
                    temp_obj.Type == 'PDD'
                elif temp_obj.Type == 'Profile at Angle':

                    temp_obj.Type = 'Diagonal(++/--)' if Data[0,0] * Data[0,1] > 0 else 'Diagonal(-+/+-)'


                temp_obj.Data = np.zeros((Data.shape[0],2))
                if temp_obj.Type == 'Crossline':
                    temp_obj.Data[:,0] = Data[:,0]
                    temp_obj.Data[:,1] = Data[:,3]
                elif temp_obj.Type == 'Inline':
                    temp_obj.Data[:, 0] = Data[:, 1]
                    temp_obj.Data[:, 1] = Data[:, 3]
                elif temp_obj.Type == 'PDD':
                    temp_obj.Data[:, 0] = Data[:, 2]
                    temp_obj.Data[:, 1] = Data[:, 3]
                else:
                    temp_obj.Data[:, 0] = np.sqrt(Data[:, 0]**2 +Data[:, 1]**2)

                    boundary = np.argmin( temp_obj.Data[:, 0] - 0 )
                    temp_obj.Data[0:boundary, 0] = -temp_obj.Data[0:boundary, 0]

                    temp_obj.Data[:, 1] = Data[:, 3]


                temp_obj.SSD = float(strfindlr(cure_infor, 'Source to Surface Distance (cm)	', '\n')[0])*10
                temp_obj.Date = strfindlr(cure_infor, 'Scan Date	', '\t\n')

                temp_obj.Depth = 0 if temp_obj.Type=='PDD' else np.round(Data[0,2]*10,1)

                all_cure_list.append(temp_obj)

    return all_cure_list

# def selectCure(all_cures,conditions):
#
#     flag = []
#     for key in conditions:
#         flag.append( [eval(f'cure.{key}') in conditions[key] for cure in all_cures])
#
#     if len(flag)==0:
#         temp = np.array(flag[0])
#     else:
#         for id in range(len(flag)-1):
#             if id == 0:
#                 temp = np.array(flag[0]) & np.array(flag[1])
#             else:
#                 temp = temp & np.array(flag[id+1])
#
#
#     return np.array(all_cures)[temp]


def renormalized(select_cures):
    a = 95.6  # 6MV FFF
    b = 0.6595
    c = 0.1255
    d = -0.0099
    e = 0.0013

    # # 10MV  FFF
    # a = 89.08
    # b = 2.4826
    # c = 0.1152
    # d = -0.0078
    # e = 0.0011

    for id in range(len(select_cures)):
        if select_cures[id].Type == 'PDD':
            continue
        else:
            if select_cures[id].Type == 'Inline' :
                FS = select_cures[id].FS[1]/10  # unit cm
                Depth = select_cures[id].Depth
                center = np.interp(0, select_cures[id].Data[:,0],select_cures[id].Data[:,1] )
                normal_value =( a + b * FS  + c * Depth) / (1 + d * FS + e * Depth )/100
                select_cures[id].Data[:,1] = select_cures[id].Data[:,1]  /center  * normal_value

            elif select_cures[id].Type == 'Crossline' :
                FS = select_cures[id].FS[0]/10  # unit cm
                Depth = select_cures[id].Depth
                center = np.interp(0, select_cures[id].Data[:,0],select_cures[id].Data[:,1] )
                normal_value =( a + b * FS  + c * Depth) / (1 + d * FS + e * Depth )/100
                select_cures[id].Data[:,1] = select_cures[id].Data[:,1]  /center  * normal_value

    return select_cures








if __name__ =='__main__':

    path_1 = r"E:\1工作列表\2022年\20220607-测试中寻找电子线各个影响因素B3\20220603-Bunker3\12MeV-APP20-JAW25-PirmaryFoil12MeV-SecondaryFoil12-15MeV.snctxt"
    path_2 = r"E:\1工作列表\2022年\20220607-测试中寻找电子线各个影响因素B3\20220603-Bunker3\12MeV-APP20-JAW28-PirmaryFoil12MeV-SecondaryFoil12-15MeV.snctxt"
    # path_3 = r"E:\1工作列表\2022年\20220607-测试中寻找电子线各个影响因素B3\20220603-Bunker3\12MeV-APP20-JAW32-PirmaryFoil12MeV-SecondaryFoil12-15MeV.snctxt"

    all_cures_1 = readSNCfile(path_1)
    all_cures_2 = readSNCfile(path_2)
    # all_cures_3 = readSNCfile(path_3)


    conditions = {
                  'FS': [[250, 250]],
                'Type': 'Crossline',
                'Depth':[21.2]
                  }
    #  'Type': 'Diagonal(++/--)',
    method = 'varian'    #  inflection  or  varian
    beam_mode = 'FFF'
    select_cures_1 = selectCure(all_cures_1,conditions)
    select_cures_2 = selectCure(all_cures_2,conditions)
    # select_cures_3 = selectCure(all_cures_3,conditions)
    # if method == 'varian':
    #     select_cures =  renormalized(select_cures)
    # / cure1.Data[int(len(cure1.Data[:, 1]) / 2)
    for cure1,cure2 in zip(select_cures_1,select_cures_2):
        plt.plot(cure1.Data[:, 0], cure1.Data[:, 1]/cure1.Data[int(len(cure1.Data[:,1])/2) , 1],label=f'{os.path.split(path_1)[-1]}')
        plt.plot(cure2.Data[:, 0], cure2.Data[:, 1]/cure2.Data[int(len(cure2.Data[:,1])/2) , 1],label=f'{os.path.split(path_2)[-1]}')
        # plt.plot(cure3.Data[:, 0], cure3.Data[:, 1]/cure3.Data[int(len(cure3.Data[:,1])/2) , 1],label=f'{os.path.split(path_3)[-1]}')

    plt.legend(loc=0,prop={'size':10})
    plt.xlabel('off-axis (cm)')
    plt.ylabel('Relative dose')
    plt.title(f'12MeV束流，不同初级散射箔，Crossline 影响（Depth={cure1.Depth}mm） ')
    # plt.grid()
    plt.show()

    # for cure1 in select_cures_2:
    #     plt.plot(cure1.Data[:, 0], cure1.Data[:, 1] ,label=f'{os.path.split(path_1)[-1]}')
    #
    # plt.show()