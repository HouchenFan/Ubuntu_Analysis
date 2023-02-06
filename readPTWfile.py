
# read ptw file(*.mcc)
# input  dir or full file name
# output: list  , each item  is obj.  the obj contain cures informationin .like ,depth,ssd,type,ect..

import os
import numpy as np
import matplotlib.pyplot as plt
from commonFunction import *

class PTW_Cure(object):  # 32102.703
    pass

def getCureType(cure_infor):

    if 'PDD' in strfindlr(cure_infor,'SCAN_CURVETYPE=','\n') :
        return 'PDD'

    elif 'CROSSPLANE_PROFILE' in strfindlr(cure_infor, 'SCAN_CURVETYPE=', '\n') and 'NOT_DIAGONAL' in strfindlr(cure_infor,'SCAN_DIAGONAL=','\n'):
        return 'Crossline'


    elif 'INPLANE_PROFILE' in strfindlr(cure_infor,'SCAN_CURVETYPE=','\n')  and 'NOT_DIAGONAL' in strfindlr(cure_infor,'SCAN_DIAGONAL=','\n' ):
        return 'Inline'

    elif 'INPLANE_PROFILE' in strfindlr(cure_infor, 'SCAN_CURVETYPE=', '\n') and 'FIRST_DIAGONAL' in strfindlr(cure_infor,
                                                                                                             'SCAN_DIAGONAL=',
                                                                                                             '\n'):
        return 'Diagonal(++/--)'

    elif 'CROSSPLANE_PROFILE' in strfindlr(cure_infor, 'SCAN_CURVETYPE=', '\n') and 'FIRST_DIAGONAL' in strfindlr(cure_infor,
                                                                                                             'SCAN_DIAGONAL=',
                                                                                                             '\n'):
        return 'Diagonal(-+/+-)'


def readPTWfile(dir_or_file):
    ptw_file_list = []
    if os.path.isdir(dir_or_file):
        ptw_file_list = findEveryFolder(dir_or_file,'.mcc')

    else:
        ptw_file_list.append(dir_or_file)

    all_cure_list  = []

    for ptw_file in ptw_file_list:

        with open(ptw_file,'r') as f:
            all_str = f.read()
            all_cure_infor_list = strfindlr(all_str,'BEGIN_SCAN ','END_SCAN ')  # 得到每一条cure的 分段信息
            for cure_infor in all_cure_infor_list:
                tem_obj = PTW_Cure()

                tem_obj.Type =  getCureType(cure_infor)
                X_FS = int(float(strfindlr(cure_infor, '\tFIELD_CROSSPLANE=', '\n')[0]))
                Y_FS = int(float(strfindlr(cure_infor, '\tFIELD_INPLANE=', '\n')[0]))
                tem_obj.FS = [X_FS,Y_FS]
                if not tem_obj.Type == 'PDD':
                    tem_obj.Depth = float(strfindlr(cure_infor, '\tSCAN_DEPTH=', '\n')[0])
                else:
                    tem_obj.Depth = 0  # PDD 的深度字段置为0

                temp = strfindlr(cure_infor, 'BEGIN_DATA', 'END_DATA')
                temp = re.sub('#\d','',temp[0])
                temp = re.sub('\s+', ' ', temp)
                temp = temp.split(' ')
                temp = [float(i)  for i in temp if i]

                if np.mod(len(temp),2) == 0 :
                    tem_obj.Data  = np.array(temp).reshape(-1,2)
                elif np.mod(len(temp),3) == 0 :
                    tem_obj.Data  = np.array(temp).reshape(-1,3)

                tem_obj.SSD = float(strfindlr(cure_infor, 'SSD=', '\n')[0])
                tem_obj.Date = strfindlr(cure_infor, 'FILE_CREATION_DATE=', ' ')
                X_offset = float(strfindlr(cure_infor, 'COLL_OFFSET_CROSSPLANE=', '\n')[0])
                Y_offset = float(strfindlr(cure_infor, 'COLL_OFFSET_INPLANE=', '\n')[0])
                tem_obj.offset = [X_offset,Y_offset]

                all_cure_list.append(tem_obj)

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

            elif select_cures[id].Type == 'Crossline' :
                FS = select_cures[id].FS[0]/10  # unit cm
                Depth = select_cures[id].Depth


            center = np.interp(0, select_cures[id].Data[:,0],select_cures[id].Data[:,1] )
            normal_value =( a + b * FS  + c * Depth) / (1 + d * FS + e * Depth )/100
            select_cures[id].Data[:,1] = select_cures[id].Data[:,1]  /center  * normal_value

    return select_cures








if __name__ =='__main__':

    path = r"E:\1工作列表\2022年\20220527-半影评估\Unity\600026_STH_Final_20200529 2\600026_STH_Final_20200529\STH600026 pdd and profile\G0_AllScans.mcc"

    # path =r"E:\1工作列表\2022年\20220527-半影评估\Unity\ZS_unity\Raw Data\Day 1\G000_MicroDiamond_Processed\X07 OPEN 10X10 IN TBA 190611 07'47'02.mcc"

    # path = r"E:\1工作列表\2022年\20220527-半影评估\Unity\ZS_unity\Raw Data\Day 1\G000_MicroDiamond_Processed\X07 OPEN 10X10 IN TBA 190611 07'47'02.mcc"
    all_cures = readPTWfile(path)


    conditions = {'Type':['Crossline','Inline'],
                  'Depth':[100,200],
                  'FS': [[100, 100]],
                  }
    # # method = 'inflection'    #  inflection  or  varian
    # method = 'varian'    #  inflection  or  varian
    # beam_mode = 'FFF'

    select_cures = selectCure(all_cures,conditions)
    # if method == 'varian':
    #     select_cures =  renormalized(select_cures)

    for cure in select_cures:
        plt.plot(cure.Data[:, 0], cure.Data[:, 1])
        # lr_20 = calFieldSize(beam_mode, cure.Data[:, 0], cure.Data[:, 1], method=method, field_size_percent=0.2)
        # lr_80 = calFieldSize(beam_mode, cure.Data[:, 0], cure.Data[:, 1], method=method, field_size_percent=0.8)
        #
        # pen_left = np.round( lr_80[1] - lr_20[1] ,2)
        # pen_right = np.round( lr_20[2] - lr_80[2] ,2 )
        # print(f'Fieldsize {cure.FS} left penubra is {pen_left} mm ------ right penubra is {pen_right} mm  ----- average:{np.round((pen_right+pen_left)/2,2)}mm')

    plt.show()


