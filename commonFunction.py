
import os
import re
import  numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt
from astropy.modeling import models,fitting
plt.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
plt.rcParams['axes.unicode_minus'] = False  # 坐标轴负号显示正常


def findEveryFolder(path,kew_word):
    result_list = []
    for root,dirs,files in os.walk(path):
        for file in files:
            if  kew_word in file:
                result_list.append( os.path.join(root,file) )

    return result_list


def findCurrentFolder(path,kew_word):
    result_list = []
    for file  in os.scandir( path):
        if file.is_file()  and (kew_word in file.name):
            result_list.append(file.path )

    return result_list


def strfindlr(str,left_str,right_str):
    all_parts = []


    left_str = left_str.replace('(','\(')
    left_str = left_str.replace(')','\)')
    #
    # right_str_re = right_str.replace('(','\(')
    # right_str_re = right_str_re.replace(')','\)')
    #
    # start_positon_list =   [i.start() for i in re.finditer(left_str_re,str)]
    # end_positon_list   =   [i.start() for i in re.finditer(right_str_re,str)]
    #
    # if len(start_positon_list) == len(end_positon_list):
    #     for start_index,end_index in zip(start_positon_list,end_positon_list):
    #         all_parts.append(str[start_index+len(left_str):end_index] )
    #
    # elif  len(start_positon_list) > len(end_positon_list):
    #     # start str > end str
    #     new_start_positon_list = []
    #     for ide,end_index in enumerate(end_positon_list)  :
    #         for ids,start_index in enumerate(start_positon_list):
    #             if int(start_index) > int(end_index):
    #                 new_start_positon_list.append( start_positon_list[ids-1] )
    #                 break
    #
    #     for start_index,end_index in zip(new_start_positon_list,end_positon_list):
    #         all_parts.append(str[start_index+len(left_str):end_index] )
    #
    # else:
    #     # start str < end str
    #     new_end_positon_list = []
    #     for ids, start_index in enumerate(start_positon_list):
    #         for ide, end_index in enumerate(end_positon_list):
    #             if int(end_index) > int(start_index):
    #                 new_end_positon_list.append(end_index)
    #                 break
    #     for start_index, end_index in zip(start_positon_list, new_end_positon_list):
    #         all_parts.append(str[start_index + len(left_str):end_index])

    contents = re.finditer(r"(?<=" + left_str + ").*?" + "(?=" + right_str + ")", str, flags=re.DOTALL)
    all_parts = [i.group() for i in contents]


    return all_parts



def calFieldSize(beamMode,x,y, method='iba',field_size_percent=0.5):
    # from scipy.interpolate import make_interp_spline
    # cubic = si.interp1d(x,y)
    # xx = np.arange( np.min(x), np.max(x)+0.1, 0.1)
    # yy = cubic(xx)
    # xx = np.arange(np.min(x), np.max(x) + 0.1, 0.1)
    # yy = make_interp_spline(x, y)(xx)
    # #
    #
    yy = y
    xx = x

    center = yy[int(len(yy)/2)]
    halfLen = int(len(yy)/2)

    rg = 5


    if beamMode=='FFF'  and method == 'inflection':

        err = yy[1:] - yy[:-1]




        leftBundary_temp = np.argmax(err)
        rightBundary_temp = np.argmin(err)
        range =5

        left_reflection_point = getPeakGaussPoint(xx[leftBundary_temp-range:leftBundary_temp+range+2],abs(err[leftBundary_temp-range:leftBundary_temp+range+2]))
        right_reflection_point = getPeakGaussPoint(xx[rightBundary_temp-range:rightBundary_temp+range+2],abs(err[rightBundary_temp-range:rightBundary_temp+range+2]))

        # left_reflection_point = -54.9
        # right_reflection_point = 55.1

        # left_reflection_point = (xx[np.argmax(err)] + xx[np.argmax(err)+1])/2
        # right_reflection_point = (xx[np.argmin(err)] + xx[np.argmin(err)+1])/2
        #
        # left_reflection_point = xx[np.argmax(err)]
        # right_reflection_point = xx[np.argmin(err)]

        left_reflection_vlue =  np.interp( left_reflection_point, xx[leftBundary_temp - rg:leftBundary_temp + rg],yy[leftBundary_temp - rg:leftBundary_temp + rg])
        right_reflection_vlue = np.interp(right_reflection_point, xx[rightBundary_temp - rg:rightBundary_temp + rg],
                                    yy[rightBundary_temp - rg:rightBundary_temp + rg])



        left_renormalized_value  = yy / left_reflection_vlue * 0.5 # 按照50% 点进行归一
        right_renormalized_value = yy / right_reflection_vlue * 0.5  # 按照50% 点进行归一

        leftBundary_temp = np.argmin((abs(left_renormalized_value[0:halfLen + 1] - field_size_percent)))  # 得到20% 附近的点
        leftBundary = np.interp(field_size_percent, left_renormalized_value[leftBundary_temp - rg:leftBundary_temp + rg],
                                xx[leftBundary_temp - rg:leftBundary_temp + rg])


        rightBundary_temp = np.argmin((abs(right_renormalized_value[halfLen + 1:] - field_size_percent))) + halfLen + 1
        rightBundary = np.interp(field_size_percent, right_renormalized_value[rightBundary_temp - rg:rightBundary_temp + rg],
                                 xx[rightBundary_temp - rg:rightBundary_temp + rg])


    elif  beamMode=='FFF' and  method == 'varian':

        leftBundary_temp = np.argmin((abs(yy[0:halfLen + 1] - field_size_percent)))
        rightBundary_temp = np.argmin((abs(yy[halfLen + 1:] - field_size_percent))) + halfLen + 1

        leftBundary = np.interp(field_size_percent, yy[leftBundary_temp - rg:leftBundary_temp + rg],xx[leftBundary_temp - rg:leftBundary_temp + rg])  # 将拐点 归一化到50%点

        rightBundary = np.interp(field_size_percent, np.flip(yy[rightBundary_temp - rg:rightBundary_temp + rg]),np.flip(xx[rightBundary_temp - rg:rightBundary_temp + rg]))



    if  beamMode=='FF' :

        data = yy / center

        leftBundary_temp = np.argmin((abs(data[0:halfLen + 1] - field_size_percent)))
        rightBundary_temp = np.argmin((abs(data[halfLen + 1:] - field_size_percent))) + halfLen + 1

        leftBundary = np.interp(field_size_percent,
                                data[leftBundary_temp - rg:leftBundary_temp + rg], x[leftBundary_temp - rg:leftBundary_temp + rg])  # 将拐点 归一化到50%点

        rightBundary = np.interp(field_size_percent,
                                 np.flip(data[rightBundary_temp - rg:rightBundary_temp + rg]),np.flip(x[rightBundary_temp - rg:rightBundary_temp + rg]))




    return round( np.abs(rightBundary - leftBundary),3),  leftBundary,  rightBundary


def selectCure(all_cures,conditions):

    flag = []
    for key in conditions:
        if key == 'Depth':
            conditions[key] = [conditions[key]] if not type(conditions[key])==list else  conditions[key]
            expand = [i-0.1 for i  in conditions[key]] + [i+0.1 for i  in conditions[key]] + [i for i  in conditions[key]]
            flag.append([eval(f'cure.{key}') in  expand   for cure in all_cures])
        else:
            flag.append([eval(f'cure.{key}') in conditions[key] for cure in all_cures])

    if len(flag)==1:
        temp = np.array(flag[0])
    elif len(flag)>1:
        for id in range(len(flag)-1):
            if id == 0:
                temp = np.array(flag[0]) & np.array(flag[1])
            else:
                temp = temp & np.array(flag[id+1])


    return np.array(all_cures)[temp]


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


            if select_cures[id].Data[0,0] > select_cures[id].Data[-1,0] :
                select_cures[id].Data = np.flip( select_cures[id].Data,0)


            center = np.interp(0, select_cures[id].Data[:,0],select_cures[id].Data[:,1] )
            normal_value =( a + b * FS  + c * Depth) / (1 + d * FS + e * Depth )/100
            select_cures[id].Data[:,1] = select_cures[id].Data[:,1]  /center  #* normal_value   FFF

    return select_cures





import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib.pyplot as plt



def getPeakGaussPoint(x,y):
    plt.figure()
    plt.plot(x, y, 'ro')

    def func(params):
        alpha, mu, sigma = tuple(params)
        return alpha * 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-(x - mu) ** 2 / sigma ** 2)

    def func_residual(params):
        return func(params) - y

    def final(x, params):
        alpha, mu, sigma = tuple(params)

        return alpha * 1.0 / np.sqrt(2 * np.pi) / sigma * np.exp(-(x - mu) ** 2 / sigma ** 2)

    fit_result = least_squares(
        func_residual,
        [10.0,x[int(len(x)/2)], 10],
        jac='3-point',
        bounds=[[0.001, x[0], 0.1], [10000, x[-1], 200]],
        loss='linear',
        tr_solver='exact',
        verbose=1
    )

    print(fit_result)
    xfine = np.linspace(x.min(), x.max(), 1001)
    plt.plot(xfine, final(xfine, fit_result.x), 'b-')
    return fit_result.x[1]
