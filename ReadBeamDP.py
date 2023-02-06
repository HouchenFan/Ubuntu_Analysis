import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
plt.rcParams['axes.unicode_minus'] = False  # 坐标轴负号显示正常



class BeamDp():
    def __init__(self,file_path):
        self.file_path = file_path
        self.readbeamdp()



    def readbeamdp(self):
        with open(self.file_path) as f:
            all_lines = f.readlines()
            data = []

            for line in all_lines:
                if not(line.startswith('@') or line.startswith('&')):
                    for item in line.split('     '):
                        if item:
                            data.append(item)

        data = np.array(data, dtype='float64').reshape(-1, 3)

        self.file_name = os.path.split(self.file_path)[1]
        self.x = data[:, 0]
        self.y = data[:, 1]
        self.y_norm_by_area = self.y/np.sum(self.y)
        self.y_norm_by_max = self.y/np.max(self.y)
        # self.y_norm_by_first_point = self.y/self.y[0]
        self.y_norm_by_center_point = self.y / self.y[int(len(self.y)/2)]

        self.center = self.y[0]
        self.average = np.sum(self.x*self.y)/np.sum(self.y)
        self.total = np.sum(self.y)



def get_beamdp_list(dir_or_filelist,plot_type):
    beamdp_obj_list = []

    if isinstance(dir_or_filelist, list):
        #dir_or_filelist 是一个文件路径 列表
        for file in dir_or_filelist:
            beamdp_obj_list.append(BeamDp(file))

        plt.figure()
        for beamdp_obj in beamdp_obj_list:
            if not plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y,label=beamdp_obj.file_name)
            elif 'firt' in plot_type:
                plt.plot(beamdp_obj.x,beamdp_obj.y_norm_by_first_point)
            elif 'center' in plot_type:
                plt.plot(beamdp_obj.x,beamdp_obj.y_norm_by_center_point)
            elif 'max' in plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_max)
            elif 'area' in plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_area)
        plt.legend(loc=0)

    else:
        if os.path.isfile(dir_or_filelist):
            #dir_or_filelist 是一个文件
            beamdp_obj_list.append(BeamDp(dir_or_filelist))
        else:
            # dir_or_filelist 是一个文件夹
            for file in os.scandir(dir_or_filelist):
                if file.is_file and not file.name.endswith('jpg'):
                    beamdp_obj_list.append(BeamDp(file.path))

        plt.figure()
        for beamdp_obj in beamdp_obj_list:

            if not plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y,label=beamdp_obj.file_name)

                plt.xlabel('off-axis （cm） ')
                plt.ylabel('能注量')
                plt.title('不同屏蔽结构能注量分布')

            elif 'first' in plot_type:
                # plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_first_point,label=beamdp_obj.file_name)
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_area,label=beamdp_obj.file_name)
                plt.xlabel('Energy(MeV) ')
                plt.ylabel('Weight')
                plt.title('不同结构能能谱比较（SSD=10cm平面）')
                plt.legend()
            elif 'center' in plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_center_point)
            elif 'max' in plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_max)
            elif 'area' in plot_type:
                plt.plot(beamdp_obj.x, beamdp_obj.y_norm_by_area)

        plt.legend(loc=0)
        #
        # plt.plot(x, cal[plot_index[0], :], label='cal')
        # ax[0].plot(x, mea[plot_index[0], :], label='mea')
        # ax[0].set_title('Crossline')
        # ax[0].legend(loc='best')
        # ax[0].text(x=0, y=cal[256, 256] / 2, ha='center', va='baseline',
        #            s=f'left_err:{b1_cal - b1_mea:.3f}\nright_err:{b2_cal - b2_mea:.3f}')
        # ax[0].set_xlim(b1_cal - 3, b2_cal + 3)
        # if os.path.isfile(dir_or_filelist):
        #     dir = os.path.split(dir_or_filelist)[0]
        #     plt.savefig(os.path.join(dir, 'result.jpg'), dpi=200)
        #
        # else:
        #     plt.savefig(os.path.join(dir_or_filelist, 'result.jpg'), dpi=200)
        plt.show()
    return beamdp_obj_list











if __name__ == '__main__':
    # # dir_or_filelist = r"D:\1"
    # dir_or_filelist = r"D:\EGSnrc\egsnrc\BEAM_shielding\E-pb"
    # dir_or_filelist = r"D:\EGSnrc\egsnrc\BEAM_shielding\energy"

    plot_type = 'center_piont'
    plot_type = ''

    dir_or_filelist = ['/home/uih/AirWater/Air_IBL/EF_dos_air.agr',
                       '/home/uih/AirWater/Water_IBL/EF_dos_water.agr',
                       '/home/uih/AirWater/FS27_Sigma0.07_AR40_D15/EF_beam_air.agr',
                       '/home/uih/AirWater/FS27_Sigma0.07_SW40_D15/EF_beam_water.agr']

    get_beamdp_list(dir_or_filelist, plot_type)
    print()

