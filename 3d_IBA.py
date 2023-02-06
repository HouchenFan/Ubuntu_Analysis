from  read3ddose import *
from  readIBAfile import *

# unity_conditions = {'FS':[[570,220]],
#                  'Type':  ['Crossline'],
#                  'Depth': [100]
#               }

versaHD_conditions = {'FS':[[100,100]],
                 'Type': ['Crossline'],
                'Depth': [50]
              }


unity_path = r"D:\Measure\08102022\D5_6MV_FF_FS10_SSD100.3ddose"
# all_unity_cures = read3ddose(unity_path)
# select_cures = selectCure(all_unity_cures, unity_conditions)
# unity_fs10 = renormalized(select_cures)
cure_list =  read3ddose(unity_path)




versa_path = r"\\dataserver03\rt\06_PH\107-装机医院\110033-北京协和医院\01-水箱数据\协和asc\FF\2-Profile&PDD测量-MLC&Jaw开野-CC13探头"
all_versa_cures = readIBAfile(versa_path)
select_cures = selectCure(all_versa_cures, versaHD_conditions)
versa_fs10 = renormalized(select_cures)


for data,versa in zip(cure_list,versa_fs10):

    for data in cure_list:
        if len(data.PDD) > 2:
            plt.plot(data.cz, data.PDD*100, label=data.file_name)

    for data in cure_list:
        if len(data.profileXnorm) > 2:
            linestyle = '-'
            if 'Double' in data.file_name:
                linestyle = '--'
                # data.cx = data.cx - 0.6
            plt.plot(data.cx, data.profileXnorm, linestyle=linestyle, label=data.file_name)

    for data in cure_list:
        if len(data.profileYnorm) > 2:
            plt.plot(data.cy, data.profileYnorm, label=data.file_name)

    plt.plot(versa.Data[:, 0], versa.Data[:, 1],label='Versa 10x10 ,SSD=100, Depth=5 [CC13]')

plt.legend(loc=0,prop={'size':8})
plt.xlabel('off-axis (cm)')
plt.ylabel('Relative dose')
plt.title('3ddose vs VersaHD comparison ')
plt.show()


#
# for cure in iba_cures:
#         plt.plot(cure.Data[:, 0], cure.Data[:, 1])
#         lr_20 = calFieldSize('FFF', cure.Data[:, 0], cure.Data[:, 1], method='ptw', field_size_percent=0.2)
#         lr_80 = calFieldSize('FFF', cure.Data[:, 0], cure.Data[:, 1], method='ptw', field_size_percent=0.8)
#
#         pen_left = np.round( lr_80[1] - lr_20[1] ,2)
#         pen_right = np.round( lr_20[2] - lr_80[2] ,2 )
#         print(f'Fieldsize {cure.FS} left penubra is {pen_left} mm ------ right penubra is {pen_right} mm  ----- average:{np.round((pen_right+pen_left)/2,2)}mm')
#
#     plt.show()
