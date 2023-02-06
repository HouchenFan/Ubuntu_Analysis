from  readPTWfile import *
from  readIBAfile import *

unity_conditions = {'FS':[[570,220]],
                 'Type':  ['Crossline'],
                 'Depth': [100]
              }

versaHD_conditions = {'FS':[[400,400]],
                 'Type': ['Crossline'],
                'Depth': [100]
              }


unity_path = r"E:\1工作列表\2022年\20220527-半影评估\Unity\600026_STH_Final_20200529 2\600026_STH_Final_20200529\STH600026 pdd and profile\G0_AllScans.mcc"
all_unity_cures = readPTWfile(unity_path)
select_cures = selectCure(all_unity_cures, unity_conditions)
unity_fs10 = renormalized(select_cures)




versa_path = r"M:\10.8.61.47\versaHD-FFF-F10.asc"
all_versa_cures = readIBAfile(versa_path)
select_cures = selectCure(all_versa_cures, versaHD_conditions)
versa_fs10 = renormalized(select_cures)


for unity,versa in zip(unity_fs10,versa_fs10):
    plt.plot(unity.Data[:, 0], unity.Data[:, 1],label='Unity 10x10 ,SSD=133.5,Depth=10cm,[Diamond]')
    plt.plot(versa.Data[:, 0], versa.Data[:, 1],label='Versa 10x10 ,SSD=100,Depth=10cm [CC13]')

plt.legend(loc=0,prop={'size':8})
plt.xlabel('off-axis (cm)')
plt.ylabel('Relative dose')
plt.title('Unity vs VersaHD comparison ')
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
