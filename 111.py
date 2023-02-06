import numpy as np
import matplotlib.pyplot as plt

x= np.linspace(-4,9,14)
print(x)
print( np.where(x<0 ,x**2,x**3 )  )
print([x<3,x<5,x<9] )


# fig,ax = plt.subplots()
# ax.plot(x,y,c='r',lw=1,alpha=0.8)
# ax.set_xlabel('off-axis (cm)',fontsize=9,color='g')
# ax.set_ylabel('value')
# plt.ioff()
# plt.show()