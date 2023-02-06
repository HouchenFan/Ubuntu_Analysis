import os
import numpy as  np
import matplotlib.pyplot as plt

def read_egsphantom(egs_phantom_path):

    with open(egs_phantom_path, 'r') as f:
        material_num = int(f.readline())
        material_list = []
        for i in range(material_num):
            material_list.append(f.readline().split()[0])

        f.readline()

        vox_num = [int(i) for i in f.readline().split()]
        x_boundary = [float(i) for i in f.readline().split()]
        y_boundary = [float(i) for i in f.readline().split()]
        z_boundary = [float(i) for i in f.readline().split()]

        vox_material_id = np.zeros((vox_num[2], vox_num[1], vox_num[0]))

        for id_z in range(vox_num[2]):
            for id_y in range(vox_num[1]):
                xlist = f.readline().split()[0]
                for id_x, x in enumerate(xlist):
                    vox_material_id[id_z, id_y, id_x] = int(x)

            f.readline()  # 每一层有个空行

        vox_density = np.zeros((vox_num[2], vox_num[1], vox_num[0]))

        for id_z in range(vox_num[2]):
            for id_y in range(vox_num[1]):
                xlist = f.readline().split()
                for id_x, x in enumerate(xlist):
                    vox_density[id_z, id_y, id_x] = float(x)
            f.readline()  # 每一层有个空行

        fig,axes = plt.subplots(1,3,figsize=(18,6))
        (dim_z,dim_x,dim_y) = vox_density.shape
        i_z = int(dim_z/2)
        i_x = int(dim_x/2)
        i_y = int(dim_y/2)

        # x = np.linspace(-, 5, 500)
        # y = np.linspace(-5, 5, 500)
        X, Y = np.meshgrid(x_boundary[0:-1], y_boundary[0:-1])
        axes[0].contourf(X, Y, vox_density[70, :, :])
        axes[0].imshow(vox_density[70,:,:])
        axes[0].set_title('XYplan')
        axes[1].imshow(vox_density[:, i_x, :])
        axes[1].set_title('YZplan')
        axes[2].imshow(vox_density[:, :, i_y])
        axes[2].set_title('XZplan')

        plt.show()


if __name__ == '__main__':
    egs_phantom_path = r"/home/uih/Head/Head_default.egsphant"
    read_egsphantom(egs_phantom_path)



