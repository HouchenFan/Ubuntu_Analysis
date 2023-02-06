import os
import numpy as np


def WritePhantom(phantom_path, dim_x, dim_y, spacing_x, spacing_y, material_list, material_density_list, start_z, thickness_list, layer_index):
    half_x = dim_x / 2
    half_y = dim_y / 2
    border_x_list = [spacing_x * (i - half_x) for i in range(dim_x + 1)]
    border_y_list = [spacing_y * (i - half_y) for i in range(dim_y + 1)]
    border_z_list = [start_z]
    for idx, thickness in enumerate(thickness_list):
        new_border = start_z + np.sum(thickness_list[:idx + 1])
        border_z_list.append(new_border)
    with open(phantom_path, "wt") as fid:
        fid.write(str(len(material_list))+"\n")
        for key in material_list.keys():
            fid.write(key+"\n")
        for i in range(len(material_list)):
            fid.write('0\n')
        fid.write("%d, %d, %d\n" % (dim_x, dim_y, len(layer_index)))
        for x in border_x_list:
            fid.write("%.6f " % x)
        fid.write('\n')
        for y in border_y_list:
            fid.write("%.6f " % y)
        fid.write("\n")
        for border_z in border_z_list:
            fid.write("%.4f " % border_z)
        fid.write("\n")
        for layer in layer_index:
            for i in range(dim_x):
                for j in range(dim_y):
                    fid.write("%d" % layer)
                fid.write("\n")
            fid.write("\n")

        for layer in layer_index:
            for i in range(dim_x):
                for j in range(dim_y):
                    fid.write("%.4f " % material_density_list[layer-1])
                fid.write("\n")
            fid.write("\n")


def EGS_PH(phantom_dir):
    dim_x = 99
    dim_y = dim_x
    spacing_x = 0.04
    spacing_y = 0.04
    start_z = 144.017
    material_list = {"170C": 1.7000, "226C": 2.2600, "462GD2O2S": 4.6200, "AIR521ICRU": 0.0012, "AL": 2.6990, "CELLULOSE_ACETATE": 1.4200, "CU": 8.9600, "POLYESTER_FILM": 1.3596, "SI": 2.3300}
    material_density_list = [1.7000, 2.2600, 4.6200, 0.0012, 2.6990, 1.4200, 8.9600, 1.3596, 2.3300]
    layer_index = [5, 4, 2, 8, 3, 6, 9, 1, 5, 4, 7, 4, 5]
    thickness_list = [0.0760, 0.7530, 0.0520, 0.0180, 0.12, 0.0070, 0.1100, 0.3000, 0.2000, 1.5210, 0.0100, 1.0290, 0.3000]
    save_path = os.path.join(phantom_dir, "EPID_Phantom_DimX%d_DimY%d_T%.2f.egsphant" % (dim_x, dim_y, thickness_list[4]))
    WritePhantom(save_path, dim_x, dim_y, spacing_x, spacing_y, material_list, material_density_list, start_z, thickness_list, layer_index)


def EGS_IMG(phantom_dir):
    dim_x = 128
    dim_y = dim_x
    # spacing_x = 0.04
    # spacing_y = 0.04
    spacing_x = 40.96 / dim_x
    spacing_y = 40.96 / dim_y
    start_z = 144.0
    material_list = {"170C": 1.7000, "226C": 2.2600, "462GD2O2S": 4.6200, "AIR521ICRU": 0.0012, "AL": 2.6990, "CELLULOSE_ACETATE": 1.4200, "ABS": 1.05, "POLYESTER_FILM": 1.3596, "SI": 2.3300}
    material_density_list = [1.7000, 2.2650, 4.6200, 0.0012, 2.6990, 1.4200, 1.05, 1.3596, 2.3300]
    layer_index = [7, 4, 5, 4, 2, 8, 3, 6, 9, 1, 5, 4]
    # thickness_list = [0.5, 1.017, 0.076, 0.753, 0.052, 0.018, 0.040, 0.007, 0.110, 0.300]
    thickness_list = [0.3, 2.45, 0.075, 0.1, 0.052, 0.018, 0.0423, 0.007, 0.110, 0.300, 0.2000, 1.5210]
    save_path = os.path.join(phantom_dir, "EPID_Phantom_IMG_DimX%d_DimY%d_T%.2f.egsphant" % (dim_x, dim_y, thickness_list[4]))
    WritePhantom(save_path, dim_x, dim_y, spacing_x, spacing_y, material_list, material_density_list, start_z, thickness_list, layer_index)

def EGS_SW(phantom_dir):
    dim_x = 80
    dim_y = dim_x
    # spacing_x = 0.04
    # spacing_y = 0.04
    spacing_x = 40 / dim_x
    spacing_y = 40 / dim_y
    start_z = 0
    material_list = {"SolidWater": 1.0450, "AIR521ICRU": 0.0012}
    material_density_list = [1.0450, 0.0012]
    for i in [0, 1, 5, 10, 15, 20, 25, 30]:
        layer_index = list(np.ones(i*2, dtype=int))
        # thickness_list = [0.5, 1.017, 0.076, 0.753, 0.052, 0.018, 0.040, 0.007, 0.110, 0.300]
        thickness_list = list(np.ones(i*2, dtype=float) * 0.5) #mm
        save_path = os.path.join(phantom_dir, "SW_Phantom_SW_DimX%d_DimY%d_T%.2f.egsphant" % (dim_x, dim_y, sum(thickness_list)))
        WritePhantom(save_path, dim_x, dim_y, spacing_x, spacing_y, material_list, material_density_list, start_z, thickness_list, layer_index)


def EGS_AP16(phantom_dir):
    # dim_x = 128
    dim_x = 128
    dim_y = dim_x
    #spacing_x = 0.01
    #spacing_y = 0.01
    spacing_x = 40.96 / dim_x
    spacing_y = 40.96 / dim_y
    # start_z = 144.0
    start_z = 0
    # material_list = {"AlM": 1.7000, "226C": 2.2600, "462GD2O2S": 4.6200, "AIR521ICRU": 0.0012, "AL": 2.6990,
    #                  "CELLULOSE_ACETATE": 1.4200, "ABS": 1.05, "POLYESTER_FILM": 1.3596, "SI": 2.3300}
    material_list = {"ABS": 1.05, "AL521ICRU": 2.702, "AIR521ICRU": 0.0012, "CU521ICRU": 8.933, "LD45": 0.045, "226C": 2.265, "PET": 1.38, "462GD2O2S": 4.62, "SI521ICRU": 2.33}
    material_density_list = [x for x in material_list.values()]
    material_name_list = [x for x in material_list.keys()]
    layer_index = [1, 3, 2, 3, 5, 2, 6, 2, 7, 8, 7, 9, 5, 2, 3, 4, 3, 2]
    # layer_index = [1, 3, 2, 3, 5, 2, 6, 2, 8, 7, 9, 5, 2]
    thickness_list = [3.0, 24.5, 0.75, 1.0, 2.5, 0.025, 0.9, 0.025, 0.009, 0.436, 0.188, 1.1, 1.5, 1.5, 4, 2.5, 3, 3]  # mm
    # thickness_list = [3.0, 24.5, 0.75, 1.0, 2.5, 0.025, 0.9, 0.025, 0.436, 0.188, 1.1, 1.5, 1.5]  # mm
    thickness_list = [x / 10.0 for x in thickness_list]  # mm to cm
    assert len(layer_index) == len(thickness_list)
    for i in range(len(layer_index)):
        print("Layer %2d: Material: %12s, Density: %.3f g/cm^3, Thickness: %.3f cm" % (i+1, material_name_list[layer_index[i] - 1], material_density_list[layer_index[i] - 1], thickness_list[i]))
    save_path = os.path.join(phantom_dir, "EPID_Phantom_AP16_DimX%d_DimY%d_T%.2f.egsphant" % (dim_x, dim_y, thickness_list[9]))
    WritePhantom(save_path, dim_x, dim_y, spacing_x, spacing_y, material_list, material_density_list, start_z, thickness_list, layer_index)


def Main(phantom_dir):
    material_list = {"170C": 1.7000, "226C": 2.2600, "462GD2O2S": 4.6200, "AIR": 0.0012, "AL": 2.6990, "CELLULOSE_ACETATE": 1.4200, "CU": 8.9600, "POLYESTER_FILM": 1.3596, "SI": 2.3300}
    material_density_list = [1.7000, 2.2600, 4.6200, 0.0012, 2.6990, 1.4200, 8.9600, 1.3596, 2.3300]
    dim_x = 99
    dim_y = dim_x
    half_x = dim_x/2
    half_y = dim_y/2
    # spacing_x = 40.96 / dim_x
    # spacing_y = 40.96 / dim_y
    spacing_x = 0.04
    spacing_y = 0.04
    layer_index = [5, 4, 2, 8, 3, 6, 9, 1, 5, 4, 7, 4, 5]
    start_z = 144.017
    # density_list = [2.6990, 0.0012, 0.0012, 0.0012, 2.2600, 1.3596, 4.6200, 1.4200, 2.3300, 1.7000, 2.6990, 0.0012, 8.9600, 0.0012, 2.6990]
    # thickness_list = [144.017000, 144.093000, 144.662000, 144.762000, 144.846000, 144.898000, 144.916000, 144.945000, 144.952000, 145.062000, 145.362000, 145.562000, 147.083000, 147.093000, 148.122000, 148.422000]
    # thickness_list = [0.0759999999999934, 0.5690000000000168, 0.09999999999999432, 0.08400000000000318, 0.0519999999999925, 0.018000000000000682, 0.028999999999996362, 0.007000000000005002, 0.11000000000001364, 0.29999999999998295, 0.20000000000001705, 1.5209999999999866, 0.009999999999990905, 1.0290000000000248, 0.29999999999998295]
    thickness_list = [0.0760, 0.7530, 0.0520, 0.0180, 0.12, 0.0070, 0.1100, 0.3000, 0.2000, 1.5210, 0.0100, 1.0290, 0.3000]
    border_z_list = [start_z]
    for idx, thickness in enumerate(thickness_list):
        new_border = start_z + np.sum(thickness_list[:idx+1])
        border_z_list.append(new_border)
    border_x_list = [spacing_x * (i - half_x) for i in range(dim_x + 1)]
    border_y_list = [spacing_y * (i - half_y) for i in range(dim_y + 1)]
    with open(os.path.join(phantom_dir, "EPID_Phantom_DimX%d_DimY%d_T%.2f.egsphant"%(dim_x, dim_y, thickness_list[4])), "wt") as fid:
        fid.write(str(len(material_list))+"\n")
        for key in material_list.keys():
            fid.write(key+"\n")
        for i in range(len(material_list)):
            fid.write('0\n')
        fid.write("%d, %d, %d\n" % (dim_x, dim_y, len(layer_index)))
        for x in border_x_list:
            fid.write("%.6f " % x)
        fid.write('\n')
        for y in border_y_list:
            fid.write("%.6f " % y)
        fid.write("\n")
        for border_z in border_z_list:
            fid.write("%.4f " % border_z)
        fid.write("\n")
        for layer in layer_index:
            for i in range(dim_x):
                for j in range(dim_y):
                    fid.write("%d" % layer)
                fid.write("\n")
            fid.write("\n")

        for layer in layer_index:
            for i in range(dim_x):
                for j in range(dim_y):
                    fid.write("%.4f " % material_density_list[layer-1])
                fid.write("\n")
            fid.write("\n")


if __name__ == "__main__":
    # Main(r"F:\Data\EGSnrc\EPID_Phantom")
    # EGS_IMG(r"F:\Data\EGSnrc\EPID_Phantom")
    #EGS_AP16(R"D:\Data_fhc\EGSnrc\EPID_Phantom")
    EGS_SW(R"D:\Data_fhc\EGSnrc\EPID_Phantom")