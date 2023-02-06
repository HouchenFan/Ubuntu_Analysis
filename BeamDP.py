import subprocess
import os, re


def Replace(line, idx, new_info):
    assert isinstance(line, str)
    new_line = line.split(',')
    assert idx < len(new_line)
    new_line[idx] = new_info
    return ','.join(new_line)


class Info(object):
    def __init__(self, line, position, info):
        self.line = line
        self.position = position
        self.info = info


def ModifyFile(path, new_path, Info_list):
    with open(path, 'r') as fid:
        txt = fid.readlines()
    # Modify
    for i in Info_list:
        assert isinstance(i, Info)
        txt[i.line] = Replace(txt[i.line], i.position, i.info)
        if txt[i.line][-1] != '\n':
            txt[i.line] += '\n'
    with open(new_path, 'w') as fid:
        fid.writelines(txt)


def Inputs(input_info):
    if isinstance(input_info, list):
        txt = input_info
    elif os.path.isfile(input_info):
        with open(input_info, 'r') as fid:
            txt=fid.readlines()
    else:
        raise ValueError("Unvalid input.")
    inputs = ''.join(txt)
    inputs = bytes(inputs, 'utf-8')
    return inputs


def Thickness(filename):
    return float(re.search("\d+\.\d+cm", filename)[0][:-2])


def BeamDP(input_info):
    print("BeamDP")
    p = subprocess.Popen(['/home/uih/EGSnrc/HEN_HOUSE/bin/linux64/beamdp'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    inputs = Inputs(input_info)
    output, err = p.communicate(inputs)
    print('Output: ', output.decode('utf-8'))
    print('Error: ', err.decode('utf-8'))
    print('Exit code:', p.returncode)


def BeamDPSet(phsp_dir, input_path):
    phsp_files = list(os.listdir(phsp_dir))
    phsp_files = list(filter(lambda x: os.path.splitext(x)[1][:-1] == '.egsphsp', phsp_files))
    # phsp_files = sorted(phsp_files, key=Thickness)
    # phsp_files = phsp_files[50:]
    print(phsp_files)
    for idx,pf in enumerate(phsp_files):
        full_path = os.path.join(phsp_dir, pf)
        i_phsp = Info(3, 0, full_path)
        i_agr = Info(6, 0, os.path.splitext(full_path)[0] + '_Spe.agr')
        i_radius = Info(2, 2, str(0.5))
        i_energy = Info(3, 2, phsp_dir.split(os.sep)[-2].split('MV')[0])
        script = os.path.join(phsp_dir, 'beamdp.script')
        ModifyFile(input_path, script, [i_phsp, i_agr, i_radius, i_energy])
        BeamDP(script)


if __name__ == "__main__":
    # phsp_dir = [r"/home/uih/Data/15MV/W_new", r"/home/uih/Data/20MV/W_new", r"/home/uih/Data/30MV/W_new", r"/home/uih/Data/50MV/W_new", r"/home/uih/Data/100MV/W_new"]
    phsp_dir = [r"/home/uih/Data/150MV/W_new"]
    inputs_path = r"/home/uih/Data/15MV/W_new/beamdpguitemp.script"
    for ph in phsp_dir:
        BeamDPSet(ph, inputs_path)
    print("Finished.")
