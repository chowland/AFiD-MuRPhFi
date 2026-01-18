import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
#import imageio

# 设置绘图样式
plt.style.use("ggplot")
xy_trans = False  # 设置是否转置坐标系
folder = r"."
path = folder + "/outputdir/flowmov/png/"

# 检查输出文件夹是否存在
if not os.path.exists(path):
    os.makedirs(path)
else:
    print(f"The folder '{path}' already exists.")

# 读取网格数据
class Grid:
    def __init__(self, xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr):
        self.xm = xm
        self.xmr = xmr
        self.xc = xc
        self.xcr = xcr
        self.ym = ym
        self.ymr = ymr
        self.yc = yc
        self.ycr = ycr
        self.zm = zm
        self.zmr = zmr

def read_grid(folder):
    with h5py.File(folder+"/outputdir/cordin_info.h5","r") as f:
        xm = f["xm"][()]
        xmr = f["xmr"][()]
        xc = f["xc"][()]
        xcr = f["xcr"][()]
        ym = f["ym"][()]
        ymr = f["ymr"][()]
        zm = f["zm"][()]
        zmr = f["zmr"][()]
    dy, dyr = ym[1] - ym[0], ymr[1] - ymr[0]
    yc, ycr = np.arange(0, dy*(ym.size+1), dy), np.arange(0, dyr*(ymr.size+1), dyr)
    return Grid(xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr)

# 读取 z-cut 数据
def read_zcut(folder, var, idx):
    varname = var+"/"+"%05d" % idx
    with h5py.File(folder+"/outputdir/flowmov/movie_zcut.h5","r") as f:
        A = f[varname][()]
    return A

# 计算样本数量
def count_samples(folder):
    with h5py.File(folder+"/outputdir/flowmov/movie_zcut.h5","r") as f:
        Nsamp = len(f["temp"].keys())
    return Nsamp

# 读取输入参数
class InputParams:
    def __init__(self, folder):
        fname = folder+"/bou.in"
        with open(fname,"r") as f:
            bou = f.readlines()
        if (True):  # 当前格式
            self.nxm, self.nym, self.nzm, self.nsst = [int(n) for n in bou[2].split()]
            self.multires, self.nxmr, self.nymr, self.nzmr = [int(n) for n in bou[6].split()]
            self.flagsal, self.flagpf = [n!="0" for n in bou[10].split()]
            self.tout, self.tframe = [float(n) for n in bou[22].split()[:2]]
            self.save_3D = float(bou[22].split()[-1])
            self.alx3, self.ylen, self.zlen = [float(n) for n in bou[26].split()]
            self.istr3 = int(bou[30].split()[0])
            self.str3 = float(bou[30].split()[1])
            self.istr3r = int(bou[30].split()[2])
            self.RayT, self.PraT, self.RayS, self.PraS = [float(n) for n in bou[34].split()[:-1]]
            self.FFscaleS = bou[34].split()[-1]=="1"
            self.inslwS, self.inslwN, self.TfixS, self.TfixN, self.SfixS, self.SfixN = [n!="0" for n in bou[42].split()]
            self.active_T, self.active_S = [n!="0" for n in bou[46].split()[:-1]]
            self.gAxis = int(bou[46].split()[-1])
            self.xplusU, self.xminusU, self.dPdy, self.dPdz = [float(n) for n in bou[50].replace("d","e").split()[:-1]]
            self.pf_D, self.pf_A, self.pf_S, self.pf_Tm = [float(n) for n in bou[54].split()[:4]]
            self.IBM = bou[54].split()[4]!="0"
            self.pf_IC = int(bou[54].split()[-1])

grid = read_grid(folder)
inputs = InputParams(folder)

# 读取流场数据
num = 0


# 计算样本数量
Nsamp = count_samples(folder)

# 输出 `tecplot` 格式数据
numls = list(range(0, Nsamp, 10))
path = folder + "/outputdir/flowmov/tecplot/"
if not os.path.exists(path):
    os.makedirs(path)

nxc = len(grid.xm)
nyc = len(grid.xm)
xc = grid.xc
xm = grid.xm
yc = grid.yc
ym = grid.ym

for num in numls:
    time = num
    T = read_zcut(folder, "temp", num)
    vx = read_zcut(folder, "vx", num)
    vy = read_zcut(folder, "vy", num)
    phi = read_zcut(folder, "phi", num)

    my_list = vx.T
    result = [[a + b for a, b in zip(sublist1, sublist2)] for sublist1, sublist2 in zip(my_list[:-1], my_list[1:])]
    vx = np.array(result).T

    filename = path + f'course_grid_{(1000000 + time * inputs.tframe)}.plt'
    with open(filename, 'w') as file:
        file.write('TITLE = "Flow Field Data"\n')
        file.write('VARIABLES = "X", "Y", "U", "V", "temp", "phi"\n')
        file.write(f'ZONE T="Flow Field Data", I={nxc}, J={nyc}, F=POINT\n')

        for i in range(nxc):
            for j in range(nyc):
                file.write(f'{float(xm[i])}\t  {float(ym[j])} \t {float(vy[i][j])} \t' +
                           f'{float(vx[i][j])}\t {float(T[i][j])}\t {float(phi[i][j])}\n')

    print(f"Finished saving Tecplot file for time step {time}.")

print("Tecplot files for all time steps have been generated.")
