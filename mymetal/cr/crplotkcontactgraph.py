import re
import os 
import numpy as np
import matplotlib.pyplot as plt

# 创建工具类
class Utils:
    
    """
    读取参数
    """
    def readParams(self, root):
        # 拼接参数文件路径
        path = root + "/properties.txt"
        params = {}
        # 读取文件
        with open(path, 'r') as f:
            contents = f.readlines()
            # 获取 kc
            kc = contents[0].split("=")[1].split("\\")[0]
            params["kc"] = float(kc)

            # 获取Wn
            Wn = contents[1].split("=")[1].split("\\")[0]
            params["wn"] = float(Wn)
            
            # 获取空气中的自由共振频率 FreqAir（kHz）
            freqAir = contents[2].split("=")[1].split("\\")[0]
            params["airFreq"] = float(freqAir)
            
            # 获取针尖的相对位置
            tipPosition = contents[3].split("=")[1].split("\\")[0]
            params["tipPosition"] = float(tipPosition)
            
            # 获取deflv修正值
            XLineWidth = contents[4].split("=")[1].split("\\")[0]
            params["XLineWidth"] = float(XLineWidth)
            
            # 获取计算模态
            n = contents[5].split("=")[1].split("\\")[0]
            params["n"] = int(n)

        return params
    
    # 读取根路径中的所有文件路径
    def readPath(self,root):
        # 1：读取原始数据名称
        listPaths = os.listdir(root)

        # 2：过滤数据文件，只保留csv文件
        dataPaths = []
        for listPath in listPaths:
            # 以任意字符开始，以.csv结束
            dataPath = re.search("[0-9]+\.txt$", listPath)
            # 如果没有csv文件，则继续，如果有csv文件，则添加到self.dataPaths中
            if dataPath == None :
                continue
            else:
                dataPaths.append(dataPath.group())
        return dataPaths
    
    """
    读取数据，并是否对数据进行旋转
    @param path : 带读取的文件的路径
    @param rotNum : 旋转次数
    """
    def readData(self,path, rotNum = 0):
        # 读取数据
        rowData = np.loadtxt(path)
        rowData = rowData / 1000 # 转化为kHz
        # 对数据进行旋转,默认不旋转
        if rotNum != 0:
            # 旋转
            rowData = np.rot90(rowData, rotNum)
        
        return rowData
    
    
    # 计算Kcontact数据
    """
    计算Kcontact
    @param data : 带计算数据
    @param airFreq : AFM探针在空气中的共振频率
    @param tipPosition : 针尖的位置参数,默认为1.0
    @param kc : AFM探针的弹性系数
    """
    def calcKcontact(self, data, airFreq, tipPosition, kc):
        # 计算波数
        knLcont = 1.8751 * np.sqrt( data / airFreq )
        knLcontR = knLcont * tipPosition
        knLcont1_R = knLcont * (1 - tipPosition)
        A = 1 + np.cos(knLcont) * np.cosh(knLcont)
        D1 = np.sin(knLcont1_R) * np.cosh(knLcont1_R) - np.cos(knLcont1_R) * np.sinh(knLcont1_R)
        D2 = 1 - np.cos(knLcontR) * np.cosh(knLcontR)
        D3 = np.sin(knLcontR) * np.cosh(knLcontR) - np.cos(knLcontR) * np.sinh(knLcontR)
        D4 = 1 + np.cos(knLcont1_R) * np.cosh(knLcont1_R)
        D = D1 * D2 - D3 * D4
        Kcont = 2 / 3 * kc * np.power(knLcontR ,3) * A / D
        return Kcont
    
    # 计算模量
    """
    计算模量数据
    """
    def calcE(self):
        pass
    
    # 画图并保存
    """
    huaKcontact图并保存
    @param data : 需要保存的数据
    @param path : 图片需要保存的路径
    @param vmin : imshow中colorbar范围的最小值
    @param vmax : imshow中colorbar范围的最大值
    """
    def plotKcontact(self, data, path, vmin=0, vmax=0):
        # 设置字体
        plt.rc('font',family='Arial')
        # 创建一个画布
        # fig = plt.figure()
        # ax = fig.add_subplot()
        fig,(ax) = plt.subplots()
        
        # 判断是否参数中有colorbar的坐标值
        if vmin==0 and vmax==0:
            pos = ax.imshow(data)
        else:
            pos = ax.imshow(data, vmin=vmin, vmax=vmax)
        fig.colorbar(pos,ax=ax)
        # 隐藏坐标轴
        plt.xticks([])
        plt.yticks([])
        plt.show()
        # 保存
        fig.savefig(path, dpi=200, bbox_inches = 'tight')

    # 保存数据
    def saveData(self, data, path):
        # 组合路径
        np.savetxt(path, data, fmt="%8.3f", delimiter=",")
    
# 实体类
class Calc:
    # 初始化
    def __init__(self):
        self.root = r"D:\OneDrive - zju.edu.cn\works\master\课题组\experience data\2023-05-03\OXFORD 6 AC200 TIP1\Freq"
        
        # 创建工具类
        self.utils = Utils()
        
        # 读取参数
        self.params = self.utils.readParams(self.root)
        
        # 读取路径
        self.paths = self.utils.readPath(self.root)
        
    def calc(self):
        for i in range(len(self.paths)):
            # 组合路径
            path = self.root + "/" + self.paths[i]
            print(path)
            # 数据保存路径
            dataPath = path[:-4] + "_Kcontact.csv"
            # 图片保存路径
            figPath = path[:-4] + "_Kcontact_fig.jpg"

            # 读取数据
            #     def readData(self,path, rotNum = 0):
            data = self.utils.readData(path,1)

            # 计算Kcontact
            #     def calcKcontact(self, data, airFreq, tipPosition=1.0， kc):
            kcont = self.utils.calcKcontact(data, self.params["airFreq"],self.params["tipPosition"],  self.params["kc"])

            # 保存数据
            #     def saveData(self, data， path):
            self.utils.saveData(kcont, dataPath)

            # 画图
            #     def plotKcontact(self, data, path):
            self.utils.plotKcontact(kcont, figPath)



####################### 画高度图，频率图，刚度图，模量图 ############################
import numpy as np
import matplotlib.pyplot as plt

def get_path(path):
    """
    对传入的路径字符串进行处理，得到 数据保存路径， 图片保存路径
    @return saveDataPath : 数据保存路径
            savePicName : 图片保存路径
    """
    # 对文件名进行分割
    pathSplit = path.split("/")
    # 获取文件名
    fileName = pathSplit[-1]
    
    # 组合数据保存名称
    saveDataName = fileName.split(".")[0] + "_handled.txt"
    saveDataPath = path.replace(fileName, saveDataName)
    
    # 组合保存图片路径
    savePicName = path.replace("txt","jpg")
    
    return [saveDataPath, savePicName]
    
def get_data(path):
    """
    读取数据
    @param path : 数据路径
    """
    # 读取原始数据，原始数据中 高度是 米 的单位
    data = np.loadtxt(path)
    return data

def H_handle(data, rotNum, savePath):
    """
    对形貌数据进行处理
    @param data : 读取的原始数据
    @param rotNum : 是否对数据进行逆时针旋转， 0 是不旋转，1 是旋转90度，
    """
    # 转换单位，m -> nm
    data = data * 10**9
    # 如果是90度扫描 将数据旋转90度
    data = np.rot90(data, rotNum)
    
    # 保存数据, 以逗号作为分隔符
    np.savetxt(savePath, data, fmt="%8.3f", delimiter=",")
    return data

def Freq_handle(data, rotNum, savePath):
    """
    对频率数据进行处理
    @param data : 原始数据
    @param rotNum : 是否对数据进行逆时针旋转， 0 是不旋转，1 是旋转90度，
    @param savePath : 保存数据的路径
    """
    # 单位转化 Hz -> kHz
    data = data / 1000
    # 如果是90度扫描 将数据旋转90度
    data = np.rot90(data, rotNum)
    
    # 保存数据, 以逗号作为分隔符
    np.savetxt(savePath, data, fmt="%8.3f", delimiter=" ")
    return data

def K_handle(data,rotNum, paramPath, savePath):
    """
    计算刚度图
    @param data : 原始的频率数据
    @param paramPath : 参数路径
    @param 
    
    """
    # 单位转化 Hz -> kHz
    data = data / 1000
    # 如果是90度扫描 将数据旋转90度
    data = np.rot90(data, rotNum)
    
    # 读取参数
    params = {}
    # 读取文件
    with open(paramPath, 'r') as f:
        contents = f.readlines()
        # 获取 kc
        kc = contents[0].split("=")[1].split("\\")[0]
        params["kc"] = float(kc)

        # 获取Wn
        Wn = contents[1].split("=")[1].split("\\")[0]
        params["wn"] = float(Wn)

        # 获取空气中的自由共振频率 FreqAir（kHz）
        freqAir = contents[2].split("=")[1].split("\\")[0]
        params["airFreq"] = float(freqAir)

        # 获取针尖的相对位置
        tipPosition = contents[3].split("=")[1].split("\\")[0]
        params["tipPosition"] = float(tipPosition)

        # 获取deflv修正值
        XLineWidth = contents[4].split("=")[1].split("\\")[0]
        params["XLineWidth"] = float(XLineWidth)

        # 获取计算模态
        n = contents[5].split("=")[1].split("\\")[0]
        params["n"] = int(n)
    
    ## 计算刚度
    knLcont = 1.8751 * np.sqrt( data / params["airFreq"] )
    knLcontR = knLcont * params["tipPosition"]
    knLcont1_R = knLcont * (1 - params["tipPosition"])
    A = 1 + np.cos(knLcont) * np.cosh(knLcont)
    D1 = np.sin(knLcont1_R) * np.cosh(knLcont1_R) - np.cos(knLcont1_R) * np.sinh(knLcont1_R)
    D2 = 1 - np.cos(knLcontR) * np.cosh(knLcontR)
    D3 = np.sin(knLcontR) * np.cosh(knLcontR) - np.cos(knLcontR) * np.sinh(knLcontR)
    D4 = 1 + np.cos(knLcont1_R) * np.cosh(knLcont1_R)
    D = D1 * D2 - D3 * D4
    Kcont = 2 / 3 * params["kc"] * np.power(knLcontR ,3) * A / D
    
    # 保存数据
    savePath = savePath.replace(".txt", "_K.txt")
    np.savetxt(savePath, Kcont, fmt="%8.3f", delimiter=" ")
    return Kcont

def E_handle(data, Ksub, savePath):
    """
    模量处理
    @params data : 刚度矩阵
    @params Ksub : 平均刚度
    """
    Eref = 75 # GPa 基地的模量
    Uref = 0.17 # 基地的泊松比

    Etip = 165 # GPa 针尖的模量
    Utip = 0.17 # 针尖的泊松比

    Usam = 0.186 # 石墨烯的泊松比

    Eref_star = (1 - np.power(Uref, 2)) / ( 1/Eref + (1-np.power(Utip,2))/Etip)
    Esam_star = Eref_star * np.power((data / Ksub), 1.5)
    Esam = (1 - np.power(Usam, 2)) / ( 1/Esam_star - (1-np.power(Utip,2))/Etip)
    
    # 保存数据
    savePath = savePath.replace(".txt", "_E.txt")
    np.savetxt(savePath, Esam, fmt="%8.3f", delimiter=" ")
    
    return Esam
    

def plot(data, savePath, cmap, vmin=0, vmax=0):
    """
    huaKcontact图并保存
    @param data : 需要保存的数据
    @param path : 图片需要保存的路径
    @param vmin : imshow中colorbar范围的最小值
    @param vmax : imshow中colorbar范围的最大值
    """
    # 设置字体
    plt.rc('font',family='Arial')
    # 创建一个画布
    fig,(ax) = plt.subplots()

    # 判断是否参数中有colorbar的坐标值
    if vmin==0 and vmax==0:
        pos = ax.imshow(data, cmap=cmap)
    else:
        pos = ax.imshow(data,cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(pos,ax=ax)
    # 隐藏坐标轴
    plt.xticks([])
    plt.yticks([])
    plt.show()
    # 保存
    fig.savefig(savePath, dpi=200, bbox_inches = 'tight')

    