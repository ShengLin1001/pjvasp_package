# created by wang xiao meng.


######################################################################################################################################
######################################################### 计算每一条曲线的刚度 ########################################################
######################################################################################################################################  
# 1 - 导入相关包
import pandas as pd
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 2 - 赫兹公式曲线拟合
def func(x, a, b):
    return a*(x**b)

# 3 - CR_FREQ 计算
class CR_FREQ:
    
    """
    ############################################################
    ##################           初始化          ###############
    ############################################################
    """
    def __init__(self, K_root, saveMaxF):
        # 成员变量
        
        ###################################### 1：根路径相关参数 ######################################
        # 根路径
        #################
        ################  请在此处更换根路径
        ###############
        self.root = K_root
        
        
        #################################### 2 ：力相关参数  ######################################
        ########## 2.1 最大力相关参数
        # saveMaxF ： 想计算并保存的最大的力的数据，该值要大于等于计算的中间力段的最大值，单位：nN
        self.saveMaxF = saveMaxF
        # sumAllData : 保存所有数据DateFrame
        self.sumAllData = pd.DataFrame()  # 保存所有的数据
        
        ########## 2.2 中间力段相关参数
        # saveMediumForceMin ： 保存的中间力段的最小值 单位：nN
        # saveMediumForceMax ： 保存的中间力段的最大值 单位：nN
        self.saveMediumForceMin = 0
        self.saveMediumForceMax = 200
        # allMediumFroceData ： 保存所有的中间力段的数据
        self.allMediumFroceData = pd.DataFrame()
        # mediumForceData：用于保存单根曲线的中间力段数据的DataFrame
#         self.mediumForceData = pd.DataFrame({"deflv":[],
#                                                "freq":[],
#                                                "f":[],
#                                                "knLcont":[],
#                                                "Kcont":[],
#                                                "1/Kcont":[],
#                                                "integration":[]})
        
        
        #################################### 3 ：图坐标轴参数  ######################################
        # 注意：坐标轴的值都设置为0，则是自适应调整坐标轴范围
        
        # axisFMin ：力的坐标轴的最小值
        # axisFMax ：力的坐标轴的最大值
        self.axisFMin = -10
        self.axisFMax = self.saveMaxF + 10
        
        # axisFreqMin ： 频率的坐标轴的最小值
        # axisFreqMax ： 频率的坐标轴的最大值
        self.axisFreqMin = 0
        self.axisFreqMax = 0
        
        # axisKcontMin ： Kcont的坐标轴的最小值
        # axisKcontMax ： Kcont的坐标轴的最大值
        self.axisKcontMin = 0
        self.axisKcontMax = 0
        
        # axisIntegrationMin ： 积分的delta的最小值
        # axisIntegrationMax ： 积分的delta的最大值
        self.axisIntegrationMin = 0
        self.axisIntegrationMax = 0
        
        # axisSaveMediumForceMin ： 中间力段的坐标轴的最小值
        # axisSaveMediumForceMax ： 中间力段的坐标轴的最大值
        self.axisSaveMediumForceMin = self.saveMediumForceMin - 10
        self.axisSaveMediumForceMax = self.saveMediumForceMax + 10

        
        #################################### 4 ：计算相关参数  ######################################
        # Kc ： 探针的弹性系数 单位：N/m
        self.kc = 0.0
        
        # Wn ： defl变化的距离和电压的转化系数 单位：nm/V
        self.wn = 0.0
        
        # n : 计算模态
        self.n = 0 # 默认第一阶模态
        
        # XLineWidth ： deflv的修正值
        self.XLineWidth = 0
        
        # freqAir ： 探针在空气中的振动频率
        self.freqAir = 0.0
        
        # tipPosition ：针尖的位置参数
        self.tipPosition = 0.0
        
        #################################### 5 ：计算相关参数  ######################################
        # params : 保存计算过程中的一些参数
        self.params = {}
        
        #################################### 6 ：读取文件路径和properties文件数据  ######################################
        self.__getPath(self.root)
        
    
    """
    ############################################################
    ##################        读取文件路径       ###############
    ############################################################
    # 注意：这里直接读取了propertise文件的数据，因为这个数据只需要读取一遍
    @Param root : 文件所在根路径
    """
    def __getPath(self, root):
        ##### 读取properties文件 ######
        # 读取参数信息
        propertiesData = self.root + "\\properties.txt"
        self.params["filePath"] = propertiesData
        with open(propertiesData, 'r') as f:
            contents = f.readlines()
            # 获取 kc
            kc = contents[0].split("=")[1].split("\\")[0]
            self.kc = float(kc)
            self.params["kc"] = self.kc

            # 获取Wn
            Wn = contents[1].split("=")[1].split("\\")[0]
            self.wn = float(Wn)
            self.params["wn"] = self.wn
            
            # 获取空气中的自由共振频率 FreqAir（kHz）
            freqAir = contents[2].split("=")[1].split("\\")[0]
            self.freqAir = float(freqAir)
            self.params["freqAir"] = self.freqAir
            
            # 获取针尖的相对位置
            tipPosition = contents[3].split("=")[1].split("\\")[0]
            self.tipPosition = float(tipPosition)
            self.params["tipPosition"] = self.tipPosition
            
            # 获取deflv修正值
            XLineWidth = contents[4].split("=")[1].split("\\")[0]
            self.XLineWidth = float(XLineWidth)
            self.params["XLineWidth"] = self.XLineWidth
            
            # 获取计算模态
            n = contents[5].split("=")[1].split("\\")[0]
            self.n = int(n)
            self.params["n"] = self.n
        
        
        
        ###### 读取文件路径 ######
        # 1：读取原始数据名称
        self.listPaths = os.listdir(self.root)

        # 2：过滤数据文件，只保留csv文件
        self.dataPaths = []
        for listPath in self.listPaths:
            # 以任意字符开始，以.csv结束
            dataPath = re.search("[0-9]+\.txt$", listPath)
#             dataPath = re.search("[0-9]+\.csv$", listPath)
            # 如果没有csv文件，则继续，如果有csv文件，则添加到self.dataPaths中
            if dataPath == None :
                continue
            else:
                self.dataPaths.append(dataPath.group())
    
    
    
    """
    ############################################################
    ##################        读取csv文件        ###############
    ############################################################
    # 注意 ：这里读取数据，只去读csv文件的数据
    @Param path : 文件路径，从__getPath()获得
    """
    def __readData(self, path):
        ###### 读取csv文件数据 #####
        rawData = pd.read_csv(path, sep="\t")
#         rawData = pd.read_csv(path, sep=",")
        # 数据过滤,
        data = rawData.iloc[:,0:2];
#         print(data)
        # 数据处理
        data["deflv"] = data.iloc[:,0] # delfv
        data["freq"] = data.iloc[:,1] / 1000;  # 频率(kHz)
        
        return data
        
    """
    ############################################################
    ##################        计算数据           ###############
    ############################################################
    """
    def calculateData(self):
        self.num = 1  # 计数
        
        ##### 开启循环计算每一条csv数据的曲线
        for path in self.dataPaths:
            # 每计算一条曲线救记录一个标志
            self.params[f"={self.num}="] = f"======================  {self.num}  ===================="
            
            # 1 - 获原始数据
            self.dataPath = self.root + "\\" + path
            self.oldData = self.__readData(self.dataPath)
            
            # 2 - 找到零点的坐标 # 
            # 2.1 找到Freq最大值及其坐标
            self.freqMax = self.oldData["freq"].max()
            self.freqMaxIndex = self.oldData["freq"].idxmax()
            self.params[f"freqMax{self.num}"] = self.freqMax
            self.params[f"freqMaxIndex"] = self.freqMaxIndex
            
            # 2.2 找到deflv的最小值
            self.deflvMin = self.oldData["deflv"].min() 
            self.deflvMinIndex = self.oldData["deflv"].idxmin()
            self.params[f"deflvMin{self.num}"] = self.deflvMin
            self.params[f"deflvMinIndex{self.num}"] = self.deflvMinIndex
            # 判断是否有deflv的修正值
            if (self.XLineWidth == 0):
                # deflv的修正值为0，则直接将deflv的最小值当作当作零点
                self.zeroPointIndex = self.deflvMinIndex
                self.params[f"zeroPointIndex{self.num}"] = self.deflvMinIndex
            else:
                # deflv 的修正值不为0
                deflvMinToLarge = self.deflvMin + self.XLineWidth / 2
                deflvMinToLargeDataFrame = pd.DataFrame()  # 存储比deflv最小的稍微大一点的数字
                deflvMinToLargeDataFrame = self.oldData[self.oldData["deflv"] < deflvMinToLarge ];
                
                self.zeroPointIndex = deflvMinToLargeDataFrame["freq"].idxmin()
                self.params[f"zeroPointIndex{self.num}"] = self.zeroPointIndex
            
            # 2.3 重新整理数据
            self.sumData = pd.DataFrame();  # 重置self.sumData,以便记录新的力曲线的数据
            self.sumData[f"deflv{self.num}"] = self.oldData["deflv"].iloc[self.freqMaxIndex:self.zeroPointIndex]
            self.sumData[f"freq{self.num}"] = self.oldData["freq"].iloc[self.freqMaxIndex:self.zeroPointIndex]
            # 重置索引
            self.sumData = self.sumData.reset_index(drop=True)
            
            # 2.4 计算相关数值
            # 2.4.1 计算f
            self.sumData[f"f{self.num}"] = self.sumData[f"deflv{self.num}"] * self.kc * self.wn
            self.sumData[f"f{self.num}"] = self.sumData[f"f{self.num}"] - self.sumData[f"f{self.num}"].iloc[-1] # 将力置零
            # 截取数据
            self.sumData = self.sumData[self.sumData[f"f{self.num}"] < self.saveMaxF]
            
            # 2.4.2 计算knLcont
            if(self.n == 1):
                # 按照第一阶模态计算
                self.sumData[f"knLcont{self.num}"] = 1.8751 * np.sqrt( self.sumData[f"freq{self.num}"] / self.freqAir )
            elif(self.n == 2):
                # 按照第二阶模态计算
                self.sumData[f"knLcont{self.num}"] = 4.6941 * np.sqrt( self.sumData[f"freq{self.num}"] / self.freqAir )
            
            # 2.4.3 计算Kcont
            dataTemp = pd.DataFrame() # 存储临时数据
            dataTemp["knLcont"] = self.sumData[f"knLcont{self.num}"]
            dataTemp["knLcontR"] = self.sumData[f"knLcont{self.num}"] * self.tipPosition
            dataTemp["knLcont1_R"] = self.sumData[f"knLcont{self.num}"] * (1 - self.tipPosition)
            dataTemp["A"] = 1 + np.cos(dataTemp["knLcont"]) * np.cosh(dataTemp["knLcont"])
            dataTemp["D1"] = np.sin(dataTemp["knLcont1_R"]) * np.cosh(dataTemp["knLcont1_R"]) - np.cos(dataTemp["knLcont1_R"]) * np.sinh(dataTemp["knLcont1_R"])
            dataTemp["D2"] = 1 - np.cos(dataTemp["knLcontR"]) * np.cosh(dataTemp["knLcontR"])
            dataTemp["D3"] = np.sin(dataTemp["knLcontR"]) * np.cosh(dataTemp["knLcontR"]) - np.cos(dataTemp["knLcontR"]) * np.sinh(dataTemp["knLcontR"])
            dataTemp["D4"] = 1 + np.cos(dataTemp["knLcont1_R"]) * np.cosh(dataTemp["knLcont1_R"])
            dataTemp["D"] = dataTemp["D1"] * dataTemp["D2"] - dataTemp["D3"] * dataTemp["D4"]
            self.sumData[f"Kcont{self.num}"] = 2 / 3 * self.kc * np.power(dataTemp["knLcontR"] ,3) * dataTemp["A"] / dataTemp["D"]
            
            # 2.4.4 计算 1 / Kcont
            self.sumData[f"1/Kcont{self.num}"] = 1 / self.sumData[f"Kcont{self.num}"]
            
            # 2.4.5 积分
            # 反转数据
            self.sumData = self.sumData[::-1]
            self.sumData = self.sumData.reset_index(drop=True)
            self.sumData[f"integration{self.num}"] = self.sumData[f"1/Kcont{self.num}"].cumsum()
            self.sumData[f"integration{self.num}"] = self.sumData[f"integration{self.num}"] * self.sumData[f"f{self.num}"].max() / len(self.sumData)
            
            # 保存csv数据计算出来的曲线数据
            self.saveData(self.sumData, 2)
            
            print("图 : " , self.num)
            print(self.dataPath)
            
            # 画全力的图片
            self.plot([self.sumData[f"f{self.num}"],self.sumData[f"f{self.num}"]],[self.sumData[f"freq{self.num}"],self.sumData[f"Kcont{self.num}"]],[("F(nN)", "Freq(kHz)"),("F(nN)","Kcont(N/m)")], 0)
            
#             if (self.saveMediumForceMax == 0 and self.saveMediumForceMin == 0):
#                 continue
#             else:
#                 # 保存中间数据
#                 self.mediumForceData = pd.DataFrame({f"deflv{self.num}":[],
#                                        f"freq{self.num}":[],
#                                        f"f{self.num}":[],
#                                        f"knLcont{self.num}":[],
#                                        f"Kcont{self.num}":[],
#                                        f"1/Kcont{self.num}":[],
#                                        f"integration{self.num}":[]})
                
#                 self.mediumForceData = self.sumData[(self.sumData[f"f{self.num}"] >= self.saveMediumForceMin) & (self.sumData[f"f{self.num}"] <= self.saveMediumForceMax) ]
                
#                 # 将中间力段参数保存到总数据表中
#                 self.allMediumFroceData = pd.concat([self.allMediumFroceData, self.mediumForceData], axis=1)
                
                
#                 # 画中间力段的图
#                 self.plot([self.mediumForceData[f"f{self.num}"],self.mediumForceData[f"f{self.num}"]],[self.mediumForceData[f"freq{self.num}"],self.mediumForceData[f"Kcont{self.num}"]],[("MediumForce(N/m)", "Freq(kHz)"),("MediumForce(N/m)","Kcont(N/m)")], 1)
                   
#                 # 保存数据
#                 self.saveData(self.mediumForceData, 3)
            
            # 4 - 将数据保存到sumAllData中
            self.sumAllData = pd.concat([self.sumAllData, self.sumData], axis=1)
                
            # 更新序号
            self.num = self.num + 1
        
#         # 保存中间力段参数
#         self.saveData(self.allMediumFroceData, 1)
        
        # 保存所有数据
        self.saveData(self.sumAllData, 0)
        
        
    """
    ############################################################
    ##################      设置图形的坐标轴     ###############
    ############################################################
    @param ax ： 需要设置图层
    @param axisNum : 需要设置的坐标轴， 0 x轴， 1 y轴
    @param axisMin : 坐标轴的最小值
    @param axisMax : 坐标轴的最大值
    @param axisLabel : 坐标轴标签
    """
    def __setLim(self, ax, axisNum, axisMin, axisMax, axisLabel):
        if (axisMin == 0 and axisMax == 0):
            # 只设置坐标轴的标签
            if axisNum == 0:
                ax.set_xlabel(axisLabel)
            elif axisNum == 1:
                ax.set_ylabel(axisLabel)
        else:
            # 设置坐标轴标签和范围
            if (axisNum == 0): # 设置x轴
                ax.set_xlim(axisMin, axisMax)
                ax.set_xlabel(axisLabel)
            elif (axisNum == 1): #设置y轴
                ax.set_ylim(axisMin, axisMax)
                ax.set_ylabel(axisLabel)
    
    """
    ############################################################
    ##################        设置坐标轴         ###############
    ############################################################
    @Param ax : 需要设置坐标轴的图
    
    @Param axisNum : 坐标轴代号
                 0 : x轴
                 1 : y轴
    @Param axisLabel : 坐标轴标签
                "F(nN)" : 
                "Freq(kHz)" : 
                "Kcont(N/m)" : 
                "delta(nm)" :
    """
    def __setAxis(self,ax, axisNum, axisLabel):
        # 设置代表力的坐标轴
        if axisLabel == "F(nN)":
            self.__setLim(ax, axisNum, self.axisFMin, self.axisFMax, axisLabel)
        
        # 设置中间力段
        elif axisLabel == "MediumForce(N/m)":
            self.__setLim(ax, axisNum, self.axisSaveMediumForceMin, self.axisSaveMediumForceMax, "MediumForce(N/m)")
        
        # 设置代表频率的坐标轴
        elif axisLabel == "Frequency(kHz)":
            self.__setLim(ax, axisNum, self.axisFreqMin, self.axisFreqMax, axisLabel)
        
        # 设置代表频率的坐标轴
        elif axisLabel == "Kcontact(N/m)":
            self.__setLim(ax, axisNum, self.axisKcontMin, self.axisKcontMax, axisLabel)
            
        # 设置代表Kcont的坐标轴
        elif axisLabel == "delta(nm)":
            self.__setLim(ax, axisNum, self.axisIntegrationMin, self.axisIntegrationMax, axisLabel)
    
    
    """
    ############################################################
    ##################          herzt拟合        ###############
    ############################################################
    @ x : 需要拟合曲线的x轴数据
    @ y : 需要拟合曲线的y轴数据
    """
    def __fit(self, x, y):
        # 参数： maxfev，最大拟合次数
        popt, pcov = curve_fit(func, x, y)
        
        # 根据拟合出来的参数，获取拟合曲线数据
        x_pred = pd.DataFrame(np.arange(0,x.max(),0.005))
        y_pred = func(x_pred, popt[0], popt[1])
        # 转化为 np.array() 不然会出错
        x_pred = np.array(x_pred)
        y_pred = np.array(y_pred)
        
        # 计算 R_squared
        mean = np.mean(y)  # 1.计算均值
        ss_tot = np.sum((y - mean) ** 2)  # 2.total sum of squares
        ss_res = np.sum((y - func(x, popt[0], popt[1])) ** 2)  # 3.residual sum of squares
        r_squared = 1 - (ss_res / ss_tot)  # 4.r squared
        
        # 保存拟合参数
        self.params[f"a{self.num}"] = popt[0]  # 保存a
        self.params[f"b{self.num}"] = popt[1]  # 保存b
        self.params[f"r2{self.num}"] = r_squared # 保存r2
        
        return x_pred, y_pred, popt, r_squared
                
    
    """
    ############################################################
    ##################             画图          ###############
    ############################################################
    @Params *params : 给参数有三部分组成
                第一部分 ：[x轴数据]
                第二部分 ：[y轴数据]
                第三部分 ：[图的坐标]，其中每个都是一个元组，每个元组包含两个数据，第一个数据是x轴的标签，第二个数据是y轴的标签
                            例如：[("F(nN)","Freq(kHz)"),("F(nN)","Kcont(N/m)"),("delta(nm)","F(nN)")]
                第四部分 ：图片所画的图所使用的原始数据
                       0 ：全力段的图
                       1 ：中间力段的图片
    """
    def plot(self, *params):
        # 设置字体
        plt.rc('font',family='Arial')
        
        # 提取参数
        figNum = len(params[0])  # 需要画的图的数量
        
        # 画图
        fig = plt.figure(figsize=(figNum*5,4), dpi=200)  # 创建一个画布
        
        for i in range(figNum):
            ax = fig.add_subplot(1, figNum, i+1) 
            ax.scatter(params[0][i], params[1][i],c="black", s=2)
            
            # 设置标题
            title = params[2][i][1] + "-" + params[2][i][0]
            ax.set_title(title)
            
            # 设置坐标轴范围及标签
            # def __setAxis(self, ax, axisNum, axisLabel):
            self.__setAxis(ax, 0, params[2][i][0])  # 设置x轴坐标范围及标签
            self.__setAxis(ax, 1, params[2][i][1])
            
        # 展示图片
        plt.show()
        
        # 判断图片类型
        if params[3] == 0:
            # 全段力图片
#             figPath = self.dataPath.split(".")[0] + "_fig.jpg"
            figPath = self.dataPath[:-4] + "_fig.jpg"
            fig.savefig(figPath, dpi=200, bbox_inches = 'tight')
        elif params[3] == 1:
            # 中间力段图片
#             figPath = self.dataPath.split(".")[0] + "_MeduimForceFig.jpg"
            figPath = self.dataPath[:-4] + "_MeduimForceFig.jpg"
            fig.savefig(figPath, dpi=200, bbox_inches = 'tight')
    
    
    """
    ############################################################
    ##################          画汇总图         ###############
    ############################################################
    @param data : 需要画图的数据
    @param figFlag : 保存的图片
                0 ：全段力的图片
                1 ：中间力段的图片
    """
    def plotAll(self, data, figFlag, Fad):

        
        # 设置字体
        plt.rc('font',family='Arial')
        
        
        # 画图的数量
        figNum = 2;
        
        fig = plt.figure(figsize=(figNum*4,3), dpi=200) # 创建画布
        ax1 = fig.add_subplot(1,figNum,1)
        ax2 = fig.add_subplot(1,figNum,2)
#         ax3 = fig.add_subplot(1,3,3)
        colors = ["black","red","blue","purple","springgreen", 
                  "orange", "bisque", "forestgreen","darkolivegreen", "pink" ,
                  "evergreen", "bright turquoise","umber","denim blue" ]
        
        ### 设置坐标轴范围
        if Fad > 0:  # 由Fad
            self.axisFMin = self.axisFMin - Fad
            self.axisFMax = self.axisFMax - Fad
            
    
    
        for i in range(1,self.num):
            ########################################### Ktot - F 
            
            ### 对data进行处理，减去Fad
            data[f"f{i}"] = data[f"f{i}"] - Fad
            
            
            ax1.scatter(data[f"f{i}"], data[f"freq{i}"],s=2, c=colors[i-1], label=f"{i}")
            if figFlag == 0:
                # 全力段图
                self.__setAxis(ax1, 0, "F(nN)") # x轴
            elif figFlag == 1:
                # 中间力段图
                self.__setAxis(ax1, 0, "MediumForce(N/m)") # x轴
            self.__setAxis(ax1, 1, "Frequency(kHz)") # y轴
            ax1.set_title("Freq - F")
            ax1.legend(fontsize=10,frameon=False)
            
            ax2.scatter(data[f"f{i}"], data[f"Kcont{i}"],s=2, c=colors[i-1], label=f"{i}")
            if figFlag == 0:
                # 全力段图
                self.__setAxis(ax2, 0, "F(nN)") # x轴
            elif figFlag == 1:
                # 中间力段图
                self.__setAxis(ax2, 0, "MediumForce(N/m)") # x轴
            self.__setAxis(ax2, 1, "Kcontact(N/m)") # y轴
            ax2.set_title("Kcont - F")
            ax2.legend(fontsize=10,frameon=False)
            
#             ax3.scatter(data[f"integration{i}"], data[f"f{i}"],s=2, c=colors[i-1], label=f"{i}")
#             if figFlag == 0:
#                 # 全力段图
#                 self.__setAxis(ax3, 0, "delta(nm)") # x轴
#             elif figFlag == 1:
#                 # 中间力段图
#                 self.__setAxis(ax3, 0, "delta(nm)") # x轴
#             self.__setAxis(ax3, 1, "F(nN)") # y轴
#             ax3.set_title("F - integration")
#             ax3.legend(fontsize=10)
            
        plt.rc('font',family='Arial')
        plt.show()
        
        # path[:-4]
        # 判断图片的数据来源
        if figFlag == 0:
            # 全段力的图片
            figPath = self.root + "/sumAllFig.jpg"
            fig.savefig(figPath, dpi=200, bbox_inches = 'tight')
            
        elif figFlag == 1:
            # 中间力段的图片
            figPath = self.root + "//allMediumForceFig.jpg"
            fig.savefig(figPath, dpi=200, bbox_inches = 'tight')
            
    
    """
    ############################################################
    ##################         保存计算数据      ###############
    ############################################################
    @Param data : 需要保存的数据
    @Param num : 要保存的数据的类型
             0 : 保存计算出来的每一条曲线的汇总数据
             1 : 保存中间力段的汇总数据
             2 : 保存每一个csv计算出来的数据
             3 : 保存每一个csv计算出来的中间力段的数据
    """
    def saveData(self, data, num):
        # 0 : 保存计算出来的每一条曲线的汇总数据
        if num == 0:
            allDataPath =self.root + "\AllData.csv"
            data.to_csv(allDataPath, encoding="utf_8_sig")
        
        # 1 : 保存中间力段的汇总数据
        elif num == 1:
            allMediumForceDataPath =self.root + "\\allMediumForceDataPath.csv"
            data.to_csv(allMediumForceDataPath, encoding="utf_8_sig")
        
        # 2 : 保存每一个csv计算出来的数据
        elif num == 2:
            # path[:-4]
#             csvDataPath =self.dataPath.split(".")[0]  + "_data.csv"
            csvDataPath =self.dataPath[:-4]  + "_data.csv"
            data.to_csv(csvDataPath, encoding="utf_8_sig")
        
        # 3 : 保存每一个csv计算出来的中间力段的数据
        elif num == 3:
            mediumForceDataPath =self.dataPath[:-4]  + "_MediumFirceData.csv"
            data.to_csv(mediumForceDataPath, encoding="utf_8_sig")
            
    
    """
    ############################################################
    ##################           保存参数        ###############
    ############################################################
    """
    def saveParams(self):
        paramsPath = self.root + "/params.txt"
        with open(paramsPath, 'w') as f:
            for k ,v in self.params.items():
                f.write(k)
                f.write(" : ")
                f.write(str(v))
                f.write("\n") 

                
######################################################################################################################################
################################################################### 计算K均值 ########################################################
######################################################################################################################################  
# 导入相关包
import pandas as pd
import os


# 读取路径
def get_average_K_Path(dataRoot):
    # 读取相关数据
    fileList = os.listdir(dataRoot)
    # 创建字典
    pathDict = {} # 用于存储路径
    # 筛选文件夹
    for f in fileList:
        if "." not in f:
            # 文件路径
            pathDict[f] = []  # 待处理文件路径
            dataPath  = dataRoot + "\\" + f  # 文件路径
            pathList = os.listdir(dataPath)
            for path in pathList:
                # 判断是否是需处理文件
                if path.endswith("_data.csv"):
                    pathDict[f].append(path)
    return pathDict

# 处理数据
def draw_point(sumData, data, fluctuation, file):
    # 创建每一个曲线文件的列
    sumData[file] = 0
    for i in range(len(sumData["F"])):
        # 注意 这里不能是 and 只能是 &
        averageFluctuationValue = data.iloc[:,5][((sumData["F"].loc[i]-fluctuation)<data.iloc[:,3]) & (data.iloc[:,3]<(sumData["F"].loc[i]+fluctuation))].mean()
        sumData[file].loc[i] = averageFluctuationValue

# 处理数据
def handle_data(dataRoot, pathDict, startPoint, endPoint, Fgap, fluctuation, curveNum):
    # 路径列表
    paths = {}
    legends = []
    # 创建汇总数据表
#     sumData = pd.DataFrame({"F":[i for i in range(startPoint, endPoint, Fgap)]})
    sumData = pd.DataFrame({"F":[i for i in np.arange(startPoint, endPoint, Fgap)]})
    # 遍历字典
    for folder, files in pathDict.items():
        for file in files:
            curData = pd.read_csv(dataRoot + "\\" + folder + "\\" + file)
            draw_point(sumData, curData, fluctuation, file)
        
        # 计算average_K
        sumData["average_K"] = 0
        for i in range(len(sumData)):
            sumData["average_K"].loc[i] = sumDa111111ta.iloc[i,1:(curveNum+1)].mean()
        
        # 计算完一个文件夹中的数据，保存，重置
        savePath = dataRoot + "\\" + folder + f"\\{folder}_average_K.csv"
        paths[folder] = savePath
        legends.append(folder)
        sumData.to_csv(savePath, encoding="utf-8")
#         sumData = pd.DataFrame({"F":[i for i in range(startPoint, endPoint, Fgap)]})
        sumData = pd.DataFrame({"F":[i for i in np.arange(startPoint, endPoint, Fgap)]})
    
#     paths.reverse()
    legends.reverse()
    
    return paths, legends
        
        
        
######################################################################################################################################
#################################################################### 刚度汇总 ########################################################
######################################################################################################################################                
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# def get_path(root, legends):
#     """拼接数据路径"""
#     paths = {}
#     # 获取根路径下的所有文件
#     allFiles = os.listdir(root)
#     # 遍历所有文件，提取出文件夹
#     for file in allFiles:
#         file = root + "\\" + file
#         # 判断是否是文件夹
#         if os.path.isdir(file):
#             # 提取文件加中所有的文件
#             currentForderFiles = os.listdir(file)
#             # 遍历该文件夹
#             for i in range(len(legends)):
#                 for cFile in currentForderFiles:
#                     curAllFile = file + "\\" + cFile
#                     # 判断该文件夹中是否有需要的文件
#                     if cFile.endswith("_average_K.csv"):
#                         paths[legends[i]] = curAllFile
#     for k, v in paths.items():
#         print(k, " >>> ", v)
#     return paths
            
    
def get_data_sum_K(path, Fad):
    """读取数据"""
    newData = pd.DataFrame() # 存放新数据
    # 拼接数据
    data = pd.read_csv(path)
    # 数据过滤
    # 去掉Fad
#     newData["F"] = data.iloc[:,3]-Fad  # 力， nN
#     newData["Kcontact"] = data.iloc[:,5] # 刚度 N/m
    newData["F"] = data.iloc[:,1]-Fad  # 力， nN
    newData["Kcontact"] = data.iloc[:,-1] # 刚度 N/m
    return newData

def plot(datas,legends,root,flag):
    """
    画图
    @params flag : 是否减去Fad, 0, 不减去Fad， 1，减去Fad
    """
    # 设置字体
    plt.rc('font',family='Arial')
    # 创建画布
    fig = plt.figure(figsize=(5,4), dpi=200) # 创建画布
    # 刚度汇总
    ax = fig.add_subplot()
    # 设置颜色
    colors = ["red", "#1a6840","#00FF7F", "#168888", "#11a9cf"]
    
    # 添加数据
    for i in range(len(datas)):
        ax.scatter(datas[i].iloc[:,0], datas[i].iloc[:,1], c=colors[i], label=legends[i], s=2)
        
    # 设置label
    ax.set_xlabel("F(nN)")
    ax.set_ylabel("Kcontact(N/m)")
    ax.set_title("F - K")
    ax.legend(fontsize=10)
    # 去掉legend的边框
    ax.legend(frameon=False)
    
    plt.show()
    
    # 保存
    if flag == 0:
        # 全段力的图片
        figPath = root + "\\K_Fad.jpg"
        fig.savefig(figPath, dpi=200, bbox_inches = 'tight')

    elif flag == 1:
        # 中间力段的图片
        figPath = root + "\\K_no_Fad.jpg"
        fig.savefig(figPath, dpi=200, bbox_inches = 'tight')


######################################################################################################################################
#################################################################### 力曲线 ########################################################
######################################################################################################################################   

# 力曲线
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import scipy
from scipy.integrate import simps # 用于计算积分

# def get_paths(root, nums, legends):
#     """拼接数据路径"""
#     paths = {}
#     # 获取根路径下的所有文件
#     allFiles = os.listdir(root)
#     # 遍历所有文件，提取出文件夹
#     for file in allFiles:
#         file = root + "\\" + file
#         # 判断是否是文件夹
#         if os.path.isdir(file):
#             # 提取文件加中所有的文件
#             currentForderFiles = os.listdir(file)
#             # 遍历该文件夹
#             for i in range(len(nums)):
#                 for cFile in currentForderFiles:
#                     cFile = file + "\\" + cFile
#                     # 判断该文件夹中是否有需要的文件
#                     if (str(nums[i]) + "_data") in cFile:
#                         paths[legends[i]] = cFile
#     for k, v in paths.items():
#         print(k, " >>> ", v)
#     return paths

def get_data(path):
    """获取数据"""
    newData = pd.DataFrame() # 存放新数据
    # 拼接数据
    data = pd.read_csv(path)
    # 数据过滤
    # 去掉Fad
    newData["F"] = data.iloc[:,1]  # 力， nN
    newData["Kcontact"] = data.iloc[:,-1] # 刚度 N/m
    return newData


def calc_data(data, Fad, Fmax):
    """处理数据"""
    # 减去粘附力
    data["F"] = data["F"] - Fad # 减去Fad
    # 截取最大的力
#     newData = data
    # 求1/K
    data["1/K"] = 1 / data["Kcontact"]
    data["delta"] = 0
    # 积分
#     data["1/K_cumsum"] = data["1/K"].cumsum()
#     data["delta"] = (data["F"].max() - data["F"].min()) * data["1/K_cumsum"] / len(data)
    for i in range(len(data["F"])): # 计算梯形的面积，由于是累加，所以是切片"i+1"
        data["delta"][i] = scipy.integrate.trapz(data["1/K"][:i + 1], data["F"][:i + 1])
#         integrals.append(scipy.integrate.trapz(y[:i + 1], x[:i + 1]))

    dataFmaxBefore = data[data["F"] <= Fmax]  # Fmax 之前的数值
    dataFmaxAfter = data[data["F"] >= Fmax]  # Fmax 之后的数值
    return dataFmaxBefore,dataFmaxAfter

def saveData(root,data, legend, flag, Fmax):
    """
    保存数据
    @param flag : before 表示 Fmax 之前的数值
                  after  表示之后的数值
    """
    # 拼接路径
    savaPath = root + "\\FC"
    # 判断文件夹是否存在
    if not os.path.exists(savaPath):
        # 如果不存在则创建
        os.makedirs(savaPath)
    # 保存文件
    # 判断flag
    savaPath = savaPath + f"\\FC_{legend}_{Fmax}_{flag}.csv"
    data.to_csv(savaPath, encoding="utf_8_sig")
#     print(f">>> {num} 数据保存成功")
    

def plot_fc(datas,legends,root, flag, Fmax):
    """
    画图
    @param flag : before 表示 Fmax 之前的数值
                  after  表示之后的数值
    """
    
    # 设置字体
    plt.rc('font',family='Arial')
    # 创建画布 单个数据图
    fig = plt.figure(figsize=(5,4), dpi=200) # 创建画布
    # 刚度汇总
    ax = fig.add_subplot()
    
    # 创建画布，汇总数据图
    figSum = plt.figure(figsize=(5,4), dpi=200)
    axSum = figSum.add_subplot()
    
    # 设置颜色
    colors = ["red", "#1a6840","#00FF7F", "#00D2FC", "#11a9cf"]
    
    # 添加数据
    for i in range(len(datas)):
        ax.cla()
        ax.scatter(datas[i].iloc[:,-1], datas[i].iloc[:,0], c="red", label=legends[i], s=0.2)
        # 设置label
        ax.set_xlabel("ΔZ(nm)")
        ax.set_ylabel("F(nN)")
        ax.set_title("F - K")
        ax.legend(fontsize=10)
        # 去掉legend的边框
        ax.legend(frameon=False)
        
        axSum.scatter(datas[i].iloc[:,-1], datas[i].iloc[:,0], c=colors[i], label=legends[i], s=0.2)
        
        # 判断flag
        figPath = root + f"\\FC\\FC_{legends[i]}_{Fmax}_{flag}.jpg"
        fig.savefig(figPath, dpi=200, bbox_inches = 'tight')
    
    # 设置label
    axSum.set_xlabel("ΔZ(nm)")
    axSum.set_ylabel("F(nN)")
    axSum.set_title("F - K")
    axSum.legend(fontsize=10)
    # 去掉legend的边框
    axSum.legend(frameon=False)
    
    plt.show()
    
    # 保存汇总图
    figPath = root + f"\\FC\\FC_sum_{Fmax}_{flag}.jpg"
    figSum.savefig(figPath, dpi=200, bbox_inches = 'tight')