import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy import stats
from scipy.optimize import curve_fit
import os
import re

from PIL import Image
import matplotlib.image as mpimg

tickSize=25
sns.set_style("ticks", {"xtick.major.size": tickSize, "ytick.major.size": tickSize})

figure_norm = 12 #convert to 17.1cm
figure_len_factor=4/3

figure_factor=17.1/8.3#ratio for sm journal

Fontsize_Title=45
Fontsize_Sub = 40
Fontsize_Leg = 30

dotSize=5
lineWidth=3

colorChoice = sns.xkcd_rgb["black"]
def exp_shift(x, a, b, c, d, e):
     return a * np.exp(b*(x)**e + c) - d
 
def exp(x, a, b, c):
     return a * np.exp(b*(x)) - c

def fun_line(x1,y1,z1,x2,y2,z2,t):
    return [float(x1 + t*(x2-x1)),float(y1 + t*(y2-y1)), float(z1 + t*(z2-z1))];

def poly4(x, a, b, c, d,e):
    return (a*x*x*x*x + b*x*x*x + c * x * x + d * x + e)


def poly3(x, a, b, c, d):
    return (a*x*x*x + b * x * x + c * x + d)

def poly2(x, a, b, c):
    return a * x * x + b * x + c

def alpha_law(x, a, b, c):
    return a * (x**b) + c

class DataForceStrain: 

    
    def __init__(self, direct_, filename_):
        self.direct = direct_ 
        self.filename = filename_
        self.force=[]
        self.strain=[]
        self.minY=[]
        self.maxY=[]
        self.minX=[]
        self.maxX=[]
        self.MPascals=[]
        
        
class DataFile:  
    FiberArea = np.pi*(0.05**2)#assumes 50nm radius or 0.05micron radius
    def __init__(self, direct_,file_):
        self.direct = direct_
        self.file=file_
        self.netWorkStrain=0
        self.force_upper=0
        self.force_lower=0
        self.originalNodeCount=0
        self.nodeCountDiscretize=0
        self.originalEdgeCount=0
        self.numFibersInMiddleSlice=0
        self.areaCoveredByFibers = 0
        self.minX = 0.0
        self.maxX = 0.0
        self.minY = 0.0
        self.maxY = 0.0
        self.minZ = 0.0
        self.maxZ = 0.0
        self.MPascals_xy=0.0
        self.MPascals_xz=0.0
        self.central_x_zone_buffer=3
        self.central_y_zone_buffer=3
        self.central_z_zone_buffer=500
        
        self.is_node_in_central_zone = []
        
        self.xPos=[]
        self.yPos=[]
        self.zPos=[]          
        self.originEdgeLeft=[]
        self.originEdgeRight=[]
        self.addedEdgeLeft=[]
        self.addedEdgeRight=[]
        self.originalEdgeStrain=[] 
        self.originalEdgeMidPoint=[]
        self.originalEdgeMidPointScaledForStrain=[]
        self.originalEdgeAlignment=[]
        self.originalNodeForce=[]
        
        
        self.addedEdgeStrain=[]
        self.bindSitesPerNode=[]
        self.curvaturePerEdge=[]
         
        self.curvaturePerEdgeMiddle=[]
        self.curvatureXPointPerEdgeMiddle=[]
        self.curvatureYPointPerEdgeMiddle=[]
        
        self.noDivisionEdgesLeft=[] 
        self.noDivisionEdgesRight=[] 
        self.noDivisionEdgesStrain=[]
        self.noDivisionEdgesForce=[] 
        self.noDivisionEdgesAlignment=[]
        self.noDivisionEdgesBindSites=[]
        self.noDivisionEdgesUniqueBindSites=[]
        self.noDivisionEdgesAverageForce=[]#averaged over all nodes on the edge
        self.noDivisionEdgesMidPoint=[]
        self.noDivisionEdgesMidPointScaledForStrain=[]


# =============================================================================
# #take the datafile class and use it's file/direct to modify it.
# def parserForceStrain(DataForceStrain):
#     direct = DataForceStrain.direct
#     filename=DataForceStrain.filename
#     with open(os.path.join(direct, filename), "r") as f:
#         
#         for line in f: 
#             string_split = re.split("[=|\n]",line)
#             for elem in string_split:
#                 if (elem == 'network_strain ' or elem == 'network_strain'):
#                     DataForceStrain.strain.append(float(string_split[1]) )
#                 if (elem == 'force_upper ' or elem == 'force_upper'):
#                     DataForceStrain.force.append(float(string_split[1]) )
#                 if (elem == 'maxY ' or elem == 'maxY'):
#                     DataForceStrain.maxY.append(float(string_split[1]) )
#                 if (elem == 'minY ' or elem == 'minY'):
#                     DataForceStrain.minY.append(float(string_split[1]) )
#                 if (elem == 'maxX ' or elem == 'maxX'):
#                     DataForceStrain.maxX.append(float(string_split[1]) )
#                 if (elem == 'minX ' or elem == 'minX'):
#                     DataForceStrain.minX.append(float(string_split[1]) )
#         #now that all the files have been read in, 
#         for i in range(0,len(DataForceStrain.force)):
#             area = (float(DataForceStrain.maxY[i]) - float(DataForceStrain.minY[i])) *  (float(DataForceStrain.maxX[i]) - float(DataForceStrain.minX[i]))
#             DataForceStrain.MPascals.append((10**-6)*(10**3)*DataForceStrain.force[i] / area)
#             #multiply by 100 since units are in nN and microns
# =============================================================================

def parserForceStrain(DataFile):
    direct = DataFile.direct
    file = DataFile.file
    with open(os.path.join(direct, file), "r") as f:
        for line in f:
            string_split = line.split()
            for elem in string_split:
                if elem == 'network_strain':
                    DataFile.netWorkStrain = float(string_split[1])
                if (elem == 'force_upper ' or elem == 'force_upper'):
                    DataFile.force_upper = float(string_split[1])
                
                if (elem == 'force_lower ' or elem == 'force_lower'):
                    DataFile.force_lower = float(string_split[1])
                    
                if elem == 'original_node_count':
                    DataFile.originalNodeCount = float(string_split[1])
                if elem == 'node_count_discretize':
                    DataFile.nodeCountDiscretize = float(string_split[1])
                if elem == 'original_edge_count':
                    DataFile.originalEdgeCount = float(string_split[1])
                    
                if elem == 'node':
                    try:
                        x=float(string_split[1])
                        y=float(string_split[2])
                        z=float(string_split[3])
                        
                        DataFile.xPos.append(x)
                        DataFile.yPos.append(y)
                        DataFile.zPos.append(z)
                        DataFile.is_node_in_central_zone.append(False)
                    except ValueError:
                        pass 
                                
        DataFile.maxX = np.max(DataFile.xPos)
        DataFile.minX = np.min(DataFile.xPos)
        DataFile.maxY = np.max(DataFile.yPos)
        DataFile.minY = np.min(DataFile.yPos)
        DataFile.maxZ = np.max(DataFile.zPos)
        DataFile.minZ = np.min(DataFile.zPos)
        
        area_xy = (DataFile.maxX - DataFile.minX) * (DataFile.maxY - DataFile.minY)
        area_xz = (DataFile.maxX - DataFile.minX) * (DataFile.maxZ - DataFile.minZ)
            
        DataFile.MPascals_xy = ((10**-6)*(10**3)*( DataFile.force_upper + DataFile.force_lower)) / area_xy
        DataFile.MPascals_xz = ((10**-6)*(10**3)*( DataFile.force_upper + DataFile.force_lower)) / area_xz
            #multiply by 100 since units are in nN and microns



def parser(DataFile):
    direct = DataFile.direct
    file = DataFile.file
    b1X = 0
    b2X = 0
    b1Y = 0
    b2Y = 0
    with open(os.path.join(direct, file), "r") as f:
        for line in f:
            string_split = line.split()
            for elem in string_split:
                if elem == 'network_strain':
                    DataFile.netWorkStrain = float(string_split[1])
                if (elem == 'force_upper ' or elem == 'force_upper'):
                    DataFile.force_upper = float(string_split[1])
                
                if (elem == 'force_lower ' or elem == 'force_lower'):
                    DataFile.force_lower = float(string_split[1])
                    
                if elem == 'original_node_count':
                    DataFile.originalNodeCount = float(string_split[1])
                if elem == 'node_count_discretize':
                    DataFile.nodeCountDiscretize = float(string_split[1])
                if elem == 'original_edge_count':
                    DataFile.originalEdgeCount = float(string_split[1])
                   
                                
                if (elem == 'maxX ' or elem == 'maxX'):
                    b2X = (float(string_split[1]) )
                    
                if (elem == 'maxY ' or elem == 'maxY'):
                    b2Y = (float(string_split[1]) )
                    
                if (elem == 'minX ' or elem == 'minX'):
                    b1X = (float(string_split[1]) )
                     
                if (elem == 'minY' or elem == 'minY'):
                    b1Y = (float(string_split[1]) )
                    
                    
                if elem == 'node':
                    try:
                        x=float(string_split[1])
                        y=float(string_split[2])
                        z=float(string_split[3])
                        
                        DataFile.xPos.append(x)
                        DataFile.yPos.append(y)
                        DataFile.zPos.append(z)
                        DataFile.is_node_in_central_zone.append(False)
                    except ValueError:
                        pass            
                
                if elem == 'force_on_node':
                    force = float(string_split[1])
                    DataFile.originalNodeForce.append(force)
                    
                if elem == 'original_edge_discretized': 
                    eL = float(string_split[1])
                    eR = float(string_split[2])
                    DataFile.originEdgeLeft.append(eL)
                    DataFile.originEdgeRight.append(eR)
                    
                if elem == 'added_edge': 
                    eL = float(string_split[1])
                    eR = float(string_split[2])
                    DataFile.addedEdgeLeft.append(eL)
                    DataFile.addedEdgeRight.append(eR)
                    
                if elem == 'original_edge_strain':
                    DataFile.originalEdgeStrain.append(float(string_split[1]))
                    
                if elem == 'original_edge_alignment':
                    DataFile.originalEdgeAlignment.append(2 * float(string_split[1]) - 1)
                    
                if elem == 'added_edge_strain':
                    DataFile.addedEdgeStrain.append(float(string_split[1]))
                
                if elem == 'bind_sites_per_node':
                    DataFile.bindSitesPerNode.append(float(string_split[1]))
                
        b1 = np.maximum(b1X,b1Y)#min
        b2 = np.minimum(b2X,b2Y)#max       
        #now take original fibers and generate 
                        
        DataFile.maxX = np.max(DataFile.xPos)
        DataFile.minX = np.min(DataFile.xPos)
        DataFile.maxY = np.max(DataFile.yPos)
        DataFile.minY = np.min(DataFile.yPos)
        DataFile.maxZ = np.max(DataFile.zPos)
        DataFile.minZ = np.min(DataFile.zPos)
        
        buffer_x_left = (DataFile.maxX + DataFile.minX)/2.0 - DataFile.central_x_zone_buffer
        buffer_x_right = (DataFile.maxX + DataFile.minX)/2.0 + DataFile.central_x_zone_buffer
        buffer_y_left = (DataFile.maxY + DataFile.minY)/2.0 - DataFile.central_y_zone_buffer
        buffer_y_right = (DataFile.maxY + DataFile.minY)/2.0 + DataFile.central_y_zone_buffer
        buffer_z_left = (DataFile.maxZ + DataFile.minZ)/2.0 - DataFile.central_z_zone_buffer
        buffer_z_right = (DataFile.maxZ + DataFile.minZ)/2.0 + DataFile.central_z_zone_buffer
        print("bufferxl")
        print(buffer_x_left)        
        print("bufferxr")
        print(buffer_x_right)
        trueEdgeLeft = DataFile.originalNodeCount
        trueEdgeRight = DataFile.originalNodeCount
        subEdgeCounter=0
        undividedEdgeCounter=0
        undividedEdgeStrainTemp = 0.0
        undividedEdgeForceTemp = 0.0
        bindSitesPerEdge = 0
        aveAlignment = 0
        aveMidPoint = 0
        aveMidPointScaled = 0
        distcurvetemp=0.0;
        edgeCurveWasInMiddle = False
        xMiddle = 0.0
        yMiddle = 0.0
        for edge in range(0,len(DataFile.originEdgeLeft)):
            
            eL = int(DataFile.originEdgeLeft[edge])
            eR = int(DataFile.originEdgeRight[edge])
            if ((eL < DataFile.nodeCountDiscretize) and ( eR < DataFile.nodeCountDiscretize)):
                zPosL = DataFile.zPos[eL]
                zPosR = DataFile.zPos[eR]
                yPosL = DataFile.yPos[eL]
                yPosR = DataFile.yPos[eR]
                xPosL = DataFile.xPos[eL]
                xPosR = DataFile.xPos[eR]
                
                #fill central volume with nodes
                if ((zPosL > buffer_z_left) and (zPosL < buffer_z_right) and
                    (yPosL > buffer_y_left) and (yPosL < buffer_y_right) and 
                    (xPosL > buffer_x_left) and (xPosL < buffer_x_right)):
                    #print(eR)
                    #then left point is in central zone
                    DataFile.is_node_in_central_zone[eL] = True #default is false
                    
                if ((zPosR > buffer_z_left) and (zPosR < buffer_z_right) and
                    (yPosR > buffer_y_left) and (yPosR < buffer_y_right) and 
                    (xPosR > buffer_x_left) and (xPosR < buffer_x_right)):
                    #then right point is in central zone
                    #print(eL)
                    DataFile.is_node_in_central_zone[eR] = True
                    
                    
                DataFile.originalEdgeMidPoint.append(( zPosL+ zPosR)/2.0 )  
                
                #now scale midpoint to strain
                edgeScaled=DataFile.netWorkStrain*(DataFile.originalEdgeMidPoint[edge])/((DataFile.maxZ-DataFile.minZ))
                DataFile.originalEdgeMidPointScaledForStrain.append(edgeScaled)
                subEdgeCounter+=1

                #midpoint slice fiber count
                midpoint = DataFile.netWorkStrain/2
                

                if ((zPosL > midpoint) & (zPosR < midpoint)):
                    t_ = (midpoint - zPosL) / (zPosR - zPosL);
                    [xTemp,yTemp,zTemp] = fun_line(xPosL,yPosL,zPosL,xPosR,yPosR,zPosR,t_);
                    
                    if (((b1<xTemp) and (xTemp<b2) and (b1<yTemp) and (yTemp < b2))):
                        DataFile.numFibersInMiddleSlice +=1
                        edgeCurveWasInMiddle=True
                        xMiddle = xTemp
                        yMiddle = yTemp
                        
                elif ((zPosL < midpoint) & (zPosR > midpoint)):    
                    t_ = (midpoint - zPosL) / (zPosR - zPosL);
                    [xTemp,yTemp,zTemp] = fun_line(xPosL,yPosL,zPosL,xPosR,yPosR,zPosR,t_);
                    
                    if (((b1<xTemp) and (xTemp<b2) and (b1<yTemp) and (yTemp < b2))):
                        DataFile.numFibersInMiddleSlice +=1
                        edgeCurveWasInMiddle=True
                        xMiddle = xTemp
                        yMiddle = yTemp
                
                
                #add bind sites for all left sides
                bindSitesPerEdge += DataFile.bindSitesPerNode[eR]
                
                #count alignment
                aveAlignment += DataFile.originalEdgeAlignment[edge]

                #add strain for subEdges
                undividedEdgeStrainTemp += DataFile.originalEdgeStrain[edge]
                
                #add force for subEdges
                undividedEdgeForceTemp += (DataFile.originalNodeForce[eR] + 
                                            DataFile.originalNodeForce[eL])
                

                distcurvetemp += np.linalg.norm(np.array([xPosL,yPosL,zPosL])-np.array([xPosR,yPosR,zPosR]))
                #average midoint
                aveMidPoint += DataFile.originalEdgeMidPoint[edge]
                aveMidPointScaled += edgeScaled
                #trial left and right
                if (eL < DataFile.originalNodeCount):
                    trueEdgeLeft = eL
                    bindSitesPerEdge += DataFile.bindSitesPerNode[eL]#last get bind sites for edge right

                if (eR < DataFile.originalNodeCount):
                    trueEdgeRight = eR
                    

                if ( (trueEdgeLeft < DataFile.originalNodeCount) & (trueEdgeRight < DataFile.originalNodeCount) ) :
                    
                    #                
                    zPosL = DataFile.zPos[trueEdgeLeft]
                    zPosR = DataFile.zPos[trueEdgeRight]
                    yPosL = DataFile.yPos[trueEdgeLeft]
                    yPosR = DataFile.yPos[trueEdgeRight]
                    xPosL = DataFile.xPos[trueEdgeLeft]
                    xPosR = DataFile.xPos[trueEdgeRight]
                    length_0 = dist=np.linalg.norm(np.array([xPosL,yPosL,zPosL])-np.array([xPosR,yPosR,zPosR]))
                    
                    DataFile.curvaturePerEdge.append(float(length_0/distcurvetemp));
                    
                    if (edgeCurveWasInMiddle==True):           
                        DataFile.curvaturePerEdgeMiddle.append(float(length_0/distcurvetemp))
                        DataFile.curvatureXPointPerEdgeMiddle.append(float(xMiddle/(b2-b1)))
                        DataFile.curvatureYPointPerEdgeMiddle.append(float(yMiddle/(b2-b1)))
                    
                    DataFile.noDivisionEdgesLeft.append(int(trueEdgeLeft))
                    DataFile.noDivisionEdgesRight.append(int(trueEdgeRight))
                    
                    
                    #average alignment
                    DataFile.noDivisionEdgesMidPoint.append(aveMidPoint/subEdgeCounter)
                    DataFile.noDivisionEdgesMidPointScaledForStrain.append(aveMidPointScaled/subEdgeCounter)

                    

                    #now that we have the correct points, get strain
                    DataFile.noDivisionEdgesStrain.append(undividedEdgeStrainTemp/(subEdgeCounter))
                    
                    #average force. Note: edges count 2 node forces, so double edge count
                    DataFile.noDivisionEdgesForce.append(undividedEdgeForceTemp/( 2 * subEdgeCounter))
                                                 
                    #count bind sites for original undivided edge
                    DataFile.noDivisionEdgesBindSites.append(bindSitesPerEdge)

                    #average alignment 
                    DataFile.noDivisionEdgesAlignment.append(aveAlignment / subEdgeCounter)
                    #then reset for next edge iteration
                    undividedEdgeCounter += 1
                    trueEdgeLeft = 2*(DataFile.originalNodeCount)
                    trueEdgeRight = 2*(DataFile.originalNodeCount)
                    undividedEdgeStrainTemp = 0.0
                    undividedEdgeForceTemp = 0.0
                    bindSitesPerEdge = 0
                    aveAlignment = 0
                    aveMidPoint = 0
                    aveMidPointScaled = 0
                    subEdgeCounter = 0
                    distcurvetemp=0.0
                    edgeCurveWasInMiddle=False
        #now everything is done, so we rescale the area
        area=(((b2-b1)*(b2-b1)))
        
        #area = ((DataFile.maxX-DataFile.minX) * (DataFile.maxY-DataFile.minY))
        DataFile.areaCoveredByFibers = 100*DataFile.numFibersInMiddleSlice * DataFile.FiberArea / area
        minNumOriginal = np.abs(np.min(DataFile.originalEdgeMidPointScaledForStrain))
        minNumNoDivision = np.abs(np.min(DataFile.noDivisionEdgesMidPointScaledForStrain))
        for i in range (len(DataFile.originalEdgeMidPointScaledForStrain)):    
            DataFile.originalEdgeMidPointScaledForStrain[i] += minNumOriginal
        for i in range (len(DataFile.noDivisionEdgesMidPointScaledForStrain)):    
            DataFile.noDivisionEdgesMidPointScaledForStrain[i] += minNumNoDivision
        
        print("original edge count")
        print(len(DataFile.originEdgeRight))
        #find unique binding sites
        #for i in range(len(DataFile.noDivisionEdgesLeft)):
        #think of a way to track unique edges that are bound to a given edge.     

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def data_generate(data_,rootdir):
    for subdir, dirs, files in os.walk(rootdir):
        colChoice = 0
        rowChoice = 0
        files.sort(key=natural_keys)

        for file in files:
            
            temp = DataFile(subdir,file)
            parserForceStrain(temp)
            data_.append(temp)
            print(file)

                    
#%%


#%%
LinkDir1 = "C:/Users/samuel_britton/Documents/collagen_elastin/eps=12/Params__max_nodes_29277_max_z_9_max_x_9_axis_0_axis_1_epsfactor_112"
NoLinkDir1 = "C:/Users/samuel_britton/Documents/collagen_elastin/eps=12/Params__max_nodes_29277_max_z_9_max_x_9_axis_1_axis_1_epsfactor_112"

#dataLink1 = pd.DataFrame((([numStrains*[DataFile] for i in range(numSamples)])),columns=colStrains)
#dataNoLink1 = pd.DataFrame((([numStrains*[DataFile] for i in range(numSamples)])),columns=colStrains)
data_axis0=[]
data_axis1=[]
data_generate(data_axis0,LinkDir1) 
data_generate(data_axis1,NoLinkDir1) 
#%%
mpa_axis0=[]
len_axis0=[]
force_axis0=[]
for i in range(0,len(data_axis0)):
    len_axis0.append(data_axis0[i].maxZ - data_axis0[i].minZ)
    mpa_axis0.append(data_axis0[i].MPascals_xz)
    force_axis0.append(data_axis0[i].force_upper  +data_axis0[i].force_lower)

mpa_axis1=[]
len_axis1=[]
force_axis1=[]
for i in range(0,len(data_axis1)):
    len_axis1.append(data_axis1[i].maxX - data_axis1[i].minX)
    mpa_axis1.append(data_axis1[i].MPascals_xy)
    force_axis1.append(data_axis1[i].force_upper + data_axis1[i].force_lower)
    
#%%
    
fig = plt.figure(figsize=(figure_norm, figure_norm))

ax1 = fig.add_subplot(111)

ax1.scatter(len_axis0/len_axis0[0] - 1.0, mpa_axis0, label='Axis Aligned')
ax1.scatter(len_axis1/len_axis1[0] - 1.0, mpa_axis1, label='Axis Not Aligned')

ax1.set_ylabel('Stress $\sigma$, MPa',fontsize = Fontsize_Sub)
#ax1.set_ylabel('Force nN',fontsize = Fontsize_Sub)
ax1.set_xlabel('Network Strain, $\Gamma$',fontsize=Fontsize_Sub)
ax1.legend(fontsize=Fontsize_Leg)


#%%
axis1_collagen = "C:/Users/samuel_britton/Documents/collagen_elastin/collagen/Params__max_nodes_29277_max_z_9_max_x_9_axis_1_pullwidth_1_epsfactor_1130_collagen_30"
axis0_collagen = "C:/Users/samuel_britton/Documents/collagen_elastin/collagen/Params__max_nodes_29277_max_z_9_max_x_9_axis_0_pullwidth_1_epsfactor_1130_collagen_30"

#dataNoLink1 = pd.DataFrame((([numStrains*[DataFile] for i in range(numSamples)])),columns=colStrains)
data_collagen_axis0=[]
data_collagen_axis1=[]
data_generate(data_collagen_axis0,axis0_collagen) 
data_generate(data_collagen_axis1,axis1_collagen) 
#%%
axis1_collagen50 = "C:/Users/samuel_britton/Documents/collagen_elastin/collagen/Params__max_nodes_29277_max_z_9_max_x_9_axis_1_pullwidth_1_epsfactor_1150_collagen_30"
axis0_collagen50 = "C:/Users/samuel_britton/Documents/collagen_elastin/collagen/Params__max_nodes_29277_max_z_9_max_x_9_axis_0_pullwidth_1_epsfactor_1150_collagen_30"

#dataNoLink1 = pd.DataFrame((([numStrains*[DataFile] for i in range(numSamples)])),columns=colStrains)
data_elastin_axis0=[]
data_elastin_axis1=[]
data_generate(data_elastin_axis0,axis0_collagen50) 
data_generate(data_elastin_axis1,axis1_collagen50) 


#%%

mpa_collagen_axis0=[]
len_collagen_axis0=[]
force_collagen_axis0=[]
for i in range(0,len(data_collagen_axis0)):
    len_collagen_axis0.append(data_collagen_axis0[i].maxZ - data_collagen_axis0[i].minZ)
    mpa_collagen_axis0.append(data_collagen_axis0[i].MPascals_xz)
    force_collagen_axis0.append(data_collagen_axis0[i].force_upper  +data_collagen_axis0[i].force_lower)

mpa_collagen_axis1=[]
len_collagen_axis1=[]
force_collagen_axis1=[]
for i in range(0,len(data_collagen_axis1)):
    len_collagen_axis1.append(data_collagen_axis1[i].maxX - data_collagen_axis1[i].minX)
    mpa_collagen_axis1.append(data_collagen_axis1[i].MPascals_xy)
    force_collagen_axis1.append(data_collagen_axis1[i].force_upper + data_collagen_axis1[i].force_lower)

mpa_elastin_axis0=[]
len_elastin_axis0=[]
force_elastin_axis0=[]
for i in range(0,len(data_elastin_axis0)):
    len_elastin_axis0.append(data_elastin_axis0[i].maxZ - data_elastin_axis0[i].minZ)
    mpa_elastin_axis0.append(data_elastin_axis0[i].MPascals_xz)
    force_elastin_axis0.append(data_elastin_axis0[i].force_upper  +data_elastin_axis0[i].force_lower)

mpa_elastin_axis1=[]
len_elastin_axis1=[]
force_elastin_axis1=[]
for i in range(0,len(data_elastin_axis1)):
    len_elastin_axis1.append(data_elastin_axis1[i].maxX - data_elastin_axis1[i].minX)
    mpa_elastin_axis1.append(data_elastin_axis1[i].MPascals_xy)
    force_elastin_axis1.append(data_elastin_axis1[i].force_upper + data_elastin_axis1[i].force_lower)
#%%
fig = plt.figure(figsize=(figure_norm, figure_norm))

ax1 = fig.add_subplot(111)

#ax1.scatter(len_axis0/len_axis0[0] - 1.0, mpa_axis0, label='Axis Aligned, mixed')
#ax1.scatter(len_axis1/len_axis1[0] - 1.0, mpa_axis1, label='Axis Not Aligned, mixed')

ax1.scatter(len_collagen_axis0[0:-30]/len_collagen_axis0[0] - 1.0, mpa_collagen_axis0[0:-30], label='Axis Aligned, eps=1.13')
ax1.scatter(len_collagen_axis1/len_collagen_axis1[0] - 1.0, mpa_collagen_axis1, label='Axis Not Aligned, eps=1.13')

ax1.scatter(len_elastin_axis0/len_elastin_axis0[0] - 1.0, mpa_elastin_axis0, label='Axis Aligned, eps=1.15')
ax1.scatter(len_elastin_axis1/len_elastin_axis1[0] - 1.0, mpa_elastin_axis1, label='Axis Not Aligned, eps=1.15')


ax1.set_ylabel('Stress $\sigma$, MPa',fontsize = Fontsize_Sub)
#ax1.set_ylabel('Force nN',fontsize = Fontsize_Sub)
ax1.set_xlabel('Network Strain, $\Gamma$',fontsize=Fontsize_Sub)
ax1.legend(fontsize=Fontsize_Leg)
