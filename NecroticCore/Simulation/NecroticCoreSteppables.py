
from PySteppables import *
import CompuCell
import sys
from copy import deepcopy
import random
from math import *
from PlayerPython import *
import os
import numpy as np


import time
import os.path



V0=16.0
S0=16.0
LBD_V0=15.0
LBD_S0=5.0
V0_nec = 4.0 
S0_nec = 4.0 
LBD_NV0 = 10.0  
LBD_NS0 = 5.0 
incvol=0.2
decvol=0.01
volmaxmit=32
Svolmaxmit=32
ktgs=4.0
maxdiv=8
probstem=0.2
probmut=0.1
PGrThr0=0.032
SGrThr0=0.032
PNeThr0=102.0
QNeThr0=204.0
SNeThr0=408.0
QSNeThr0=916.0
QPThr0=79.5
QSSThr0=79.5
GluD=0.0032
cadhstdev=2.0

def getCellCOMPoint3D(cell):
    pt=CompuCell.Point3D()
    pt.x=int(round(cell.xCOM))
    pt.y=int(round(cell.yCOM))
    pt.z=int(round(cell.zCOM))
    
    return pt

class TGrowthSteppable(SteppableBasePy):
    
    
    
    def start(self):    
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM):
            cell.targetVolume=V0
            cell.lambdaVolume=LBD_V0
            cell.targetSurface=S0
            cell.lambdaSurface=LBD_S0
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Counter"] = 0
            cellDict["Health"] = 0
            
        for cell in self.cellListByType(self.NECROTIC):
            cell.targetVolume=V0_nec
            cell.lambdaVolume=LBD_NV0
            cell.targetSurface=S0_nec
            cell.lambdaSurface=LBD_NS0
            
    def step(self,mcs):    
        glucoseField=CompuCell.getConcentrationField(self.simulator,"Glucose")
        
        for cell in self.cellList: 
          pt = getCellCOMPoint3D(cell)
          glucoseField=CompuCell.getConcentrationField(self.simulator,"Glucose")
          conc = glucoseField.get(pt)
          cellDict = CompuCell.getPyAttrib(cell)
          if cell.type == self.NECROTIC:
              cell.targetVolume-=min(decvol, cell.targetVolume)
              cell.targetSurface = ktgs*sqrt(cell.targetVolume)
          if cell.type == self.PCANCER:
              cell.targetVolume+=incvol*max(0, conc-PGrThr0)
              cell.targetSurface =ktgs*sqrt(cell.targetVolume)
          if cell.type == self.PSTEM:
              cell.targetVolume+=incvol*max(0, conc-SGrThr0)
              cell.targetSurface =ktgs*sqrt(cell.targetVolume)   
    

    
class StravHealthCalculator(SteppableBasePy):
     
    def start(self):
            
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM, self.NECROTIC):
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Starv"]=0
            cellDict["Health"]=0
            
    def MM(self,x,m,k):
        return (m*x/(x+k))
        
        for cell in self.cellList:
            if cell.type!=self.NECROTIC:
                cellDict = CompuCell.getPyAttrib(cell)
                pt = getCellCOMPoint3D(cell)
                conc = glucoseField.get (pt)
            if cell.type == self.PCANCER: 
               if conc < GluD:
                   cellDict["Starv"]+=abs(self.MM(conc,PUgMax,GluK)\
                   -self.MM(GluD,PUgMax,GluK))
               else:
                   cellDict["Health"]+=self.MM(conc,PUgMax,GluK)\
                   -self.MM(GluD,PUgMax,GluK)
            if cell.type == self.QCANCER: 
               if conc < GluD:
                   cellDict["Starv"]+=abs(self.MM(conc,QUgMax,GluK)\
                   -self.MM(GluD,QUgMax,GluK))
               else:
                   cellDict["Health"]+=self.MM(conc,QUgMax,GluK)\
                   -self.MM(GluD,QUgMax,GluK)
            if cell.type == self.PSTEM:
               if conc < GluD:
                   cellDict["Starv"]+=abs(self.MM(conc,SUgMax,GluK)\
                   -self.MM(GluD, SUgMax, GluK))
               else:
                   cellDict["Health"]+=self.MM(conc, SUgMax, GluK)\
                   -self.MM(GluD, SUgMax, GluK)
            if cell.type == self.QSTEM:
               if conc < GluD:
                   cellDict["Starv"]+=abs(self.MM(conc, QSUgMax, GluK)\
                   -self.MM(GluD, QSUgMax, GluK))
               else:
                   cellDict["Health"]+=self.MM(conc, QSUgMax, GluK)\
                   -self.MM(GluD, QSUgMax, GluK)
        
                        
class CellStateTransition(SteppableBasePy):
     
    def step(self,mcs): 
        for cell in self.cellList: 
          cellDict=CompuCell.getPyAttrib(cell)
          if cell.type == self.PCANCER: 
            if cellDict["Starv"] > PNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
          if cell.type == self. QCANCER: 
            if cellDict["Starv"] > QNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
            if cellDict["Health"] > QPThr0:
              cell.type = self.PCANCER
              cellDict["Health"]=0 
          if cell.type == self. PSTEM: 
            if cellDict["Starv"] > SNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
          if cell.type == self. QSTEM: 
            if cellDict["Starv"] > QSNeThr0:
              cell.type = self.NECROTIC
              cellDict["Health"]=0 
            if cellDict["Health"] > QSSThr0:
              cell.type = self.PSTEM
              cellDict["Health"]=0   
                     

from PySteppables import *
from PySteppablesExamples import MitosisSteppableBase
import CompuCell
import sys

from PlayerPython import *
from math import *

class MitosisSteppable(MitosisSteppableBase):
    def  start(self):
             
        for cell in self.cellListByType(self.PCANCER, self.QCANCER, self.PSTEM, self.QSTEM):
         
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["Counter"]=0    
        
    
    def step(self,mcs):
        
        cells_to_divide=[]
        
        
        for cell in self.cellList:
            if ((cell.type==self.PCANCER or cell.type==self.QCANCER)\
            and cell.volume>volmaxmit)\
            or ((cell.type==self.PSTEM or cell.type==self.QSTEM)\
            and cell.volume>Svolmaxmit):
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
                self.divideCellRandomOrientation(cell)
          
    def updateAttributes(self):
        
        parentCell = self.mitosisSteppable.parentCell
        childCell = self.mitosisSteppable.childCell
        parentCell.targetVolume=V0
        childCell.targetVolume= V0
        parentCell.targetSurface=ktgs*sqrt(parentCell.targetVolume)
        childCell.targetSurface= ktgs*sqrt(childCell.targetVolume)           
        parentCell.lambdaVolume=LBD_V0
        parentCell.lambdaSurface=LBD_S0
        childCell.lambdaVolume=LBD_V0
        childCell.lambdaSurface=LBD_S0
        parentCellDict=CompuCell.getPyAttrib(parentCell)
        childCellDict=CompuCell.getPyAttrib(childCell)
        
        if parentCell.type==self.PCANCER:
            
          temp = random.gauss(maxdiv,2)
          if parentCellDict["Counter"]<= temp:
             parentCell.type=self.QCANCER #both cells are QC after mitosis 
             childCell.type=self.QCANCER
             parentCellDict["Counter"]+=1
            
          childCellDict["Counter"] =deepcopy(parentCellDict["Counter"])
            
          if parentCellDict["Counter"]> temp:
              parentCell.type=self.NECROTIC
              childCell.type=self.NECROTIC
              
        if parentCell.type == self.PSTEM or parentCell.type == self.QSTEM: 
            parentCellDict["Counter"]+=1
            parentCell.type=self.QSTEM #one is QS
            childCell.type=self.QCANCER #the other is QC
            
            if random.random()<=probstem: #0.2 chance for a 2nd QS
                childCell.type=self.QSTEM
            childCellDict["Counter"]=0
            
        parentCellDict["Starv"]=0
        childCellDict["Starv"]=0
        parentCellDict["Health"]=0
        childCellDict["Health"]=0 
            
        jcadh = self.adhesionFlexPlugin.getAdhesionMoleculeDensityByIndex(parentCell,0)
        jint = self.adhesionFlexPlugin.getAdhesionMoleculeDensityByIndex(parentCell,1)
        jFN = self.adhesionFlexPlugin.getAdhesionMoleculeDensityByIndex(parentCell,2) 
        if parentCell.type!=self.NECROTIC:
          r = random.random()
          if r<probmut:
              new_cadh = random.gauss(jcadh,cadhstdev)
              if new_cadh>=0 and new_cadh<=16:
                jcadh = new_cadh
          r = random.random()
          if r<probmut:
              new_int = random.gauss(jint,cadhstdev)
              if new_int>=0 and new_int<=16:
                jint = new_int
        #Setting new adhesion density for product cell
          self.adhesionFlexPlugin.assignNewAdhesionMoleculeDensityVector (parentCell,[jcadh,jint,jFN])
          self.adhesionFlexPlugin.assignNewAdhesionMoleculeDensityVector (childCell,[jcadh,jint,jFN])
        
