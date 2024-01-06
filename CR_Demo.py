#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 23:03:44 2022

@author: sawa
"""


import sys
import os

from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QComboBox, QPushButton, QLineEdit, QPlainTextEdit, QWidget
from PyQt5.QtGui import QIcon, QPixmap, QFont
from PyQt5.QtCore import pyqtSlot


import pandas as pd

#%% Import data

# Get working directory for full datapath:
workdir = os.getcwd()

# Import material and gasses/medium data:
materials_data = pd.read_csv(workdir + '/MaterialProperties.csv')
gasses_data = pd.read_csv(workdir + '/GassesProperties.csv')


#%% Functions

# Gather relevant parameters into objects:
def get_matconst(mat1_name,
                 mat2_name,
                 rgh1,
                 rgh2,
                 mdm,
                 mat_list,
                 medium_list,
                 materials_data,
                 gasses_data,
                 force,
                 area):
    
    ### Material 1:
    mat1_index = mat_list.index(mat1_name)
    
    class Mat1:
        E     = float(materials_data["Young's modulus"][mat1_index + 2])
        nu    = float(materials_data["Poisson's ratio"][mat1_index + 2])
        Hk    = float(materials_data["Knoop hardness"][mat1_index + 2])
        k     = float(materials_data["Conductivity"][mat1_index + 2])
        eps   = materials_data["emissivity"][mat1_index + 2]
        Ra    = float(rgh1)
        sigma = 1.25*float(rgh1)
        
        if float(rgh1) > 0.0000016/1.25:
            m = (0.171/1.25)*(1.25*float(rgh1)*1000000)**0.402 
        else:
            m = (0.183/1.25)*(1.25*float(rgh1)*1000000)**0.743
            
    
    ### Material 2:
    mat2_index = mat_list.index(mat2_name)
    
    class Mat2:
        E     = float(materials_data["Young's modulus"][mat2_index + 2])
        nu    = float(materials_data["Poisson's ratio"][mat2_index + 2])
        Hk    = float(materials_data["Knoop hardness"][mat2_index + 2])
        k     = float(materials_data["Conductivity"][mat2_index + 2])
        eps   = materials_data["emissivity"][mat2_index + 2]
        Ra    = float(rgh2)
        sigma = 1.25*float(rgh2)
        
        if float(rgh2) > 0.0000016/1.25:
            m = (0.171/1.25)*(1.25*float(rgh2)*1000000)**0.402 
        else:
            m = (0.183/1.25)*(1.25*float(rgh2)*1000000)**0.743
            
    
    ### Harmonic averages:
    class har_avg:
        E     = 1 / ( (1 - Mat1.nu**2)/Mat1.E + (1 - Mat2.nu**2)/Mat2.E )
        nu    = ( 0.5*Mat1.nu**2 + 0.5*Mat2.nu**2 )**0.5
        Hk    = min(Mat1.Hk,Mat2.Hk)
        ks    = 2*Mat1.k*Mat2.k/(Mat1.k + Mat2.k)
        sigma = ( Mat1.sigma**2 + Mat2.sigma**2 )**0.5
        m     = ( Mat1.m**2 + Mat2.m**2 )**0.5
        
    ### Medium in gaps:
    mdm_index = medium_list.index(mdm)
        
    class cont_par:
        medium = mdm
        k      = float(gasses_data["Conductivity"][mdm_index+2])
        F      = force
        A      = area
        P      = force / area
        delta  = 1.91 * har_avg.sigma * ( 10000000 / har_avg.Hk )**( -0.097 )
        
    return(Mat1, Mat2, har_avg, cont_par)


# Mikic plastic thermal contact resistance model:
def mikic_plastic(har_avg, cont_par):
    hc = 1.25 * har_avg.ks * ( har_avg.m / har_avg.sigma ) * ( cont_par.P / har_avg.Hk )**0.95
    return(hc)

# Mikic elastic thermal contact resistance model:
def mikic_elastic(har_avg, cont_par):
    hc = 1.55 * har_avg.ks * ( har_avg.m / har_avg.sigma ) * ( (cont_par.P * 2**0.5) / (har_avg.E*har_avg.m) )**0.94
    hc1 = 2.83 * har_avg.ks * har_avg.sigma**-0.955 * ( cont_par.P / har_avg.E )**0.94
    hc2 = 1.9 * ( har_avg.ks / har_avg.sigma ) * ( cont_par.P / har_avg.E )**0.94
    return(hc, hc1, hc2)

# Yavonovich thermal contact resistance model:
def yavonovich(har_avg, cont_par):
    hc = 1.25 * har_avg.ks * ( har_avg.m / har_avg.sigma ) * ( cont_par.P / har_avg.Hk )**0.95
    return(hc)
    
    
    

#%% Create GUI

# Break-out window: This class is not relevant at the moment!
class References(QWidget):    
    
    def __init__(self):
        super().__init__()
        
        # Window geometry:
        self.setGeometry(800,200,600,600)
        self.setWindowTitle("Contact Resistance References - Demo")

# Main window:
class Window(QMainWindow):
    
    def __init__(self):
        super().__init__()
        
        # Lists of relevant stuff:
        mat_list       = list(materials_data['Solids'][2::])
        medium_list    = list(gasses_data['Gasses'][2::])
        
  
        #######################################################################
        ##### Material 1:                                                 #####
        #######################################################################
        
        # Drop-down menu material 1:
        self.ddm_mat1 = QComboBox(self)
        self.ddm_mat1.addItems(mat_list)
        self.ddm_mat1.setEditable(False)
        self.ddm_mat1.setInsertPolicy(QComboBox.NoInsert)
        self.ddm_mat1.move(30, 40)

        # Label menu material 1:
        lbl_mat1 = QLabel(self)
        lbl_mat1.setText("Select material 1:")
        lbl_mat1.setGeometry(30,10,150,40)


        # Set roughness material 1:
        self.rgh1 = QLineEdit(self)
        self.rgh1.setGeometry(30,100,100,30)

        # Label menu roughness material 1:
        lbl_rgh1 = QLabel(self)
        lbl_rgh1.setText("Roughness material 1:")
        lbl_rgh1.setGeometry(30,70,150,40)
        
        # Label unit roughness material 1:
        unit_rgh1 = QLabel(self)
        unit_rgh1.setText("Ra (m)")
        unit_rgh1.setGeometry(140,100,40,30)
        
        
        #######################################################################
        ##### Material 2:                                                 #####
        #######################################################################        
        
        # Drop-down menu material 2:
        self.ddm_mat2 = QComboBox(self)
        self.ddm_mat2.addItems(mat_list)
        self.ddm_mat2.setEditable(False)
        self.ddm_mat2.setInsertPolicy(QComboBox.NoInsert)
        self.ddm_mat2.move(230, 40)

        # Label menu material 1:
        lbl_mat2 = QLabel(self)
        lbl_mat2.setText("Select material 2:")
        lbl_mat2.setGeometry(230,10,150,40)
        
        # Set roughness material 2:
        # Input box for contact pressure:
        self.rgh2 = QLineEdit(self)
        self.rgh2.setGeometry(230,100,100,30)

        # Label menu roughness material 2:
        lbl_rgh2 = QLabel(self)
        lbl_rgh2.setText("Roughness material 2:")
        lbl_rgh2.setGeometry(230,70,150,40)
        
        # Label unit roughness material 2:
        unit_rgh2 = QLabel(self)
        unit_rgh2.setText("Ra (m)")
        unit_rgh2.setGeometry(340,100,40,30)


        #######################################################################
        ##### Contact interface parameters:                               #####
        #######################################################################

        # Input box for contact area:
        self.inp_area = QLineEdit(self)
        self.inp_area.setGeometry(30,170,190,30)

        # Label for input contact area:
        lbl_area = QLabel(self)
        lbl_area.setText("Contact area:")
        lbl_area.setGeometry(30,140,150,40)
        
        # Unit contact area:
        unit_area = QLabel(self)
        unit_area.setText("m^2")
        unit_area.setGeometry(230,175,40,20)

    
        # Input box for contact force:
        self.inp_force = QLineEdit(self)
        self.inp_force.setGeometry(30,225,190,30)

        # Label for input contact force:
        lbl_force = QLabel(self)
        lbl_force.setText("Contact force:")
        lbl_force.setGeometry(30,195,150,40)
        
        # Unit contact force:
        unit_force = QLabel(self)
        unit_force.setText("N")
        unit_force.setGeometry(230,230,40,20)
    

        # Drop-down menu medium in gaps:
        self.ddm_mdm = QComboBox(self)
        self.ddm_mdm.addItems(medium_list)
        self.ddm_mdm.setEditable(False)
        self.ddm_mdm.setInsertPolicy(QComboBox.NoInsert)
        self.ddm_mdm.move(430, 40)
        
        # Label medium in gaps:
        lbl_mdm = QLabel(self)
        lbl_mdm.setText("Medium between gaps:")
        lbl_mdm.setGeometry(430,10,150,40)
        
        
        #######################################################################
        ##### Reduction factor:                                           #####
        #######################################################################
        
        # Input box for reduction factor:
        self.reduction_factor = QLineEdit('10',self)
        self.reduction_factor.setGeometry(340,170,50,30)

        # Label for input contact force:
        lbl_redfac = QLabel(self)
        lbl_redfac.setText("Reduction factor:")
        lbl_redfac.setGeometry(340,140,150,40)
        
        # Unit contact force:
        redfac_text = QLabel(self)
        redfac_text.setText("(typically between 7-30)")
        redfac_text.setGeometry(410,175,200,20)
        
        
        warning_lbl = QLabel(self)
        warning_lbl.setText("Please use numerical values")
        warning_lbl.setGeometry(340,230,250,20)
        warning_lbl.setFont(QFont('Arial',13))

        
        #######################################################################
        ##### Image:                                                      #####
        #######################################################################
       
        # Add figure:
        im_holder = QLabel(self)
        im = QPixmap("./ContactResistance.PNG")
        im_holder.setPixmap(im)
        im_holder.setGeometry(30,280,380,150)
    
        
        #######################################################################
        ##### Display results:                                            #####
        #######################################################################
        
        # Label to show results:
        ### Row 1:
        self.CR = QLabel(self)
        self.CR.setGeometry(30,450,300,20)
        self.CR.setText("Contact Resistance:")
        
        self.lbl_theory = QLabel(self)
        self.lbl_theory.setGeometry(240,450,300,20)
        self.lbl_theory.setText("Theory")
        
        self.lbl_theory = QLabel(self)
        self.lbl_theory.setGeometry(330,450,300,20)
        self.lbl_theory.setText("Conservative")
        
        self.lbl_theory = QLabel(self)
        self.lbl_theory.setGeometry(420,450,300,20)
        self.lbl_theory.setText("Empirical")
        
        ### Row 2:
        self.result_me_original = QLabel(self)
        self.result_me_original.setGeometry(30,480,300,20)
        self.result_me_original.setText("Mikic elastic (Original):")
        
        self.res_meori_theory = QLabel(self)
        self.res_meori_theory.setGeometry(240,480,300,20)
        self.res_meori_theory.setText("")
        
        self.res_meori_conservative = QLabel(self)
        self.res_meori_conservative.setGeometry(330,480,300,20)
        self.res_meori_conservative.setText("")
        
        self.res_meori_empirical = QLabel(self)
        self.res_meori_empirical.setGeometry(420,480,300,20)
        self.res_meori_empirical.setText("")
        
        self.res_meori_lbl = QLabel(self)
        self.res_meori_lbl.setGeometry(510,480,300,20)
        self.res_meori_lbl.setText("W/m^2/K")
        
        ### Row 3:
        self.result_me_lowRq = QLabel(self)
        self.result_me_lowRq.setGeometry(30,510,300,20)
        self.result_me_lowRq.setText("Mikic elastic (Low roughness):")
    
        self.res_melrq_theory = QLabel(self)
        self.res_melrq_theory.setGeometry(240,510,300,20)
        self.res_melrq_theory.setText("")
        
        self.res_melrq_conservative = QLabel(self)
        self.res_melrq_conservative.setGeometry(330,510,300,20)
        self.res_melrq_conservative.setText("")
        
        self.res_melrq_empirical = QLabel(self)
        self.res_melrq_empirical.setGeometry(420,510,300,20)
        self.res_melrq_empirical.setText("")
        
        self.res_melrq_lbl = QLabel(self)
        self.res_melrq_lbl.setGeometry(510,510,300,20)
        self.res_melrq_lbl.setText("W/m^2/K")

        ### Row 4:
        self.result_me_simplified = QLabel(self)
        self.result_me_simplified.setGeometry(30,540,300,20)
        self.result_me_simplified.setText("Mikic elastic (Simplified):")
        
        self.res_mesim_theory = QLabel(self)
        self.res_mesim_theory.setGeometry(240,540,300,20)
        self.res_mesim_theory.setText("")
        
        self.res_mesim_conservative = QLabel(self)
        self.res_mesim_conservative.setGeometry(330,540,300,20)
        self.res_mesim_conservative.setText("")
        
        self.res_mesim_empirical = QLabel(self)
        self.res_mesim_empirical.setGeometry(420,540,300,20)
        self.res_mesim_empirical.setText("")
        
        self.res_mesim_lbl = QLabel(self)
        self.res_mesim_lbl.setGeometry(510,540,300,20)
        self.res_mesim_lbl.setText("W/m^2/K")
        
        ### Row 5:
        self.result_mikicplastic = QLabel(self)
        self.result_mikicplastic.setGeometry(30,570,300,20)
        self.result_mikicplastic.setText("Mikic plasitc:")
        
        self.res_mp_theory = QLabel(self)
        self.res_mp_theory.setGeometry(240,570,300,20)
        self.res_mp_theory.setText("")
        
        self.res_mp_conservative = QLabel(self)
        self.res_mp_conservative.setGeometry(330,570,300,20)
        self.res_mp_conservative.setText("")
        
        self.res_mp_empirical = QLabel(self)
        self.res_mp_empirical.setGeometry(420,570,300,20)
        self.res_mp_empirical.setText("")
        
        self.res_mp_lbl = QLabel(self)
        self.res_mp_lbl.setGeometry(510,570,300,20)
        self.res_mp_lbl.setText("W/m^2/K")
        
        ### Row 6:
        self.result_yovanovich = QLabel(self)
        self.result_yovanovich.setGeometry(30,600,300,20)
        self.result_yovanovich.setText("Yovanovich:")
        
        self.res_yov_theory = QLabel(self)
        self.res_yov_theory.setGeometry(240,600,300,20)
        self.res_yov_theory.setText("")
        
        self.res_yov_conservative = QLabel(self)
        self.res_yov_conservative.setGeometry(330,600,300,20)
        self.res_yov_conservative.setText("")
        
        self.res_yov_empirical = QLabel(self)
        self.res_yov_empirical.setGeometry(420,600,300,20)
        self.res_yov_empirical.setText("")
        
        self.res_yov_lbl = QLabel(self)
        self.res_yov_lbl.setGeometry(510,600,300,20)
        self.res_yov_lbl.setText("W/m^2/K")
        
        ### Row 7:
        self.result_final = QLabel(self)
        self.result_final.setGeometry(30,640,300,20)
        self.result_final.setText("Values to be used:")
        
        self.res_fin_theory = QLabel(self)
        self.res_fin_theory.setGeometry(240,640,300,20)
        self.res_fin_theory.setText("")
        
        self.res_fin_conservative = QLabel(self)
        self.res_fin_conservative.setGeometry(330,640,300,20)
        self.res_fin_conservative.setText("")
        
        self.res_fin_empirical = QLabel(self)
        self.res_fin_empirical.setGeometry(420,640,300,20)
        self.res_fin_empirical.setText("")
        
        self.res_fin_lbl = QLabel(self)
        self.res_fin_lbl.setGeometry(510,640,300,20)
        
        # Create button to compute contact resistance:
        comp_cr = QPushButton(self)
        comp_cr.setText("Compute Contact Resistance")
        comp_cr.setGeometry(30,680,200,30)
        comp_cr.clicked.connect(self.button_pressed)
        
        
        #######################################################################
        ##### References contact resistance:                              #####
        #######################################################################
        
        # This doesn't work yet!
        # self.ref_win = References()
        
        # # Button to get references:
        # ref_but = QPushButton(self)
        # ref_but.setText("References")
        # ref_but.setGeometry(280,530,80,30)
        # ref_but.clicked.connect(self.ref_win.show())
        
        
        #######################################################################
        ##### Message box:                                                #####
        #######################################################################
        
        # # Label message box:
        # msg_lbl = QLabel(self)
        # msg_lbl.setText('Messages:')
        # msg_lbl.setGeometry(30,620,300,20)
        
        # # Message box: 
        # self.msg_box = QPlainTextEdit(self)
        # self.msg_box.setGeometry(30,650,350,120)
        # self.msg_box.setReadOnly(True)
        
    
        #######################################################################
        ##### Window size:                                                #####
        #######################################################################

        # Window geometry:
        self.setGeometry(200,200,600,740)
        self.setWindowTitle("Contact Resistance Calculator - Demo")
        self.show()
   
     
    def button_pressed(self):
        content = [self.ddm_mat1.currentText(),
                   self.rgh1.text(),
                   self.ddm_mat2.currentText(),
                   self.rgh2.text(),
                   self.ddm_mdm.currentText(),
                   self.inp_area.text(),
                   self.inp_force.text(),
                   self.reduction_factor.text()]
        
        print(content)
        
        mat1 = self.ddm_mat1.currentText()
        mat2 = self.ddm_mat2.currentText()
        rgh1 = self.rgh1.text()
        rgh2 = self.rgh2.text()
        mdm  = self.ddm_mdm.currentText()
        red_fac = self.reduction_factor.text()
        
        print(type(rgh1))
        
        # Get objects with properties:
        [Mat1, Mat2, har_avg, cont_par] = get_matconst(mat1,
                                                       mat2,
                                                       rgh1,
                                                       rgh2,
                                                       mdm,
                                                       list(materials_data['Solids'][2::]),
                                                       list(gasses_data['Gasses'][2::]),
                                                       materials_data,
                                                       gasses_data,
                                                       float(self.inp_force.text()),
                                                       float(self.inp_area.text()))
        
        # Compute contact resistance:
        [meori_theory, melrq_theory, mesim_theory] = mikic_elastic(har_avg, cont_par)
        mp_theory = mikic_plastic(har_avg, cont_par)
        yov_theory = yavonovich(har_avg, cont_par)
        
        fin_theory = (meori_theory + melrq_theory + mesim_theory) / 3
        
        meori_conservative = meori_theory / float(red_fac)
        melrq_conservative = melrq_theory / float(red_fac)
        mesim_conservative = mesim_theory / float(red_fac)
        mp_conservative    = mp_theory / float(red_fac)
        yov_conservative   = yov_theory / float(red_fac)
        
        fin_conservative = (meori_conservative + melrq_conservative + mesim_conservative) / 3
        
        self.res_meori_theory.setText(f'{meori_theory:.7}')
        self.res_melrq_theory.setText(f'{melrq_theory:.7}')
        self.res_mesim_theory.setText(f'{mesim_theory:.7}')
        self.res_mp_theory.setText(f'{mp_theory:.7}')
        self.res_yov_theory.setText(f'{yov_theory:.7}')
        self.res_fin_theory.setText(f'{fin_theory:.7}')
        
        
        self.res_meori_conservative.setText(f'{meori_conservative:.7}')
        self.res_melrq_conservative.setText(f'{melrq_conservative:.7}')
        self.res_mesim_conservative.setText(f'{mesim_conservative:.7}')
        self.res_mp_conservative.setText(f'{mp_conservative:.7}')
        self.res_yov_conservative.setText(f'{yov_conservative:.7}')
        self.res_fin_conservative.setText(f'{fin_conservative:.7}')
        
        
        self.res_meori_empirical.setText("---")
        self.res_melrq_empirical.setText("---")
        self.res_mesim_empirical.setText("---")
        self.res_mp_empirical.setText("---")
        self.res_yov_empirical.setText("---")
        self.res_fin_empirical.setText("---")
        
     
    # Not used at the moment!        
    def get_references(self):
        self.ref_win.show()
   

if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = Window()
    sys.exit(app.exec_())
