# -*- coding: utf-8 -*- 

from PyQt5.QtWidgets import (QWidget, QApplication, QGroupBox,QTabWidget, QPushButton, 
QLabel, QHBoxLayout,  QVBoxLayout, QGridLayout, QFormLayout, QLineEdit, QTextEdit,QComboBox,QMessageBox,
    QDesktopWidget, QFileDialog,QTextEdit, QMessageBox, QProgressBar,QShortcut,
    QCheckBox)
from PyQt5.QtGui import QPixmap, QIcon, QFont, QGuiApplication, QTextBlockFormat, QTextCursor
from PyQt5.QtCore import QThread,pyqtSignal,Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar  
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sys
import os
from numpy import *
                        #used to store and reload fig
from TRANK import rms_error_spectrum, nk_plot, try_mkdir, functionize_nk_file, extrap, extrap_c, error_adaptive_iterative_fit_spectra
import time
from PyQt5.Qt import QCoreApplication, QKeyEvent, QEvent
import io
import re
import webbrowser

from multiprocessing import Manager, cpu_count
from multiprocessing.managers import BaseManager
from new_basic_setup import dyn_basic_setup



class MyManager(BaseManager):
            pass
MyManager.register('dyn_basic_setup',dyn_basic_setup)



class MyCustomToolbar(NavigationToolbar):
            toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Back','Forward','Pan', 'Zoom','Save')]

    #QThread for TRANK_fit
class TRANK_fit_Thread(QThread):
    
    Output_axes=pyqtSignal(tuple)
    Thickness_fit=pyqtSignal(tuple)
    Invalid_Input=pyqtSignal()
    def __init__(self,fig,fig2,fig3,TB1,t_list,R_line,R_dir,T_dir,nk_dir,if_scan,if_metal,if_e, set_list,lamda_list):
        QThread.__init__(self)
       
        self.fig=fig
        self.fig2=fig2
        self.fig3=fig3
        self.TB1=TB1
        self.nk_dir=nk_dir
        self.if_scan=if_scan
        self.if_metal=if_metal
        self.if_e=if_e
        self.t_list=t_list
        self.set_list=set_list
        self.lamda_list=lamda_list
        #from new_basic_setup  import dyn_basic_setup
        self.input_data=[int(t_list[0]),R_line,R_dir,T_dir]
        self.fig.clear()
        self.fig.canvas.draw()
        self.fig2.clear()
        self.fig2.canvas.draw()
        self.fig3.clear()
        self.fig3.canvas.draw()
        self.TB1.clear()
        #QCoreApplication.processEvents()
    
    ###### okay so this goofy stuff should allow us to update the thickness of our layer from other places
    
    def __del__(self):
        
        self.wait(3)
        
    def run(self):
        if self.if_scan :
          
            self.scan_thickness(int(self.t_list[0]),int(self.t_list[1]),int(self.t_list[2]))
        else:
            self.thickness_fit(self.if_metal,self.set_list)
    
    
    def scan_thickness(self,max_thickness,min_thickness,t_tol):
        
        from basic_setup import spectrum_list_generator,parameter_list_generator
      
        #test_setup=dyn_basic_setup(thickness=max_thickness,R_dir=self.input_data[2],T_dir=self.input_data[3])
        
        try: 
            test_setup=dyn_basic_setup(thickness=max_thickness,R_line=self.input_data[1],R_dir=self.input_data[2],
                                       T_dir=self.input_data[3],if_metal=self.if_metal)
        except:
            self.Invalid_Input.emit()
            time.sleep(1)
            return
       
        lamda_min,lamda_max=test_setup.get_lamda_range()
        if self.lamda_list[0]!='':
            lamda_min=float(self.lamda_list[0])
        if self.lamda_list[1]!='':
            lamda_max=float(self.lamda_list[1])
        
        
        ####log spacing
        #number_of_points = 10
        #film_thickness_list = logspace(log10(min_thickness),log10(max_thickness), number_of_points)
        
        ### single point thickness
        #film_thickness_list = [9.0, 10.5, 11.0, 11.5]
        
        
        dlamda_min = 5
        dlamda_max = 100
        delta_weight = 0.1/dlamda_min
        lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
        
        show_plots = True
        max_passes = 0 # limit the passes for a quicker scan, 0 means it will guess the amount needed.
        
       
        
        def fit_nk_f(lamda):                        #how the nk_guess start 
            return (1.0 + 0.1j) + 0.0*lamda
        
        
        ######## KK preconditioning is almost always a good idea, unless you have a metal, which it fails badly for
        film_thickness = max_thickness
        self.TB1.append('\n>>>>>>>> Film Thickness: %.3f<<<<<<<<\n'%film_thickness)
        parameter_list_generator.thickness = film_thickness
        data_directory = 'TRANK_nk_fit_%f_nm/'%film_thickness
        self.input_data[0]=max_thickness
        try_mkdir(data_directory)
        
        inputs=dict(nk_f_guess = fit_nk_f,
                fig=self.fig,
                fig2=self.fig2,
                fig3=self.fig3,
                TB1=self.TB1,
                spectrum_list_generator = spectrum_list_generator,
                parameter_list_generator = parameter_list_generator,
                lamda_min = lamda_min,
                lamda_max = lamda_max,
                dlamda_min = dlamda_min,
                dlamda_max = dlamda_max,
                max_passes = 1,
                delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'linear',
                adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.001,
                use_reducible_error = True,
                method='least_squares',
                KK_compliant = False,
                reuse_mode = False,
                zero_weight_extra_pass = False,
                data_directory = data_directory,
                verbose = True, make_plots = True, show_plots = show_plots,
                interpolate_to_fine_grid_at_end=False,
                QCoreApplication=QCoreApplication,
                Gui_mode=0,
                input_data=self.input_data,
                test_setup=test_setup)
        
        if self.set_list:
            inputs.update(dict(tolerance = float(self.set_list[0]),
                    adaptation_threshold_max = float(self.set_list[1]), adaptation_threshold_min = float(self.set_list[2])))
        if not self.if_metal: 
            inputs.update(dict(KK_compliant = True))
                  
        fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(**inputs)
    
        t_rms_dic={}
        #use Error Curves to find the minimum
        if not self.if_metal: 
            inputs.update(dict(KK_compliant = False, interpolation_type = 'cubic' ))
        while (max_thickness-min_thickness)>=3*t_tol :
            mid1=round((2*min_thickness+max_thickness)/3)
            mid2=round((min_thickness+2*max_thickness)/3)
            self.TB1.append ('>>>>>>>> Film Thickness: %.3f<<<<<<<<'%mid1)
            test_setup.parameter_list_generator(thickness = mid1)
            data_directory = 'TRANK_nk_fit_%f_nm/'%mid1
            try_mkdir(data_directory)
            self.input_data[0]=mid1
      
            inputs.update(dict(nk_f_guess = fit_nk_f,lamda_list = lamda_list,max_passes = 0,Gui_mode=3))
            rms_spectrum, fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(**inputs)
            t_rms_dic[mid1]=(sqrt( mean(array(rms_spectrum)**2)))
            
            self.TB1.append ('>>>>>>>> Film Thickness: %.3f<<<<<<<<'%mid2)
            test_setup.parameter_list_generator(thickness = mid2)
            data_directory = 'TRANK_nk_fit_%f_nm/'%mid2
            try_mkdir(data_directory)
            self.input_data[0]=mid2
          
            rms_spectrum, fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(**inputs)
            t_rms_dic[mid2]=(sqrt( mean(array(rms_spectrum)**2)))
            if t_rms_dic[mid1]<t_rms_dic[mid2]:
                max_thickness=mid2
            else:
                min_thickness=mid1
            print(t_rms_dic)
        
        thick=[]
        rms=[]
        rms_min=1
        thick_min=1
        for t in t_rms_dic.keys():
            thick.append(t)
            rms.append(t_rms_dic[t])
            if t_rms_dic[t]<rms_min:
                rms_min=t_rms_dic[t]
                thick_min=t
     
        self.Thickness_fit.emit((thick,rms,thick_min))
               
    def thickness_fit(self,if_metal=False,set_list=None):    
        time.sleep(1)
        show_plots = True
        data_directory = 'TRANK_nk_fit/'

        try_mkdir(data_directory)
        #if show_plots:
        #    from matplotlib.pylab import show    #use pyplot instead
    ###################### structure parameters
        from basic_setup import fit_nk_f,spectrum_list_generator,parameter_list_generator
        try: 
            test_setup=dyn_basic_setup(thickness=int(self.input_data[0]),R_line=self.input_data[1],R_dir=self.input_data[2],
                                       T_dir=self.input_data[3],if_metal=if_metal)
        except:
            self.Invalid_Input.emit()
            time.sleep(1)
            return
        #test_setup=dyn_basic_setup(thickness=int(self.input_data[0]),R_dir=self.input_data[2],
        #                                           T_dir=self.input_data[3])
        
    ###########
        from os import getcwd, walk, listdir
        from os.path import isfile

        #from numpy import arange, loadtxt, sqrt, mean, array
        if if_metal:
            dlamda_min=1
        else:
            dlamda_min=4         #how to determine the dlamda
        dlamda_max = 50
       
        lamda_max, lamda_min = test_setup.get_lamda_max_min()
        if self.lamda_list[0]!='':
            lamda_min=float(self.lamda_list[0])
        if self.lamda_list[1]!='':
            lamda_max=float(self.lamda_list[1])
        lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
        passes=0
        max_rms_cutoff = 5.0 #percentage points
        net_rms_cutoff = 1.0
        use_old_nk = False
        has_old_nk = False
        old_lamda = []
        if self.nk_dir:
        #if isfile(data_directory+'fit_nk_fine.txt') and isfile(data_directory+'fit_nk.txt'): # fine has the complete set
            old_data = loadtxt(self.nk_dir).T
        #if self.old_data!=None:
            #print('Found local data.')
            self.TB1.append('Found local data.')

            #old_data = loadtxt( data_directory+'fit_nk.txt').T
            #fit_nk_f =  functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')
            fit_nk_f =  functionize_nk_file(self.nk_dir, skiprows = 0, kind = 'cubic')   #use the user-defined file
            old_lamda = old_data[0]
            has_old_nk = True

            if has_old_nk:
               
               
                rms_spectrum = rms_error_spectrum(lamda_list = lamda_fine,
                                                  nk_f = fit_nk_f,
                                                  spectrum_list_generator = spectrum_list_generator,
                                                  parameter_list_generator = parameter_list_generator,
                                                  input_data=self.input_data,test_setup=test_setup)# different way in using spectrum_lamda_error

                net_rms = sqrt( mean( array(rms_spectrum)**2 ) ) * 100.0
                max_rms =     max(rms_spectrum) * 100.0

                #print('nk found! RMS (max): %.2f (%.2f)'%(net_rms, max_rms))
                self.TB1.append('nk found! RMS (max): %.2f (%.2f)'%(net_rms, max_rms))
                ylim = max_rms_cutoff - (max_rms_cutoff/net_rms_cutoff)*net_rms
                if max_rms  < ylim:
                    use_old_nk = True
                    passes = 2

    #use_old_nk = False
        if use_old_nk == False:
            if not if_metal:
                def fit_nk_f(lamda):
                    return 2.0+0.5j+0.0*lamda
            else:
                old_lamda = lamda_fine
    
                from numpy.random import rand
                min_n, max_n  = 0.0, 2.0
                min_k, max_k  = 0.0, 0.1
                rand_n = rand(lamda_fine.size)*(max_n - min_n) + min_n
                rand_k = rand(lamda_fine.size)*(max_k - min_k) + min_k
                fit_nk_f = extrap_c(lamda_fine, rand_n + 1.0j*rand_k)
    
                def fit_nk_f(lamda):
                        return 1.0+0.0*lamda
                   
            


        nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = old_lamda, title_string='Initial nk',
                show_nodes = False, show_plots = show_plots, fig=self.fig)
            
        self.fig.canvas.draw()
        self.error_pages=[]
        self.nk_pages=[]
        #if show_plots: self.canvas.draw()
        #QCoreApplication.processEvents()  # Force GUI to update
        inputs=dict(nk_f_guess = fit_nk_f,
                    fig=self.fig,
                    fig2=self.fig2,
                    fig3=self.fig3,
                    TB1=self.TB1,
                    spectrum_list_generator = spectrum_list_generator,
                    parameter_list_generator = parameter_list_generator,
                    lamda_min = lamda_min,
                    lamda_max = lamda_max,
                    dlamda_min = dlamda_min,
                    dlamda_max = dlamda_max,
                    max_passes=passes,
                    delta_weight = 0.02, tolerance = 1e-5, interpolation_type = 'cubic',
                    adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.001,
                    use_reducible_error = True,
                    method='least_squares',
                    KK_compliant = False,
                    reuse_mode = use_old_nk, lamda_list = old_lamda,
                    zero_weight_extra_pass = False,
                    verbose = True, make_plots = True, show_plots = show_plots,
                    #nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', 
                    #rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' ,
                    threads=cpu_count(),
                    QCoreApplication=QCoreApplication,
                    Gui_mode=1,
                    if_e=self.if_e,
                    input_data=self.input_data
                    ,test_setup=test_setup)
        
        if self.set_list:
            inputs.update(dict(tolerance = float(self.set_list[0]),
                    adaptation_threshold_max = float(self.set_list[1]), adaptation_threshold_min = float(self.set_list[2])))
        if not if_metal:
            delta_weight = 0.2/dlamda_min
            inputs.update(dict(KK_compliant = True,delta_weight=delta_weight,k_weight_fraction=0.25,interpolation_type = 'linear',
                               max_passes=2,Gui_mode=4,interpolate_to_fine_grid_at_end=False,reuse_mode = False))
            error_pages,nk_pages,fit_nk_f,lamda_list=error_adaptive_iterative_fit_spectra(**inputs)
            self.error_pages+=error_pages
            self.nk_pages+=nk_pages
            self.Output_axes.emit((self.error_pages,self.nk_pages))
            inputs.update(dict(nk_f_guess = fit_nk_f,lamda_list = lamda_list,KK_compliant = False,interpolation_type = 'cubic',
                               max_passes=1,Gui_mode=1,interpolate_to_fine_grid_at_end=True))
        error_pages,nk_pages=error_adaptive_iterative_fit_spectra(**inputs)
        self.error_pages+=error_pages
        self.nk_pages+=nk_pages
        self.Output_axes.emit((self.error_pages,self.nk_pages))
        

class Window(QWidget):
    
    
    
    def __init__(self):
        super(Window,self).__init__()
     
        #add widgets and set the layout
        self.setWindowTitle('TRANK_fit v1.1')
        self.creategridbox()   
        self.createpic() 
        self.createpic2()
        self.createpic3()
        #set the tab layout
        self.tabs=QTabWidget()
        self.tab1=QWidget()
        #self.tab2=QWidget()
        self.tabs.addTab(self.tab1, 'TRANK_fit')
        #self.tabs.addTab(self.tab2,'Scan_thickness')
        
        #main layout
        self.vbox=QVBoxLayout()
        self.vbox.addWidget(self.tabs)
        main_bot_hbox=QHBoxLayout()
        lb_help=QLabel('Press F1 for Help')
        main_bot_hbox.addWidget(lb_help)
        main_bot_hbox.addStretch(1)
        self.vbox.addLayout(main_bot_hbox)
        self.setLayout(self.vbox)
        
        #set layout for tab1
        hbox=QHBoxLayout()
        vbox=QVBoxLayout()
       
        hbox.addWidget(self.picbox1)
        hbox.addWidget(self.picbox2)
        vbox.addLayout(hbox)
        vbox.addStretch(1)
        vbox.addWidget(self.picbox3)
        #self.vbox.addStretch(1)
        vbox.addWidget(self.gridbox)
        self.tab1.setLayout(vbox)
        
        self.setGeometry(100, 50, 1090,1000)     
        self.center()                   #set the size of windows dynamicly
        
        
        self.T_dir='Transmission_10nm_CuAg_on_silica_substrate.txt'
        self.R_dir='Reflection_10nm_CuAg_on_silica_substrate.txt'
        self.nk_dir=None
      
       
    def creategridbox(self):  
        self.gridbox=QGroupBox()
        self.bot_tabs=QTabWidget()
        self.input_tab=QWidget()
        self.setting_tab=QWidget()
        self.output_tab=QWidget()
        self.bot_tabs.addTab(self.input_tab,'Input' )
        self.bot_tabs.addTab(self.setting_tab,'Setting')
        self.bot_tabs.addTab(self.output_tab,'Output')
        
        #set the layout for input tab
        grid1=QGridLayout()
        self.t_ch_box=QCheckBox('Scan Thickness')
        self.t_ch_box.setChecked(False)
        self.t_ch_box.stateChanged.connect(self.if_scan)
        self.m_ch_box=QCheckBox('Metal Material')
        self.m_ch_box.setChecked(True)
        self.e_ch_box=QCheckBox('Draw  ε')
        self.e_ch_box.setChecked(True)
        self.m_ch_box.stateChanged.connect(self.if_metal)
        
       
        self.lb_t=QLabel('Thickness  (nm)')
        self.le_t=QLineEdit('10')
        self.le_t.setFixedWidth(100)
        self.le_t.setToolTip('scan and fit thickness if leave it blank')
        self.lb_t_min=QLabel('min Thickness  (nm)')
        self.le_t_min=QLineEdit('')
        self.le_t_min.setFixedWidth(100)
        self.lb_t_tol=QLabel('Thickness Tolerance  (nm)')
        self.le_t_tol=QLineEdit('')
        self.le_t_tol.setFixedWidth(100)
        lb_R=QLabel('Reflectance  (%s)'%chr(248))       #use ASCII to get the icon
        self.le_R=QLineEdit('0 50 60')
      
        #hide sth
        self.lb_t_min.hide()
        self.le_t_min.hide()
        self.lb_t_tol.hide()
        self.le_t_tol.hide()
        
        lb_T=QLabel('Transmittance  (%s)'%chr(248))
        self.le_T=QLineEdit('0')
        lb_upload=QLabel('Upload Files')
        hbox=QHBoxLayout()
        btn_load_R=QPushButton('Reflectance')
        btn_load_R.setFixedWidth(120)
        btn_load_R.clicked.connect(self.readfile)
        btn_load_T=QPushButton('Transmittance')
        btn_load_T.setFixedWidth(120)
        btn_load_T.clicked.connect(self.readfile)
        btn_load_nk=QPushButton('nk file')
        btn_load_nk.setFixedWidth(120)
        btn_load_nk.clicked.connect(self.readfile)
       
        hbox.addWidget(btn_load_R,Qt.AlignLeft)
        hbox.addStretch(1)
        hbox.addWidget(btn_load_T,Qt.AlignHCenter)
        hbox.addStretch(1)
        hbox.addWidget(btn_load_nk,Qt.AlignRight)
        
        
        grid1.addWidget(self.t_ch_box,0,0,1,1)
        grid1.addWidget(self.m_ch_box,0,2,1,1)
        grid1.addWidget(self.e_ch_box,0,4,1,1)
        grid1.addWidget(self.lb_t,1,0,1,1)
        grid1.addWidget(self.le_t,1,1,1,1)
        grid1.addWidget(self.lb_t_min,1,2,1,1)
        grid1.addWidget(self.le_t_min,1,3,1,1)
        grid1.addWidget(self.lb_t_tol,1,4,1,1)
        grid1.addWidget(self.le_t_tol,1,5,1,1)
       
        grid1.addWidget(lb_R,2,0,1,1)
        grid1.addWidget(self.le_R,2,1,1,5)
    
        grid1.addWidget(lb_T,3,0,1,1)
        grid1.addWidget(self.le_T,3,1,1,5)
      
        grid1.addWidget(lb_upload,4,0,1,1)
        grid1.setRowStretch(3,20)
        grid1.setColumnMinimumWidth(1,200)
        grid1.setColumnMinimumWidth(3,200) 
        grid1.setColumnMinimumWidth(5,100)
        grid1.addLayout(hbox,5,0,1,6)
        
        self.input_tab.setLayout(grid1)
        
        
        #set the layout for setting tab
        grid3=QGridLayout()
        lb_tol=QLabel('tolerance')
        self.le_tol=QLineEdit('0.00001')
        self.le_tol.setFixedWidth(100)
        lb_thes_min=QLabel('threshold_min')
        self.le_thes_min=QLineEdit('0.001')
        self.le_thes_min.setFixedWidth(100)
        lb_thes_max=QLabel('threshold_max')
        self.le_thes_max=QLineEdit('0.01')
        self.le_thes_max.setFixedWidth(100)
        lb_lam_min=QLabel('lamda_min (nm)')
        self.le_lam_min=QLineEdit()
        self.le_lam_min.setFixedWidth(100)
        lb_lam_max=QLabel('lamda_max (nm)')
        self.le_lam_max=QLineEdit()
        self.le_lam_max.setFixedWidth(100)
        grid3.addWidget(lb_tol,2,0,1,1)
        grid3.addWidget(self.le_tol,2,1,1,1)
        grid3.addWidget(lb_thes_min,1,0,1,1)
        grid3.addWidget(self.le_thes_min,1,1,1,1)
        grid3.addWidget(lb_thes_max,1,2,1,1)
        grid3.addWidget(self.le_thes_max,1,3,1,1)
        grid3.addWidget(lb_lam_min,0,0,1,1)
        grid3.addWidget(self.le_lam_min,0,1,1,1)
        grid3.addWidget(lb_lam_max,0,2,1,1)
        grid3.addWidget(self.le_lam_max,0,3,1,1)
        grid3.setColumnMinimumWidth(1,200)
        grid3.setColumnMinimumWidth(3,200)
        grid3.setRowMinimumHeight(1,60)
        grid3.setColumnStretch(4,1)
        grid3.setRowStretch(2,1)
        self.setting_tab.setLayout(grid3)
        
        #set the layout for output tab
        grid2=QGridLayout()
        
        self.TB1=QTextEdit('')
        
        #set the height of QTextEdit
        text_format=QTextBlockFormat()
        text_format.setBottomMargin(0)
        text_format.setLineHeight(15,QTextBlockFormat.FixedHeight)
        text_cursor=self.TB1.textCursor()
        text_cursor.setBlockFormat(text_format)
        self.TB1.setTextCursor(text_cursor)
        
       
        self.TB1.verticalScrollBar().rangeChanged.connect(self.tb_scroll)
        #self.TB1.moveCursor(QTextCursor.End)
         
        
        self.bth=QPushButton('Help')
        self.bth.clicked.connect(self.helpfile)
        hbox2=QHBoxLayout()
        hbox2.addStretch(1)
        #hbox2.addWidget(self.bth)
       
        grid2.addWidget(self.TB1,0,0,10,2)
        grid2.addLayout(hbox2, 10,0,1,2)
        
        grid2.setVerticalSpacing(10)
        grid2.setHorizontalSpacing(15)
        self.output_tab.setLayout(grid2)
        
        vbox=QVBoxLayout()
        vbox.addWidget(self.bot_tabs)
        self.gridbox.setLayout(vbox)
       
        #self.gridbox.setWindowTitle('test')   
    
    def createpic(self):
        self.nk_pages=[]
        self.nk_page=1
        
        self.picbox1=QGroupBox()
        hbox1=QHBoxLayout()
        hbox2=QHBoxLayout()
   
        vbox=QVBoxLayout()
      
        self.fig = Figure(figsize=(4,3), dpi=100)
        self.canvas=FigureCanvas(self.fig)
        self.canvas.setFixedWidth(500)
        toolbar=MyCustomToolbar(self.canvas,self)
        
        self.btn_up1=QPushButton('Prev')
        self.btn_up1.clicked.connect(self.nk_pageup)
        self.btn_dn1=QPushButton('Next')
        self.btn_dn1.clicked.connect(self.nk_pagedn)

        self.btn_start=QPushButton('Start')
        self.btn_start.clicked.connect(self.TRANK_fit_spectra)
        hbox1.addWidget(toolbar)
        hbox1.addWidget(self.btn_start)
        hbox2.addWidget(self.btn_up1)
        hbox2.addStretch(1)
        hbox2.addWidget(self.btn_dn1)
        vbox.addLayout(hbox1)
        #vbox.addStretch(1)
        vbox.addWidget(self.canvas)
        vbox.addStretch(1)
        vbox.addLayout(hbox2)
        self.picbox1.setLayout(vbox)
    
    def createpic2(self):
        self.error_pages=[]
        self.error_page=1
        
        self.picbox2=QGroupBox()
        hbox1=QHBoxLayout()
        hbox2=QHBoxLayout()
        vbox=QVBoxLayout()
        self.fig2 = Figure(figsize=(4,3), dpi=100)
        self.canvas2=FigureCanvas(self.fig2)
        self.canvas2.setFixedWidth(500)
       
        toolbar=MyCustomToolbar(self.canvas2,self)
        self.btn_test=QPushButton('Try')
        self.btn_test.clicked.connect(self.testplot)
        self.btn_test.setVisible(False)
        self.btn_stop=QPushButton('Stop')
        self.btn_stop.setEnabled(False)
       
        self.btn_up2=QPushButton('Prev')
        self.btn_up2.clicked.connect(self.error_pageup)
        self.btn_dn2=QPushButton('Next')
        self.btn_dn2.clicked.connect(self.error_pagedn)
               
        hbox1.addWidget(toolbar)
        hbox1.addWidget(self.btn_test)
        hbox1.addWidget(self.btn_stop)
        hbox2.addWidget(self.btn_up2)
        hbox2.addStretch(1)
        hbox2.addWidget(self.btn_dn2)
        vbox.addLayout(hbox1)
        #vbox.addStretch(1)
        vbox.addWidget(self.canvas2)
        vbox.addStretch(1)
        vbox.addLayout(hbox2)
        self.picbox2.setLayout(vbox) 
    
    def createpic3(self):
       
        self.picbox3=QGroupBox()
     
        vbox=QVBoxLayout()
        self.fig3 = Figure(figsize=(8,2.3), dpi=100)
        self.canvas3=FigureCanvas(self.fig3)
        self.canvas3.setFixedWidth(1000)
       
        toolbar=MyCustomToolbar(self.canvas3,self)
       
        vbox.addWidget(toolbar)
       
        vbox.addStretch(1)
      
       
        #vbox.addStretch(1)
        vbox.addWidget(self.canvas3)
        vbox.addStretch(1)
        self.picbox3.setLayout(vbox) 
    
        
    def keyPressEvent(self,e):
        if e.key()==Qt.Key_F1:
             fo=open(webbrowser.open("read.txt"))
        
     
    def helpfile(self):
        fo=open(webbrowser.open("read.txt"))
        #f2.close()
    
    def readfile(self):
        try:
            sender=self.sender()
            filename=QFileDialog.getOpenFileName(self, 'Open File Dialog', 'C:',"Txt files(*.txt)")
            if(sender.text()=='Transmittance'):
                self.T_dir=filename[0]
            elif(sender.text()=='Reflectance'):
                self.R_dir=filename[0]
            elif(sender.text()=='nk file'):
                self.nk_dir=filename[0]
  
        except:
            QMessageBox.warning(self, 'Warning', 'Fail',QMessageBox.Yes)
            return
    
    def center(self):   #locate the windows the center of the screen
        index=QDesktopWidget().primaryScreen()
        screen = QDesktopWidget().availableGeometry(index)
        size = self.frameGeometry()
        self.move(screen.width()/2 - size.width()/1.5,    
        (screen.height() - size.height()) / 2)   
    
    
    def if_scan(self,checked):
        if checked:
            time.sleep(0.2)
            self.lb_t.setText('max Thickness  (nm)')
            self.lb_t_min.show()
            self.le_t_min.show()
            self.lb_t_tol.show()
            self.le_t_tol.show()
        else:
            self.lb_t_min.hide()
            self.le_t_min.hide()
            self.lb_t_tol.hide()
            self.le_t_tol.hide()
    
    def if_metal(self,checked):
        if checked:
            time.sleep(0.2)
            self.e_ch_box.show()
        else:
            self.e_ch_box.setChecked(False)
            self.e_ch_box.hide()    
    
    def tb_scroll(self,minV,maxV):
        self.TB1.verticalScrollBar().setValue(maxV)
       
    #These functions maybe combined but currently I cannot figure out how to pass parameters every time the button clicked  
    def error_pageup(self):
        if self.error_page>1:
            
            self.fig2.delaxes(self.error_pages[self.error_page-1])
            self.fig2.canvas.draw()
            self.error_page-=1
            time.sleep(2)
            
            self.fig2.add_axes(self.error_pages[self.error_page-1])
            self.fig2.canvas.draw()
        else:
            QMessageBox.warning(self, 'Warning', 'Reach first',QMessageBox.Yes)
            return
    def error_pagedn(self):
        if self.error_page<len(self.error_pages):
            
            self.fig2.delaxes(self.error_pages[self.error_page-1])
            self.fig2.canvas.draw()
            self.error_page+=1
            time.sleep(2)
            
            self.fig2.add_axes(self.error_pages[self.error_page-1])
            self.fig2.canvas.draw()
        else:
            QMessageBox.warning(self, 'Warning', 'Reach last',QMessageBox.Yes)
            return
    
    def nk_pageup(self):
        if self.nk_page>1:
            for i in range(0,2):
                self.fig.delaxes(self.nk_pages[self.nk_page-2+i])
            self.fig.canvas.draw()
            self.nk_page-=2
            time.sleep(2)
            for i in range(0,2):
                self.fig.add_axes(self.nk_pages[self.nk_page-2+i])
            self.fig.canvas.draw()
        else:
            QMessageBox.warning(self, 'Warning', 'Reach first',QMessageBox.Yes)
            return
               
    def nk_pagedn(self):
        if self.nk_page<len(self.nk_pages):
            for i in range(0,2):
                self.fig.delaxes(self.nk_pages[self.nk_page-2+i])
            self.fig.canvas.draw()
            self.nk_page+=2
            time.sleep(2)
            for i in range(0,2):
                self.fig.add_axes(self.nk_pages[self.nk_page-2+i])
            self.fig.canvas.draw()
        else:
            QMessageBox.warning(self, 'Warning', 'Reach last',QMessageBox.Yes)
            return
        
    def testplot(self,t_rms):
        thick=t_rms[0]
        rms=t_rms[1]
        thick_min=t_rms[2]
       
        QMessageBox.warning(self, 'Warning', 'Thickness is %s'%thick_min,QMessageBox.Yes)  
       
        fig, ax=plt.subplots()
        order=argsort(thick)
        line1,=ax.plot(array(thick)[order],array(rms)[order])
       
        plt.show()
        
        #self.fig.clear()
        
        #self.fig.add_axes(ax)
        #self.canvas2.draw()
        #self.fig.canvas.draw()
    
    def update_pages(self,pages):
        self.error_pages=pages[0]
        self.error_page=len(self.error_pages)
        self.nk_pages=pages[1]
        self.nk_page=len(self.nk_pages)
        self.btn_up1.setEnabled(True)
        self.btn_up2.setEnabled(True)
        self.btn_dn1.setEnabled(True)
        self.btn_dn2.setEnabled(True)
        
    def invalid_input(self):
        QMessageBox.warning(self, 'Warning', 'Invalid Input Files',QMessageBox.Yes)
        self.btn_stop.click()  
        
       
    def done(self):
        self.btn_start.setEnabled(True)
        self.btn_stop.setEnabled(False)
        
    
        
    def TRANK_fit_spectra(self):
        #spectrum_list_generator,parameter_list_generator=self.dyn_basic_setup()
       
        t_list=[self.le_t.text(),self.le_t_min.text(),self.le_t_tol.text()]
        set_list=[self.le_tol.text(),self.le_thes_min.text(),self.le_thes_min.text()]
        lamda_list=[self.le_lam_min.text(),self.le_lam_max.text()]
        self.get_thread=TRANK_fit_Thread(self.fig,self.fig2,self.fig3,self.TB1,t_list,self.le_R.text(),self.R_dir,self.T_dir,self.nk_dir,
                                         self.t_ch_box.isChecked(),self.m_ch_box.isChecked(),self.e_ch_box.isChecked(),set_list,lamda_list)
        self.get_thread.start()
        self.btn_start.setEnabled(False)
        self.btn_up1.setEnabled(False)
        self.btn_up2.setEnabled(False)
        self.btn_dn1.setEnabled(False)
        self.btn_dn2.setEnabled(False)
        #self.connect(self.get_thread,pyqtSignal('finished()'),self.test_done)
        
        self.btn_stop.setEnabled(True)
        time.sleep(2)
        self.bot_tabs.setCurrentWidget(self.output_tab)
        self.get_thread.finished.connect(self.done)
        self.get_thread.Output_axes.connect(self.update_pages)
        self.get_thread.Invalid_Input.connect(self.invalid_input)
        self.get_thread.Thickness_fit.connect(self.testplot)
        self.btn_stop.clicked.connect(self.get_thread.terminate)
   
        
if __name__=='__main__':
    #manager=MyManager()
    #manager.start()
    app=0
    app=QApplication(sys.argv)
    times=QFont('Times New Roman',10)
    app.setFont(times)    
    m=Window()
    m.show()
    sys.exit(app.exec_())
        