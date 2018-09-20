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
    
    def __init__(self,fig,fig2,TB1,t_list,R_line,R_dir,T_dir,if_scan):
        QThread.__init__(self)
       
        self.fig=fig
        self.fig2=fig2
        self.TB1=TB1
        self.if_scan=if_scan
        self.t_list=t_list
        #from new_basic_setup  import dyn_basic_setup
        self.input_data=[int(t_list[0]),R_line,R_dir,T_dir]
        self.fig.clear()
        self.fig.canvas.draw()
        self.fig2.clear()
        self.fig2.canvas.draw()
        self.TB1.clear()
        #QCoreApplication.processEvents()
    
    ###### okay so this goofy stuff should allow us to update the thickness of our layer from other places
    
    def __del__(self):
        
        self.wait(3)
        
    def run(self):
        if self.if_scan :
            self.scan_thickness(int(self.t_list[0]),int(self.t_list[1]),int(self.t_list[2]),)
        else:
            self.thickness_fit()
    
    def scan_thickness(self,max_thickness,min_thickness,t_tol):
        
        from basic_setup import spectrum_list_generator,parameter_list_generator
        test_setup=dyn_basic_setup(thickness=max_thickness,R_dir=self.input_data[2],T_dir=self.input_data[3])
        lamda_min,lamda_max=test_setup.get_lamda_range()
     
        
        
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
    
        fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
                                nk_f_guess = fit_nk_f,
                                fig=self.fig,
                                fig2=self.fig2,
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
                                Gui_mode=0,
                                input_data=self.input_data,
                                test_setup=None )
    
        t_rms_dic={}
        #use Error Curves to find the minimum
        while max_thickness-min_thickness>1 :
            mid1=round((2*min_thickness+max_thickness)/3)
            mid2=round((min_thickness+2*max_thickness)/3)
            self.TB1.append ('>>>>>>>> Film Thickness: %.3f<<<<<<<<'%mid1)
            parameter_list_generator.thickness = mid1
            data_directory = 'TRANK_nk_fit_%f_nm/'%mid1
            try_mkdir(data_directory)
            self.input_data[0]=mid1
            rms_spectrum, fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
                                nk_f_guess = fit_nk_f,
                                fig=self.fig,
                                fig2=self.fig2,
                                TB1=self.TB1,
                                spectrum_list_generator = spectrum_list_generator,
                                parameter_list_generator = parameter_list_generator,
                                lamda_min = lamda_min,
                                lamda_max = lamda_max,
                                dlamda_min = dlamda_min,
                                dlamda_max = dlamda_max,
                                max_passes = max_passes,
                                lamda_list = lamda_list,## reuse the old lamda points
                                delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'cubic',
                                adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.0005,
                                use_reducible_error = True,
                                method='least_squares',
                                KK_compliant = False,
                                reuse_mode = True,
                                zero_weight_extra_pass = False,
                                data_directory = data_directory,
                                verbose = True, make_plots = True, show_plots = show_plots,
                                interpolate_to_fine_grid_at_end=False,
                                Gui_mode=3,
                                input_data=self.input_data,
                                test_setup=None )
            t_rms_dic[mid1]=(sqrt( mean(array(rms_spectrum)**2)))
            
            self.TB1.append ('>>>>>>>> Film Thickness: %.3f<<<<<<<<'%mid2)
            parameter_list_generator.thickness = mid2
            data_directory = 'TRANK_nk_fit_%f_nm/'%mid2
            try_mkdir(data_directory)
            self.input_data[0]=mid2
            rms_spectrum, fit_nk_f, lamda_list = error_adaptive_iterative_fit_spectra(
                                nk_f_guess = fit_nk_f,
                                fig=self.fig,
                                fig2=self.fig2,
                                TB1=self.TB1,
                                spectrum_list_generator = spectrum_list_generator,
                                parameter_list_generator = parameter_list_generator,
                                lamda_min = lamda_min,
                                lamda_max = lamda_max,
                                dlamda_min = dlamda_min,
                                dlamda_max = dlamda_max,
                                max_passes = max_passes,
                                lamda_list = lamda_list,## reuse the old lamda points
                                delta_weight = delta_weight, tolerance = 1e-5, interpolation_type = 'cubic',
                                adaptation_threshold_max = 0.05, adaptation_threshold_min = 0.0005,
                                use_reducible_error = True,
                                method='least_squares',
                                KK_compliant = False,
                                reuse_mode = True,
                                zero_weight_extra_pass = False,
                                data_directory = data_directory,
                                verbose = True, make_plots = True, show_plots = show_plots,
                                interpolate_to_fine_grid_at_end=False,
                                Gui_mode=3,
                                input_data=self.input_data,
                                test_setup=None )
            t_rms_dic[mid2]=(sqrt( mean(array(rms_spectrum)**2)))
            if t_rms_dic[mid1]<t_rms_dic[mid2]:
                max_thickness=mid2
            else:
                min_thickness=mid1
            print(t_rms_dic)
        print(t_rms_dic)
        fig, ax=plt.subplots()
        thick=[]
        rms=[]
        for t in t_rms_dic.keys():
            thick.append(t)
            rms.append(t_rms_dic[t])
        order=argsort(thick)
        line1,=ax.plot(array(thick)[order],array(rms)[order])
        plt.show()
               
    def thickness_fit(self):    
        show_plots = True
        data_directory = 'TRANK_nk_fit/'

        try_mkdir(data_directory)
        #if show_plots:
        #    from matplotlib.pylab import show    #use pyplot instead

    ###################### structure parameters
      
        from basic_setup import fit_nk_f,spectrum_list_generator,parameter_list_generator
        #test_setup=dyn_basic_setup(thickness=int(self.input_data[0]),R_dir=self.input_data[2],
        #                                           T_dir=self.input_data[3])
       
    ###########
        from os import getcwd, walk, listdir
        from os.path import isfile

        #from numpy import arange, loadtxt, sqrt, mean, array
        dlamda_min = 1
        dlamda_max = 50
        lamda_min = 300
        lamda_max = 1000
        lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)
        passes=0
        max_rms_cutoff = 5.0 #percentage points
        net_rms_cutoff = 1.0
        use_old_nk = False
        has_old_nk = False
        old_lamda = []
        if isfile(data_directory+'fit_nk_fine.txt') and isfile(data_directory+'fit_nk.txt'): # fine has the complete set
            old_data = loadtxt( data_directory+'fit_nk.txt').T
        #if self.old_data!=None:
            #print('Found local data.')
            self.TB1.append('Found local data.')

            #old_data = loadtxt( data_directory+'fit_nk.txt').T
            #fit_nk_f =  functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')
            fit_nk_f =  functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')   #use the user-defined file
            old_lamda = old_data[0]
            has_old_nk = True

            if has_old_nk:
               
                test_setup=dyn_basic_setup(thickness=int(self.input_data[0]),R_dir=self.input_data[2],T_dir=self.input_data[3])
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
            old_lamda = lamda_fine

            from numpy.random import rand
            min_n, max_n  = 0.0, 2.0
            min_k, max_k  = 0.0, 0.1
            rand_n = rand(lamda_fine.size)*(max_n - min_n) + min_n
            rand_k = rand(lamda_fine.size)*(max_k - min_k) + min_k
            fit_nk_f = extrap_c(lamda_fine, rand_n + 1.0j*rand_k)

            def fit_nk_f(lamda):
                return 1.0+0.0*lamda



        nk_plot(fit_nk_f, lamda_fine = lamda_fine, lamda_list = old_lamda, 
                #file_name = data_directory+'initial_nk.pdf',
                title_string='Initial nk',show_nodes = True, 
                show_plots = show_plots,fig=self.fig)
            
        self.fig.canvas.draw()
             
            #if show_plots: self.canvas.draw()
        #QCoreApplication.processEvents()  # Force GUI to update
       
        self.error_pages,self.nk_pages=error_adaptive_iterative_fit_spectra(nk_f_guess = fit_nk_f,
                                                                            fig=self.fig,
                                                                            fig2=self.fig2,
                                                                            TB1=self.TB1,
                                                                            spectrum_list_generator = spectrum_list_generator,
                                                                            parameter_list_generator = parameter_list_generator,
                                                                            lamda_min = lamda_min,
                                                                            lamda_max = lamda_max,
                                                                            dlamda_min = dlamda_min,
                                                                            dlamda_max = dlamda_max,
                                                                            max_passes=passes,
                                                                            delta_weight = 0.02, tolerance = 1e-5, interpolation_type = 'cubic',
                                                                            adaptation_threshold_max = 0.01, adaptation_threshold_min = 0.001,
                                                                            use_reducible_error = True,
                                                                            method='least_squares',
                                                                            KK_compliant = False,
                                                                            reuse_mode = use_old_nk, lamda_list = old_lamda,
                                                                            zero_weight_extra_pass = False,
                                                                            verbose = True, make_plots = True, show_plots = show_plots,
                                                                            #nk_spectrum_file_format = 'TRANK_nk_pass_%i.pdf', 
                                                                            #rms_spectrum_file_format = 'rms_spectrum_pass_%i.pdf' ,
                                                                            threads=cpu_count(),
                                                                            #QCoreApplication=QCoreApplication,
                                                                            Gui_mode=1,
                                                                            input_data=self.input_data
                                                                            ,test_setup=None)
        self.Output_axes.emit((self.error_pages,self.nk_pages))
        

class Window(QWidget):
    
    
    
    def __init__(self):
        super(Window,self).__init__()
     
        #add widgets and set the layout
        self.setWindowTitle('TRANK_fit v1.1')
        self.creategridbox()   
        self.createpic() 
        self.createpic2()
        
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
        #self.vbox.addStretch(1)
        vbox.addWidget(self.gridbox)
        self.tab1.setLayout(vbox)
        
        self.setGeometry(100, 50, 1090,810)     
        self.center()                   #set the size of windows dynamicly
        
        
        self.T_dir='Transmission_10nm_CuAg_on_silica_substrate.txt'
        self.R_dir='Reflection_10nm_CuAg_on_silica_substrate.txt'
        self.nk_dir=None
      
       
    def creategridbox(self):  
        self.gridbox=QGroupBox()
        self.bot_tabs=QTabWidget()
        self.input_tab=QWidget()
        self.output_tab=QWidget()
        self.bot_tabs.addTab(self.input_tab,'Input' )
        self.bot_tabs.addTab(self.output_tab,'Output')
        
        #set the layout for input tab
        grid1=QGridLayout()
        self.t_ch_box=QCheckBox('Scan Thickness')
        self.t_ch_box.setChecked(False)
        self.t_ch_box.stateChanged.connect(self.if_scan)
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
        lb_upload=QLabel('Uploda Files')
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
        
        btn_up1=QPushButton('Prev')
        btn_up1.clicked.connect(self.nk_pageup)
        btn_dn1=QPushButton('Next')
        btn_dn1.clicked.connect(self.nk_pagedn)

        self.btn_start=QPushButton('Start')
        self.btn_start.clicked.connect(self.TRANK_fit_spectra)
        hbox1.addWidget(toolbar)
        hbox1.addWidget(self.btn_start)
        hbox2.addWidget(btn_up1)
        hbox2.addStretch(1)
        hbox2.addWidget(btn_dn1)
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
        self.btn_stop=QPushButton('Stop')
        self.btn_stop.setEnabled(False)
       
        btn_up2=QPushButton('Prev')
        btn_up2.clicked.connect(self.error_pageup)
        btn_dn2=QPushButton('Next')
        btn_dn2.clicked.connect(self.error_pagedn)
               
        hbox1.addWidget(toolbar)
        hbox1.addWidget(self.btn_test)
        hbox1.addWidget(self.btn_stop)
        hbox2.addWidget(btn_up2)
        hbox2.addStretch(1)
        hbox2.addWidget(btn_dn2)
        vbox.addLayout(hbox1)
        #vbox.addStretch(1)
        vbox.addWidget(self.canvas2)
        vbox.addStretch(1)
        vbox.addLayout(hbox2)
        self.picbox2.setLayout(vbox) 
    
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
                self.R_dir=filename[0]
  
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
            time.sleep(0.5)
            self.lb_t.setText('max Thickness')
            self.lb_t_min.show()
            self.le_t_min.show()
            self.lb_t_tol.show()
            self.le_t_tol.show()
        else:
            self.lb_t_min.hide()
            self.le_t_min.hide()
            self.lb_t_tol.hide()
            self.le_t_tol.hide()
        
    
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
        
    def testplot(self):
        t_rms_dic={8.333333333333334: 0.77112439032284141, 11.666666666666666: 0.61334518148140205, 10.555555555555555: 0.62919731117350541, 
           12.777777777777779: 0.67174999134803504, 9.814814814814815: 0.61676329995077439, 11.296296296296298: 0.54184046164374322, 
           10.802469135802468: 0.51302339882446413, 11.790123456790125: 0.51786386479506286, 10.473251028806585: 0.48083555835042946, 
           11.131687242798355: 0.45437096390409232, 10.912208504801098: 0.45430411140263094, 11.351165980795612: 0.46581275052951787, 
           10.76588934613626: 0.46198520653026631, 11.058527663465936: 0.45990535584192221, 10.960981557689378: 0.45986089144968789, 
           11.156073769242495: 0.46123127076957487}
        from matplotlib import pyplot as plt
        fig, ax=plt.subplots()
        thick=[]
        rms=[]
        for t in t_rms_dic.keys():
            thick.append(t)
            rms.append(t_rms_dic[t])
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
        
    
    def done(self):
        self.btn_start.setEnabled(True)
        self.btn_stop.setEnabled(False)
        
    
        
    def TRANK_fit_spectra(self):
        #spectrum_list_generator,parameter_list_generator=self.dyn_basic_setup()
        self.bot_tabs.setCurrentWidget(self.output_tab)
        t_list=[self.le_t.text(),self.le_t_min.text(),self.le_t_tol.text()]
        self.get_thread=TRANK_fit_Thread(self.fig,self.fig2,self.TB1,t_list,self.le_R.text(),self.R_dir,self.T_dir,
                                         self.t_ch_box.isChecked())
        self.get_thread.start()
        self.btn_start.setEnabled(False)
        #self.connect(self.get_thread,pyqtSignal('finished()'),self.test_done)
        
        self.btn_stop.setEnabled(True)
        self.get_thread.finished.connect(self.done)
        self.get_thread.Output_axes.connect(self.update_pages)
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
        