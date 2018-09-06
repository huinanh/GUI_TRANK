from PyQt5.QtWidgets import (QWidget, QApplication, QGroupBox,QTabWidget, QPushButton, 
QLabel, QHBoxLayout,  QVBoxLayout, QGridLayout, QFormLayout, QLineEdit, QTextEdit,QComboBox,QMessageBox,
    QDesktopWidget, QFileDialog,QTextEdit, QMessageBox)
from PyQt5.QtGui import QPixmap, QIcon, QFont, QGuiApplication, QTextBlockFormat, QTextCursor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar  
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sys
import os
from numpy import *
import pickle                         #used to store and reload fig
from TRANK import rms_error_spectrum, nk_plot, try_mkdir, functionize_nk_file, extrap_c,  error_adaptive_iterative_fit_spectra
import time
from PyQt5.Qt import QCoreApplication
import io
import webbrowser




class MyCustomToolbar(NavigationToolbar):
            toolitems = [t for t in NavigationToolbar.toolitems if
                 t[0] in ('Home', 'Back','Forward','Pan', 'Zoom','Save')]


class Window(QWidget):
    
    
    
    def __init__(self):
        super(Window,self).__init__()
        
        #add widgets and set the layout
        self.setWindowTitle('TRANK_fit v1.0')
        self.creategridbox()   
        self.createpic() 
        self.createpic2()
        
        #set the tab layout
        self.tabs=QTabWidget()
        self.tab1=QWidget()
        self.tab2=QWidget()
        self.tabs.addTab(self.tab1, 'TRANK_fit')
        self.tabs.addTab(self.tab2,'Scan_thickness')
        
        #main layout
        
        self.vbox=QVBoxLayout()
        self.vbox.addWidget(self.tabs)
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
        
        self.setGeometry(100, 50, 1000,750)     
        self.center()                   #set the size of windows dynamicly
       
        self.old_data=None

       
    def creategridbox(self):  
        self.gridbox=QGroupBox()
        grid=QGridLayout()
        
        self.TB1=QTextEdit('')
        
        #set the height of QTextEdit
        text_format=QTextBlockFormat()
        text_format.setBottomMargin(0)
        text_format.setLineHeight(15,QTextBlockFormat.FixedHeight)
        text_cursor=self.TB1.textCursor()
        text_cursor.setBlockFormat(text_format)
        self.TB1.setTextCursor(text_cursor)
        
        self.TB1.moveCursor(QTextCursor.End)
              
        self.LEf=QLineEdit('')
        self.btr=QPushButton('Upload')
        self.btr.clicked.connect(self.readfile)
        self.bth=QPushButton('Help')
        self.bth.clicked.connect(self.helpfile)
       
        
        #self.LEf.setFixedWidth(100)
        
        #set the layout
        grid.addWidget(self.TB1,0,0,10,2)
        grid.addWidget(self.LEf,10,0,1,1)
        grid.addWidget(self.btr,10,1,1,1)
        grid.addWidget(self.bth,11,1,1,1)
        
        grid.setVerticalSpacing(10)
        grid.setHorizontalSpacing(15)
        self.gridbox.setLayout(grid)
       
        #self.gridbox.setWindowTitle('test')   
    
    def createpic(self):
        self.nk_pages=[]
        self.nk_page=1
        
        self.picbox1=QGroupBox()
        hbox1=QHBoxLayout()
        hbox2=QHBoxLayout()
   
        vbox=QVBoxLayout()
      
        self.fig = Figure(figsize=(5,5), dpi=100)
        self.canvas=FigureCanvas(self.fig)
        toolbar=MyCustomToolbar(self.canvas,self)
        
        btup1=QPushButton('Prev')
        btup1.clicked.connect(self.nk_pageup)
        btdn1=QPushButton('Next')
        btdn1.clicked.connect(self.nk_pagedn)

        btnp=QPushButton('Start')
        btnp.clicked.connect(self.TRANK_fit_spectra)
        hbox1.addWidget(toolbar)
        hbox1.addWidget(btnp)
        hbox2.addWidget(btup1)
        hbox2.addStretch(1)
        hbox2.addWidget(btdn1)
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
        self.fig2 = Figure(figsize=(5,5), dpi=100)
        self.canvas2=FigureCanvas(self.fig2)
       
        toolbar=MyCustomToolbar(self.canvas2,self)
        btn=QPushButton('Try')
        #btn.clicked.connect(self.testplot)
       
        btup2=QPushButton('Prev')
        btup2.clicked.connect(self.error_pageup)
        btdn2=QPushButton('Next')
        btdn2.clicked.connect(self.error_pagedn)
               
        hbox1.addWidget(toolbar)
        hbox1.addWidget(btn)
        hbox2.addWidget(btup2)
        hbox2.addStretch(1)
        hbox2.addWidget(btdn2)
        vbox.addLayout(hbox1)
        #vbox.addStretch(1)
        vbox.addWidget(self.canvas2)
        vbox.addStretch(1)
        vbox.addLayout(hbox2)
        self.picbox2.setLayout(vbox) 
     
    def helpfile(self):
        fo=open(webbrowser.open("read.txt"))
        #f2.close()
    
    def readfile(self):
        try:
            self.nk_filename=QFileDialog.getOpenFileName(self, 'Open File Dialog', 'C:',"Txt files(*.txt)")
            self.old_data=loadtxt(self.nk_filename[0]).T
            self.LEf.setText(self.nk_filename[0])
       
        except:
            QMessageBox.warning(self, 'Warning', 'Fail',QMessageBox.Yes)
            return
    
    def center(self):   #locate the windows the center of the screen
        index=QDesktopWidget().primaryScreen()
        screen = QDesktopWidget().availableGeometry(index)
        size = self.frameGeometry()
        self.move(screen.width()/2 - size.width()/1.5,    
        (screen.height() - size.height()) / 2)   
       
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
    '''    
    def testplot(self):
        #self.fig.clear()
        with open('myplot.pkl','rb') as f:
             ax=pickle.load(f)
        
        if(self.count==0):
            self.fig.delaxes(self.nax)
            self.fig.delaxes(self.kax)
            self.fig.canvas.draw()
        else:
            self.fig.add_axes(self.nax)
            self.fig.add_axes(self.kax)
            self.fig.canvas.draw()
        self.count+=1
        #self.fig.add_axes(ax)
        #self.canvas2.draw()
        #self.fig.canvas.draw()
    '''   
        
    def TRANK_fit_spectra(self):
        #initialize the GUI
        self.fig.clear()
        self.fig.canvas.draw()
        self.fig2.clear()
        self.fig2.canvas.draw()
        self.TB1.clear()
        QCoreApplication.processEvents()
        
        
        show_plots = True
        data_directory = 'TRANK_nk_fit/'

        try_mkdir(data_directory)
        #if show_plots:
        #    from matplotlib.pylab import show    #use pyplot instead

    ###################### structure parameters
        from basic_setup  import fit_nk_f, spectrum_list_generator,   parameter_list_generator



    ###########
        from os import getcwd, walk, listdir
        from os.path import isfile

        #from numpy import arange, loadtxt, sqrt, mean, array
        
        dlamda_min = 1
        dlamda_max = 50
        lamda_min = 300
        lamda_max = 1000
        lamda_fine = arange(lamda_min, lamda_max + dlamda_min/2.0 , dlamda_min)



        max_rms_cutoff = 5.0 #percentage points
        net_rms_cutoff = 1.0
        use_old_nk = False
        has_old_nk = False
        old_lamda = []
        #if isfile(data_directory+'fit_nk_fine.txt') and isfile(data_directory+'fit_nk.txt'): # fine has the complete set
        if self.old_data!=None:
            #print('Found local data.')
            self.TB1.append('Found local data.')

            #old_data = loadtxt( data_directory+'fit_nk.txt').T
            #fit_nk_f =  functionize_nk_file(data_directory+'fit_nk.txt', skiprows = 0, kind = 'cubic')
            fit_nk_f =  functionize_nk_file(self.nk_filename[0], skiprows = 0, kind = 'cubic')   #use the user-defined file
            old_lamda = self.old_data[0]
            has_old_nk = True

            if has_old_nk:

                rms_spectrum = rms_error_spectrum(lamda_list = lamda_fine,
                                                  nk_f = fit_nk_f,
                                                  spectrum_list_generator = spectrum_list_generator,
                                                  parameter_list_generator = parameter_list_generator)

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
        QCoreApplication.processEvents()  # Force GUI to update
        from multiprocessing import cpu_count
        
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
                                                                            threads=cpu_count()-1,
                                                                            QCoreApplication=QCoreApplication,
                                                                            ifGui=True)
        
        self.error_page=len(self.error_pages) 
        self.nk_page=len(self.nk_pages)
    
            
if __name__=='__main__':
    app=0
    app=QApplication(sys.argv)
    times=QFont('Times New Roman',10)
    app.setFont(times)    
    m=Window()
    m.show()
    sys.exit(app.exec_())
        