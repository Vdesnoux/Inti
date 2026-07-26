
"""
INTI application interface Qt

Recontruction automatique d'une image ou de plusieurs images monochromatiques 
à partir d'un ou plusieurs fichiers video du spectre solaire
 
Author : Valérie Desnoux

Version post 6.4 en 2025

"""

import sys
import os
import yaml as yaml


from PySide6.QtUiTools import QUiLoader
from PySide6.QtWidgets import QApplication,QMainWindow,QFileDialog,QMessageBox, QListWidgetItem,QDockWidget,QWidget, QDialog
from PySide6.QtCore import QFile, QLocale, QIODevice, Qt, QPoint, QRect, QSettings, QTranslator, Signal, QObject
from PySide6 import QtGui
import pyqtgraph as pg
from pyqtgraph import ImageView, PlotWidget

from astropy.io import fits
import astropy.time
import math
import numpy as np
from PIL import Image
#import Inti_recon_bisector as sol
import Inti_recon as sol
from Inti_functions import *
import config as cfg
import stonyhurst as sth
import threading
import queue
import json
from skimage import exposure
#from scipy.signal import savgol_filter
import shutil
import re


import cv2 as cv2


from serfilesreader_vhd import Serfile

from datetime import datetime, UTC
import requests as rq
import webbrowser as web

import time # only for measure

import requests
import matplotlib.pyplot as plt #only for debug

# TODO : Angle P sur hélium



"""
version 7.0.4 WIP
- ajout detect_edge gestion crop dans detect_bord pour gerer fort tilt et diffusé
- si erreur detection inversions sur gong, affiche image gong original
- bug mode couronne y&,y2 etait en dur

Version 7.0.3
- bug fix sur nom fichier dp, cont avec decalage decimal
- affiche console valeur de la longueur d'onde selectionnée BASS2000
- ajout dans inti_finctions detect_bord test si bords trop proche pour eviter lignes sombres
- test si ligne très sombre dans image moyenne pour detection y1,y2 

Version 7.0.2
- gestion gros fichier SER
- sauve flag correction transversalium helium
- interdit +20%Hor si pas autocrop
- passage en decimal de cont et dop dans onglet doppler et free_shift1 et free_shift2

Version 7.0.1 
- replace ',' par '.' pour float (pas sure que cela suffise pour saveini et validators)
- gestion overflow fichier ser très gros

Version 7.0
- correction bug angle P annuler si flags weak et trans
- mise en route de check version
- suppression saut de ligne sur message de creation inti.yaml

Version 6.9+ 
- correction bug trad 

Version 6.9 - antibes 5 sept 2025
- ajout anti-aliasing polynomial
- recrée un fichier log au nom de l'image enregistrée
- supprime les accents dans les logs

Version 6.8 - paris 22 aout 2025
- ajout des params calculette shg autre
- save fits doppler optionel
- save settings fenetre gong
- corrige noms fits avec cont et png cont avec dp
- ajoute log valeurs de shift
- corrige trame_profil ecart de 1 dans auto
- recon avec double passe flat et meme seuils mask
- recon lecture ram
- saveas img avec seuils

Version 6.6g -6.7 - post ohp 25
- supprime sauvegarde du tab current
- corrige erreur lancement dernier ser
- amélioration rapidité inti_recon
- ajout fonction blink et ajustement taille dans gong
- bug inversion slwidth et slheight, valeur scampix et keyword shgangle
- Ha2cb passe a 6561.231
- mode doppler optionnel
- ajout mode couleur rouge

Version 6.6f - 15 juiller 2025
- compil trad
- gestion seuil et format png protus

Version 6.6e - 15 juillet 2025
- mode couronne zone calcul poly fixe
- ajuste polynome hors outliers
- ajout deux mots clefs BASS2000

Version 6.6+ - 28 juin 2025
- mode CLI un fichier avec option lance ou pas le traitement
- toujours onglet general au lancement
- angle P non auto par défaut
- ajout extension de 20% horizontal sur autocrop
- amelioration autocrop disque très decentré mais entier
- ajout save
- ajout edition crop value
- evite saturation sur clahe
- garde les seuils de crop 
- trame image entiere
- amélioration interpolation 4 points
- gestion nom _cont avec shift non nul
- BASS2000 ajouts de mots-clefs
- test deg 3 poly: aucune diff, reste en deg2, idem CB test ISIS
- accepte shift decimal

Version 6.5c - 16 avril 2025
- Inti avec nouvelle interface en Qt
- prise en compte de l'angle P sur une serie
- Correction bug formule dans angle_P_B0 et rotation de carrington
- test batch fichiers BASS2000
- correction kept_frame dans inti_recon et frame_index
- ajout d'une image mean_start de 20 trames autour de l'apparition du spectre en mode free only
- test image 8 bits
- detect_edge ajout test +:- clip et gestion hauteur de crop
- redirection flag pour thread et print
- image gong et _disk cote à cote
- choix images à enregistrer
- image cam_height avec padding pour angle P sur disque partiel
- case à cocher angle P auto ou pas dans section avancées
- bouton return par defaut declenche le traitement
- gestion inversions dans gong avec detection inversions auto
- detection du max profil sum de trame mean_start
- choix limité entre mean et mean_start et ne declenche pas calc max si mean
- bouton effacer contenu console 
- gestion constante poly négative
- ajout du gif en plus du mp4
- affiche fichier trame sans Qfiledialog et sans selecteur meanstart en mode magnet
- gestion erreurs si pas de fichier et test Flags au lancement
- ajout mask sur non-unif basse freq correction
- corrige reference tilt
- pas d'angle P pour disque partiel et helium transversallium ou magnet
- pas angle P auto au départ, toujours premier onglet, no more edition polynome
"""


# TODO : reduc bruit ajustable
# TODO : stop dans thread process
# TODO : memoriser les positions des fenetres images si double écran



original_stdout = sys.stdout

class main_wnd_UI(QMainWindow) :
    def __init__(self, fichier = None, flag = None):
    #def __init__(self, parent=None, fichier = None):
       
        #super().__init__(parent)
        super(main_wnd_UI, self).__init__()

        self.version ="7.0.3a"
        #iv = ImageView # force le load d'ImageView avant QUILoader
        # ne change rien...
        
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        loader.registerCustomWidget(ImageView)
        loader.registerCustomWidget(PlotWidget)
        ui_file_name=resource_path('inti_qt.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()
        
        
        self.bass_entete ={'observer':'', 'instru':'','wave_label':'Manual', 'lat':'0', 'long':'0',
                           'camera':'', 'spectro':'SOLEX',"pixel":'0','objcollim':'80', 'objcam':'125', 'focinstru':'0',
                            'diaminstru':'0', 'waveID':0, "binID":0, 'reseau':'2400', 'ordre':'1',
                            'contact':'', 'angle':'34','fentelong':'4.5', 'fentelarg':'10',
                            'diaph' : '0', 'dfilter' : ''}
        

        self.read_settings()
               
        # set icon application
        self.ui.setWindowIcon(QtGui.QIcon(resource_path("inti_logo.png")))
        
        # connect window close button to closeEvent
        app.aboutToQuit.connect(self.close)
        
        debug_noredir = False # True = mode debug ne redirigie pas vers le textEdit console
        
        if debug_noredir : 
            pass
        else :
            # redirection
            self.stdout_redirector = StdoutRedirector()
            self.stdout_redirector.new_text.connect(self.add_text)
            sys.stdout = self.stdout_redirector
        
        # geometry écrans
        screen_geom = QtGui.QGuiApplication.primaryScreen().availableGeometry()
        self.dpr = QtGui.QGuiApplication.primaryScreen().devicePixelRatio()
        self.myscreen_w = int(screen_geom.right()*self.dpr)
        self.myscreen_h = int(screen_geom.bottom()*self.dpr)
        
        # gestion langue

        if self.langue=='EN' :
           self.ui.lang_button.setText('EN')
           cfg.LG = 1
        else :
           self.ui.lang_button.setText('FR')
           cfg.LG = 2

        
        # init param
        self.myROI=[]
        self.pattern=''
        
        self.serfiles=[]
        self.Flags={}
        self.Shift=[]
        self.previous_serfile=''
        self.ang_tilt = 0 
        self.ratio_fixe = 1
        #self.list_wave=[['Manual','Ha','Ha2cb','Cah','Cah1v','Cak','Cak1v','HeID3'],[0,6562.762,6561.432,3968.469,3966.968,3933.663,3932.163,5877.3]]
        self.list_wave=[['Manual','Ha','Ha2cb','Cah','Cah1v','Cak','Cak1v','HeID3'],[0,6562.762,6561.232,3968.469,3966.968,3933.663,3932.163,5877.3]]
        self.list_binning=['1x1','2x2','3x3','4x4']
        self.solar_dict={}
        self.img_list =[]
        self.subwindows=[]
        self.param=[0,0,0,0]
        
        
        
        self.flag_grid_disp = True 
        self.grid_color = 'yellow'
        self.flag_ser = True
        
        self.ii=1
        
        # inti read data from ini file and update UI
        #--------------------------------------------------------------------
        self.read_ini()
        
        # connecte les signaux
        #--------------------------------------------------------------------
        self.ui.inti_open_btn.clicked.connect(self.inti_open)
        self.ui.inti_last_btn.clicked.connect(self.inti_open_last)
        
        self.ui.Exit.clicked.connect(self.exit_clicked)
        self.ui.lang_button.clicked.connect(self.lang_switch_clicked)
        self.ui.version_label.setText("Version : "+self.version)
        
        # Create a keyboard shortcut Ctrl+M
        #shortcut = QtGui.QKeySequence(Qt.CTRL + Qt.Key_M)
        """
        shortcut = QtGui.QKeySequence("Escape")
        self.shortcut = QtGui.QShortcut(shortcut, self.ui)
        self.shortcut.activated.connect(self.img_allclose)
        """
        
        #self.ui.connect(QtGui.QShortcut(QtGui.QKeySequence("Escape"), self.ui),SIGNAL('activated()'), self.img_allclose)
        # Create a keyboard shortcut Ctrl+M        #shortcut = QtGui.QKeySequence(Qt.CTRL + Qt.Key_M)
        shortcut_e = QtGui.QKeySequence("Return")
        self.shortcut_e = QtGui.QShortcut(shortcut_e, self.ui)
        self.shortcut_e.activated.connect(self.inti_go)
        self.ui.inti_go_btn.clicked.connect (self.inti_go)
        
        validator_2dec = DotOnlyDoubleValidator()
        validator_2dec.setLocale(QLocale(QLocale.English))
        validator_2dec.setDecimals(2)           # Nombre max de décimales autorisées
        
        
        # tab general
        # -------------------------------------------------------------------
        self.ui.inti_section_btn.clicked.connect(self.inti_section_show)
        self.ui.inti_database_btn.clicked.connect(self.inti_database_show)
        self.ui.inti_orientation_btn.clicked.connect(self.inti_orientation_show)
        self.ui.dock_main.hide()
        self.ui.inti_calc_btn.clicked.connect (self.inti_calc)
        self.ui.inti_go_btn.setDefault(True)
        self.ui.inti_go_btn.setAutoDefault(True)
        self.ui.console_clear_btn.clicked.connect(self.console_clear)
        self.ui.inti_shift_text.setValidator(validator_2dec)
        
        # database
        self.ui.db_wave_text.addItems(self.list_wave[0])
        self.ui.db_entete_btn.clicked.connect(self.db_edit_entete)
        self.ui.db_lat_text.setValidator(validator_2dec)
        self.ui.db_long_text.setValidator(validator_2dec)
        
        # orientation
        self.ui.ori_gong_btn.clicked.connect(self.ori_gong)
        self.ui.ori_angP_btn.clicked.connect(self.ori_angP)
        self.ui.ori_angP_text.setValidator(validator_2dec)
        
        # avancées
        self.ui.ori_grid_format_btn.clicked.connect(self.ori_grid_format)
        self.ui.cfg_files_btn.clicked.connect(self.cfg_files_to_save)
        self.ui.cfg_cropmanu_btn.clicked.connect(self.crop_manu)
        self.ui.section_autocrop_chk.clicked.connect(self.crop_clicked)
        self.ui.db_wave_text.currentIndexChanged.connect(self.db_print_wave)
        
        # doppler
        # -------------------------------------------------------------------
        self.ui.dop_shift_text.setText(str(self.dec_pix_dop))
        self.ui.dop_conti_shift_text.setText(str(self.dec_pix_cont))
        self.ui.dop_profil_btn.clicked.connect(self.dop_profil)
        self.ui.dop_calc_btn.clicked.connect(self.inti_calc)
        self.ui.dop_offset_text.setValidator(validator_2dec)
        self.ui.dop_conti_shift_text.setValidator(validator_2dec)
        
        # free
        # -------------------------------------------------------------------
        self.ui.free_trame_btn.clicked.connect(self.trame_mean_img)
        self.ui.inti_calc_btn_3.clicked.connect (self.inti_calc)
        self.ui.free_corona_chk.stateChanged.connect (self.corona_clicked)
        self.ui.free_shift_text.setValidator(validator_2dec)
        self.ui.free_shift1_text.setValidator(validator_2dec)
        self.ui.free_shift2_text.setValidator(validator_2dec)
        
        # magnet
        # -------------------------------------------------------------------
        self.ui.magnet_trame_btn.clicked.connect(self.trame_mean_img)
        self.ui.inti_calc_btn_2.clicked.connect (self.inti_calc)
        self.ui.magnet_shift_text.setValidator(validator_2dec)
        
        # on teste version du web
        self.check_version()
        
        # force antialiasing option de pyqtgraph
        pg.setConfigOptions(antialias=True)
        
        # on affiche version
        print('INTI Version : ' + str(self.version))
        
        # finalement on demarre toujours sur tab general
        self.current_tab= int(0)
        self.ui.tab_main.setCurrentIndex(int(0))
        
        # recupere si un fichier est passé
        #print((f"Fichier reçu : {fichier if fichier else 'Aucun'}"))
        #print((f"Flag reçu : {flag if flag else 'None'}"))
        if fichier != None :
            self.working_dir= self.get_dirpath(fichier)
            self.ui.inti_dir_lbl.setText(self.working_dir)
            self.current_tab= 0
            if flag != None :
                self.inti_go_one(fichier, flag)
            else :
                self.inti_go_one(fichier, False)
            
            
        self.serfiles=[]
    
    def show(self) :
        self.ui.show()
        self.ui.dock_console.show()
        self.get_Flags_tab()
       
    
    def closeEvent (self,event):
        sys.stdout=original_stdout
        QApplication.restoreOverrideCursor()
        try :
            self.save_ini()
            self.write_settings()
        except :
            pass

        
    def exit_clicked(self) :
        self.save_ini()
        self.write_settings()
        sys.stdout=original_stdout
        print('exit')
        for dock in self.ui.findChildren(QDockWidget):
            
            if dock.isFloating() :
                #print(dock.windowTitle())
                dock.close()
                self.ui.dock_console.hide()
        
        QApplication.instance().quit()
        app.quit()
        
    def register(self, window):
        self.subwindows.append(window)
        print('register')
        
 
    def lang_switch_clicked (self):
        if self.ui.lang_button.text()=='FR' :
            self.ui.lang_button.setText('EN')
            self.langue='EN'
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("Restart application to change language")
            msg.setWindowTitle("Message")
            msg.exec()
        else:
            self.ui.lang_button.setText('FR')
            self.langue='FR'
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setText("Redémarrer l'application pour changer la langue")
            msg.setWindowTitle("Message")
            msg.exec()

        
    
    
    # fonction non active
    def check_version (self):
        # contact site inti pour nouvelle version
        try :
            reponse = requests.get('http://valerie.desnoux.free.fr/inti/inti_version.html', timeout=10)
            if reponse.status_code == 200 :
                html_text=reponse.text
                pos=html_text.find('Version =')
                v=html_text[pos+10:pos+15]
                print('INTI web version : '+ v)
                if v != self.version :
                    box=False
                    if box :
                        msg = QMessageBox()
                        msg.setIcon(QMessageBox.Icon.Information)
                        msg.setText(self.tr("Une nouvelle Version est disponible"))
                        msg.setWindowTitle("Inti")
                        msg.exec()
                    # ou alors on met en rouge la version 
                    self.ui.version_label.setStyleSheet('color: red')
        except:
            pass
        
    def read_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        self.ui.restoreGeometry(settings.value("MainWindow/geometry"))
        self.ui.restoreState(settings.value("MainWindow/windowState"))
        
        if settings.value("App/lang") is not None :
            self.langue=settings.value("App/lang")
        else :
            self.langue='FR'
            
        if settings.value("Shg/foc") is not None :
            self.shg_foc_camera_other=settings.value("Shg/foc")
        else :
            self.shg_foc_camera_other='0'
        
        if settings.value("Shg/angle") is not None :
            self.shg_angle_other=settings.value("Shg/angle")
        else :
            self.shg_angle_other='0'
        
        if settings.value("App/save") is not None :
            valeur_json= settings.value("App/save", "{}")
            self.files_to_save = json.loads(valeur_json)
        else :
            self.files_to_save={'disk': True, 'raw':True, 'color' : True, "inv" : False, "mix" :False, "diff" : False, "sum": False, "x0x1x2" : False, "dop":False}
        
        if settings.value("App/db") is not None :
            valeur_json= settings.value("App/db", "{}")
            self.bass_entete = json.loads(valeur_json)
        else :
            self.bass_entete ={'observer':'', 'instru':'','wave_label':'Manual', 'lat':'0', 'long':'0',
                               'camera':'', 'spectro':'SOLEX',"pixel":'0','objcollim':'80', 'objcam':'125', 'focinstru':'0',
                                'diaminstru':'0', 'waveID':0, "binID":0, 'reseau':'2400', 'ordre':'1',
                                'contact':'', 'angle':'24','fentelong':'4.5', 'fentelarg':'10',
                                'diaph' : '0', 'dfilter' : ''}
        
        
        if self.ui.dock_console.isFloating() :
            dock_geometry = self.ui.dock_console.geometry()
            visible = any (screen.geometry().intersects(dock_geometry) for screen in QtGui.QGuiApplication.screens())
            if not visible :
                screen_geom = QtGui.QGuiApplication.primaryScreen().availableGeometry()
                new_geom = QRect (screen_geom.x()+50, screen_geom.y()+50,dock_geometry.width(), dock_geometry.height())
                self.ui.dock_console.setGeometry(new_geom)
            
    
    
    def write_settings(self) :
        # force docks non floating
        #self.ui.dock_console.setFloating(False)
        #self.addDockWidget(Qt.BottomDockWidgetArea, self.ui.dock_console)
        self.bass_entete['observer'] = self.ui.db_observer_text.text()
        self.bass_entete['instru'] = self.ui.db_instru_text.text()
        self.bass_entete['lat'] = self.ui.db_lat_text.text()
        self.bass_entete['long'] = self.ui.db_long_text.text()
        self.bass_entete['wavelength'] = self.list_wave[1][self.ui.db_wave_text.currentIndex()]
        self.bass_entete['wave_label'] = self.ui.db_wave_text.currentText()
        self.bass_entete['version']= "INTI "+self.version
        self.bass_entete['waveID']=self.ui.db_wave_text.currentIndex()
        self.bass_entete['binID']=self.ui.db_binning_combo.currentIndex()
        
        # sauve settings
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("MainWindow/geometry", self.ui.saveGeometry())
        settings.setValue("MainWindow/windowState", self.ui.saveState())
        settings.setValue("App/lang", self.langue)
        #settings.setValue("App/tab_index", self.ui.tab_main.currentIndex())
        settings.setValue("App/save", json.dumps(self.files_to_save))
        settings.setValue("App/db", json.dumps(self.bass_entete))
        #settings.setValue("Shg/foc", self.shg_foc_camera_other)
        #settings.setValue("Shg/angle", self.shg_angle_other)
        
    def read_ini(self) :
        # recuperation des parametres de configuration memorisé au dernier traitement
        # --------------------------------------------------------------------------------------

        # inti.yaml is a bootstart file to read last directory used by app
        # this file is stored in the module directory
        self.my_ini=data_path('inti.yaml')
        my_dictini={'directory':'', 'dec doppler':3, 'dec cont':15, 
                    'size_pix_cam':'0.0024', 'bin_cam':'1',
                    'poly_slit_a':0, "poly_slit_b":0,'poly_slit_c':0, 
                    'ang_tilt':0, 'ratio_sysx':1,
                    'free_autpoly':0, 'zee_autpoly':0,'poly_free_a':0,'poly_free_b':0,'poly_free_c':0,
                    'pos_free_blue':0, 'pos_free_red':0, 'free_shift':0,
                    'force_free_magn': False,
                    'win_posx':300, 'win_posy':200, 'screen_scale':0,'observer':'', 'instru':'','site_long':0, 'site_lat':0,
                    'angle P':0,'contact':'','wavelength':0, 'wave_label':'Manuel', 'inversion NS':0, 'inversion EW':0,
                    'autocrop':1,'ext20pct':0, 'pos fente min':0, 'pos fente max':0,'trame couronne1':-25, 'trame couronne2':-5,
                    "zeeman_shift":0, "reduction bruit":0, 'grid disk':'on', 'lang' :'FR', 'correction He':0, 'angP auto': 0,
                    "zoom":0, "cont_only":0, "dop_color_force":0}
    
        poly=[]
        Flags=self.Flags
        
        try:
            with open(self.my_ini, "r") as f1:
                my_dictini = yaml.safe_load(f1)
        except:
           
           print('Création de inti.yaml : '+ self.my_ini)
           
        
        self.working_dir=my_dictini['directory']
        self.dec_pix_dop=int(my_dictini['dec doppler'])
        self.dec_pix_cont=int(my_dictini['dec cont'])
        self.saved_tilt=to_float(str(my_dictini['ang_tilt']))
        self.saved_angP=to_float(str(my_dictini['angle P']))
        saved_ratio=to_float(str(my_dictini['ratio_sysx']))
        if saved_ratio == 0 :
            saved_ratio= 1
        poly.append(float(my_dictini['poly_free_a']))
        poly.append(float(my_dictini['poly_free_b']))
        poly.append(float(my_dictini['poly_free_c']))
        self.pos_free_blue=int(my_dictini['pos_free_blue'])
        self.pos_free_red=int(my_dictini['pos_free_red'])
        
        # ces deux paramètres ne sont pas mémorisés dans le fichier ini
        Flags['DOPFLIP']=0
        Flags["SAVEPOLY"]=0
        
        if 'cont_only' in my_dictini :
            Flags['Contonly']=my_dictini['cont_only']
        else :
            Flags['Contonly']=0
            
        if 'dop_color_force' in my_dictini :
            Flags['dop_color_force']=my_dictini['dop_color_force']
        else :
            Flags['dop_color_force']=0
            
        if 'ext20pc' in my_dictini :
            Flags['h20percent']=my_dictini['ext20pc']
        else :
            Flags['h20percent']=0
        
        if 'correction He' in my_dictini :
            Flags["FREE_TRANS"]=my_dictini['correction He']
        else :  
            Flags["FREE_TRANS"]=1
            
        if 'inversion EW' in my_dictini:
            Flags['FLIPRA']=my_dictini['inversion EW']
            Flags['FLIPNS']=my_dictini['inversion NS']
        else:
            Flags['FLIPRA']=0
            Flags['FLIPNS']=0
            
        if 'force_free_magn' in my_dictini :
            Flags['FORCE_FREE_MAGN']  = my_dictini['force_free_magn']
        else:
            Flags['FORCE_FREE_MAGN'] = False
        
        if 'size_pix_cam' not in my_dictini:
            # si pas dans fichier ini
            my_dictini['size_pix_cam']='0.0024'
            self.size_pix_cam='0.0024'
        else :
            self.size_pix_cam=my_dictini['size_pix_cam']
        
        if 'bin_cam' not in my_dictini:
            my_dictini['bin_cam']='1'
            self.bin_cam='1'
        else:
            self.bin_cam=my_dictini['bin_cam']
         
        if 'autocrop' not in my_dictini:
            my_dictini['autocrop']='1'
            Flags['Autocrop']=1
        else:
            Flags['Autocrop']=my_dictini['autocrop']
            
        if 'zoom' not in my_dictini:
            my_dictini['zoom']='0'
            Flags['zoom']=0
        else:
            Flags['zoom']=my_dictini['zoom']
            
        if 'free_autopoly' not in my_dictini:
            my_dictini['free_autopoly']=0
            Flags['FREE_AUTOPOLY']=0
        else:
            Flags['FREE_AUTOPOLY']=my_dictini['free_autopoly']
        
        if 'zee_autopoly' not in my_dictini:
            my_dictini['zee_autopoly']=0
            Flags['ZEE_AUTOPOLY']=0
        else:
            Flags['ZEE_AUTOPOLY']=my_dictini['zee_autopoly']
            
        if 'free_shift' not in my_dictini:
            my_dictini['free_shift'] = 0 #decalage en mode zeeman
            self.free_shift = 0
        else :
            self.free_shift = my_dictini['free_shift']
            
        if 'correction_He' not in my_dictini:
            my_dictini['correction He']=0
            self.Flags["FREE_TRANS"] = 0 
        else :
            self.Flags["FREE_TRANS"] = my_dictini['correction He']
            
            
        if 'zeeman_shift' not in my_dictini:
            my_dictini['zeeman_shift'] = 0 #decalage en mode zeeman
            self.zee_shift = 0
        else :
            self.zee_shift = my_dictini['zeeman_shift']
            
            
        if 'reduction_bruit' not in my_dictini:
            # si pas dans fichier ini
            my_dictini['reduction_bruit']=0
            Flags['NOISEREDUC'] = 0
        else :
            Flags['NOISEREDUC']=my_dictini['reduction_bruit']
            
        if 'angP auto' not in my_dictini:
            # si pas dans fichier ini
            my_dictini['angP auto']=1
            self.ui.section_angP_auto_chk.setChecked(True)
        else :
            self.ui.section_angP_auto_chk.setChecked(my_dictini['angP auto'])
            
        if 'grid disk' not in my_dictini:
            # si pas dans fichier ini
            my_dictini['grid disk']='on'
            self.grid_on=True
            Flags['Grid'] = 1
        else :
            if my_dictini['grid disk'] == 'on' :
                self.grid_on=True
                Flags['Grid'] = 1
            else :
                self.grid_on=False
                Flags['Grid'] = 0
        
        # param cas difficiles
        if 'trame couronne1' not in my_dictini:
            my_dictini['trame couronne1']='-25'
        if 'trame couronne2' not in my_dictini:
            my_dictini['trame couronne2']='-5'
        if 'pos fente min' not in my_dictini:
            my_dictini['pos fente min']='0'
        if 'pos fente max' not in my_dictini:
            my_dictini['pos fente max']='0'
            
        #gestion langue dans inti.yaml
        if 'lang' not in my_dictini :
            my_dictini['lang']='FR'
            LG_str='FR'
            cfg.LG=1
        else :
            LG_str = my_dictini['lang']
            if LG_str == 'FR' : cfg.LG=1
            if LG_str !='FR' : cfg.LG=2
            
                    
        
        self.param=[my_dictini['pos fente min'],my_dictini['pos fente max'],0,0,my_dictini['trame couronne1'],my_dictini['trame couronne2']]

        self.Flags=Flags
        self.my_dictini=my_dictini
        self.polynome = poly
        
        # met à jour l'interface UI
        self.ui.section_reduc_chk.setChecked(self.Flags['NOISEREDUC'])
        self.ui.section_autocrop_chk.setChecked(self.Flags['Autocrop'])
        self.ui.inti_h20_radio.setChecked(self.Flags['h20percent'])
        self.ui.section_zoom_chk.setChecked(self.Flags['zoom'])
        self.ui.free_tilt_chk.setChecked(self.Flags["FORCE_FREE_MAGN"])
        self.ui.magnet_tilt_chk.setChecked(self.Flags["FORCE_FREE_MAGN"])
        self.ui.section_grid_chk.setChecked( self.Flags['Grid'])
        self.ui.ori_inv_EW_chk.setChecked(self.Flags['FLIPRA'])
        self.ui.ori_inv_NS_chk.setChecked(self.Flags['FLIPNS'])
        self.ui.free_poly_auto_chk.setChecked(self.Flags['FREE_AUTOPOLY'])
        self.ui.magnet_poly_auto_chk.setChecked(self.Flags['ZEE_AUTOPOLY'])
        self.ui.inti_dir_lbl.setText(self.working_dir)
        self.ui.dop_groupbox.setChecked(self.Flags['Contonly'])
        self.ui.dop_color_force_chk.setChecked(self.Flags['dop_color_force'])
        self.ui.free_trans_helium_chk.setChecked(self.Flags["FREE_TRANS"])
        
        
        self.ui.db_observer_text.setText(self.bass_entete['observer'])
        self.ui.db_instru_text.setText(self.bass_entete['instru'])
        self.ui.db_lat_text.setText(str(self.bass_entete['lat']))
        self.ui.db_long_text.setText(str(self.bass_entete['long']))
        self.ui.db_wave_text.setCurrentText(self.bass_entete['wave_label'])
        self.ui.db_binning_combo.setCurrentIndex(self.bass_entete['binID'])
        
        self.ui.free_a_text.setText("{:5f}".format(self.my_dictini['poly_free_a']))
        self.ui.free_b_text.setText("{:5f}".format(self.my_dictini['poly_free_b']))
        self.ui.free_c_text.setText("{:5f}".format(self.my_dictini['poly_free_c']))
        self.ui.magnet_a_text.setText("{:5f}".format(self.my_dictini['poly_free_a']))
        self.ui.magnet_b_text.setText("{:5f}".format(self.my_dictini['poly_free_b']))
        self.ui.magnet_c_text.setText("{:5f}".format(self.my_dictini['poly_free_c']))
        self.ui.magnet_shift_text.setText(str(self.my_dictini['zeeman_shift']))
        self.ui.free_shift_text.setText(str(self.my_dictini['free_shift']))
        
        self.ui.dop_shift_text.setText(str(self.my_dictini['dec doppler']))
        self.ui.dop_conti_shift_text.setText(str(self.my_dictini['dec cont']))
        self.ui.free_shift1_text.setText(str(self.my_dictini['pos_free_blue']))
        self.ui.free_shift2_text.setText(str(self.my_dictini['pos_free_red']))
        #self.my_dictini['wavelength'] = self.list_wave[1][self.ui.db_wave_text.currentIndex()]
        
        self.ui.lang_button.setText(LG_str)

        self.get_Flags_tab()
        
    def save_ini(self) :
        Flags= self.Flags
        #data_entete=self.data_entete
        # recupère les données de l'UI !!! shift 1 et 2 changent si free mode
        self.Shift.append(float(conv(self.ui.inti_shift_text.text())))
        self.Shift.append(float((self.ui.dop_shift_text.text())))
        self.Shift.append(float(self.ui.dop_conti_shift_text.text()))
        self.Shift.append(float(conv(self.ui.free_shift_text.text())))
        self.Shift.append(float(conv(self.ui.magnet_shift_text.text())))
        #self.polynome =[0,0,0] # !!! 
        
        app_tab= self.ui.tab_main.currentIndex()
        # on se base sur l'index car les noms sont traduits
        # 0-General, 1-Doppler, 2-Volume, 3-Magnet, 4-Free
        if app_tab == 4 :
            self.Flags['WEAK'] =True
        else:
            self.Flags['WEAK'] =False
        if app_tab == 3 :
            self.Flags['POL'] =True
        else :
            self.Flags['POL'] =False
            
        self.Flags['NOISEREDUC']=self.ui.section_reduc_chk.isChecked()
        self.Flags['Autocrop']=self.ui.section_autocrop_chk.isChecked()
        self.Flags['zoom']=self.ui.section_zoom_chk.isChecked()
        self.Flags["FORCE_FREE_MAGN"]=self.ui.free_tilt_chk.isChecked() or self.ui.magnet_tilt_chk.isChecked() 
        self.Flags['Grid']=self.ui.section_grid_chk.isChecked()
        self.Flags['FLIPRA']=self.ui.ori_inv_EW_chk.isChecked()
        self.Flags['FLIPNS']=self.ui.ori_inv_NS_chk.isChecked()
        self.Flags['FREE_AUTOPOLY']=self.ui.free_poly_auto_chk.isChecked()
        self.Flags['ZEE_AUTOPOLY']=self.ui.magnet_poly_auto_chk.isChecked()
        self.Flags['h20percent'] = self.ui.inti_h20_radio.isChecked()
        self.Flags['Contonly'] = self.ui.dop_groupbox.isChecked()
        self.Flags['dop_color_force'] = self.ui.dop_color_force_chk.isChecked()
        self.Flags["FREE_TRANS"] = self.ui.free_trans_helium_chk.isChecked()
        
        
        self.my_dictini['lang']=self.langue
        #print(self.langue)
        self.my_dictini['directory']=self.working_dir
        if self.Flags['Grid']==True :
            self.my_dictini['grid disk']='on'
        else:
            self.my_dictini['grid disk']='off'
        """
        self.my_dictini['observer'] = self.ui.db_observer_text.text()
        self.my_dictini['instru'] = self.ui.db_instru_text.text()
        self.my_dictini['site_lat'] = float(self.ui.db_lat_text.text())
        self.my_dictini['site_long'] = float(self.ui.db_long_text.text())
        self.my_dictini['contact'] = self.ui.db_contact_text.text()
        self.my_dictini['wave_label'] = self.ui.db_wave_text.currentText()
        self.my_dictini['wavelength'] = self.list_wave[1][self.ui.db_wave_text.currentIndex()]
        """
        self.my_dictini['inversion EW']=Flags['FLIPRA']
        self.my_dictini['inversion NS']=Flags['FLIPNS']
        self.my_dictini['size_pix_cam']=self.size_pix_cam
        self.my_dictini['bin_cam']=self.bin_cam
        self.my_dictini['autocrop']=self.Flags['Autocrop']
        self.my_dictini['ext20pc']=self.Flags['h20percent']
        self.my_dictini['zoom']=self.Flags['zoom']
        #self.my_dictini['grid disk']=self.Flags['Grid']
        self.my_dictini['pos fente min']=self.param[0]
        self.my_dictini['pos fente max']=self.param[1]
        self.my_dictini['trame couronne1']=self.param[4]
        self.my_dictini['trame couronne2']=self.param[5]
        self.my_dictini['force_free_magn'] = Flags["FORCE_FREE_MAGN"]
        self.my_dictini['zeeman_shift'] = self.Shift[4]
        self.my_dictini['reduction_bruit'] = Flags['NOISEREDUC']
        self.my_dictini['angP auto']= self.ui.section_angP_auto_chk.isChecked()
        
        
        #if Flags['WEAK']:
        self.my_dictini['free_autopoly']=Flags['FREE_AUTOPOLY']
        self.my_dictini['zee_autopoly']=Flags['ZEE_AUTOPOLY']
        self.my_dictini['pos_free_blue']= float(self.ui.free_shift1_text.text())   # round(0+self.Shift[1])
        self.my_dictini['pos_free_red']= float(self.ui.free_shift2_text.text())    # round(0+self.Shift[2])
        self.my_dictini['free_shift'] = float(self.ui.free_shift_text.text()) #self.Shift[3]
       
        self.my_dictini['dec doppler']= float(self.ui.dop_shift_text.text())    # self.Shift[1]
        self.my_dictini['dec cont']= float(self.ui.dop_conti_shift_text.text()) # self.Shift[2]
        self.my_dictini['cont_only']=self.Flags['Contonly']
        self.my_dictini['dop_color_force']= self.Flags['dop_color_force']

        self.my_dictini['ang_tilt']=self.ang_tilt
        self.my_dictini['angle P']=0 # ne conserve pas angle P pour eviter confusion
        self.my_dictini['ratio_sysx']=self.ratio_fixe

        self.my_dictini['correction He']= self.Flags['FREE_TRANS']
        
        if Flags['WEAK']==False and Flags['POL']==False:
            self.my_dictini['poly_slit_a']=float(self.polynome[0])
            self.my_dictini['poly_slit_b']=float(self.polynome[1])
            self.my_dictini['poly_slit_c']=float(self.polynome[2])
        else:
            self.my_dictini['poly_free_a']=float(self.polynome[0])
            self.my_dictini['poly_free_b']=float(self.polynome[1])
            self.my_dictini['poly_free_c']=float(self.polynome[2])
            self.my_dictini['poly_slit_a']=0
            self.my_dictini['poly_slit_b']=0
            self.my_dictini['poly_slit_c']=0
            
        if Flags["SAVEPOLY"] != 0 :
            #on met ajour les params poly dans mode free et magn
            self.my_dictini['poly_free_a']=float(self.polynome[0])
            self.my_dictini['poly_free_b']=float(self.polynome[1])
            self.my_dictini['poly_free_c']=float(self.polynome[2])
            
        if Flags['WEAK']==True and Flags['FREE_AUTOPOLY']==1 :
            self.my_dictini['poly_free_a']=float(self.polynome[0])
            self.my_dictini['poly_free_b']=float(self.polynome[1])
            self.my_dictini['poly_free_c']=float(self.polynome[2])
            
        if Flags['POL']==True and Flags['FREE_AUTOPOLY']==1 :
            self.my_dictini['poly_free_a']=float(self.polynome[0])
            self.my_dictini['poly_free_b']=float(self.polynome[1])
            self.my_dictini['poly_free_c']=float(self.polynome[2])
        
            
        # sauvegarde
        try:
            with open(self.my_ini, "w") as f1:
                yaml.dump(self.my_dictini, f1, sort_keys=False)
                #print(self.my_ini)
        except :
            print (self.tr('Erreur à la sauvegarde de inti_ini.yaml : '), self.my_ini)

    def add_text(self, text) :
        self.ui.log_edit.append(text)
        self.ui.log_edit.moveCursor(QtGui.QTextCursor.End)
    
    
    
               
    def inti_open (self) :
        self.img_list =[]
        self.serfiles_short=[]
        self.serfiles=[]
        # mémorise les valeurs de crop - si 0 alors autocrop
        #self.param[2]=0
        #self.param[3]=0
        #self.ui.section_autocrop_chk.setStyleSheet('color: ivory;')
        self.ui.inti_filenames_lbl.setText('')
        posx = 5
        posy = 10
        dx = 100
        dy = 30
        serfiles = QFileDialog.getOpenFileNames(self, "Selectionner fichier(s) SER", self.working_dir, "Fichiers SER (*.ser);;Tous les fichiers png (*.png);;Fichiers disk png (*_disk.png);;Fichiers protus png (*_protus.png);;Fichiers clahe png (*_clahe.png);;Fichiers fits (*.fits)",self.pattern)
        if serfiles[1] != '' :
            if serfiles[0] :
                self.pattern = serfiles[1]
                self.serfiles=serfiles[0]
                self.working_dir= self.get_dirpath(self.serfiles[0])
                self.ui.inti_dir_lbl.setText(self.working_dir)
                self.ext= self.get_extension(self.serfiles[0])
                if self.ext != "ser" :
                    for f in self.serfiles :
                        #creer la liste de noms de fichiers short
                        self.serfiles_short.append(self.short_name(f))
                        # recupere extension
                        ext = self.get_extension (f)
                        # lecture est fonction de l'extension
                        if ext == 'png' :
                            # si png
                            img_proc=cv2.imread(f,cv2.IMREAD_UNCHANGED)
                            if len(img_proc.shape) == 3 :
                                img_proc=cv2.cvtColor(img_proc,cv2.COLOR_BGR2RGB)
                        if ext== 'fits' :
                            #si fits
                            img_proc, header=self.read_fits_image(f)
                            img_proc=np.array(img_proc,dtype='uint16')
                        
                        self.img_data = img_proc
                
                        # affiche image dans une fenetre flottante
                        myimg=img_wnd( None, self.Flags['zoom'], self)
                        
                        self.img_list.append(myimg)
                        self.img_list[-1].show()
                        self.img_list[-1].set_img(img_proc)
                        self.img_list[-1].set_title(self.short_name(f))
                        self.img_list[-1].on_ferme.connect(self.img_allclose)
                        self.img_list[-1].set_file_name(f)
                        #
                        posx=posx+dx
                        posy=posy+dy
                        self.img_list[-1].set_pos(posx, posy)
                        self.flag_ser =False
                else :
                    #print(self.working_dir)
                    self.flag_ser = True 
                    for f in self.serfiles :
                        #creer la liste de noms de fichiers short
                        self.serfiles_short.append(self.short_name(f))
                        #print(self.short_name(f))
                        
                    last_file = self.serfiles[-1]
                    fits_dateobs, nbframe, scan_size=get_data_ser(last_file)
                    
                    print(self.short_name(last_file))
                    print("H : "+str(scan_size[1])+" Nb frames : "+str(nbframe))
                    self.ui.ori_date_text.setText(fits_dateobs)
                    #print(fits_dateobs)
                    angP, paramB0, longL0, RotCarr = angle_P_B0(fits_dateobs)
                    self.ui.ori_angP_text.setText(angP)
                    self.solar_dict['B0']=paramB0
                    self.solar_dict['L0']=longL0
                    self.solar_dict['Carr']=RotCarr 

            
        else :
            print(self.tr('Aucun fichier sélectionné'))
        
        self.ui.inti_filenames_lbl.setText(','.join(self.serfiles_short))

    def inti_open_last (self) :
        #self.img_list =[]
        #self.serfiles_short=[]
        #self.serfiles=[]
        #self.ext='ser'
        files = [x for x in os.listdir(self.working_dir) if x.endswith(".ser")]
        paths = [os.path.join(self.working_dir, basename) for basename in files]
        newest = max(paths , key = os.path.getctime)
        self.inti_go_one(newest, True)
    
    def inti_go_one (self, onefile, flag_go) :
        self.img_list =[]
        self.serfiles_short=[]
        self.serfiles=[]
        self.ext='ser'
        self.ui.inti_filenames_lbl.setText(self.short_name(onefile))
        self.serfiles=[onefile]
        # prepare
        self.flag_ser = True 
        for f in self.serfiles :
            #creer la liste de noms de fichiers short
            self.serfiles_short.append(self.short_name(f))
            #print(self.short_name(f))
            
        last_file = self.serfiles[-1]
        fits_dateobs,FrameCount, scan_size=get_data_ser(last_file)
        self.ui.ori_date_text.setText(fits_dateobs)
        #print(fits_dateobs)
        angP, paramB0, longL0, RotCarr = angle_P_B0(fits_dateobs)
        self.ui.ori_angP_text.setText(angP)
        self.solar_dict['B0']=paramB0
        self.solar_dict['L0']=longL0
        self.solar_dict['Carr']=RotCarr 
        
        if flag_go :
            # go !
            print("Go")
            self.inti_go()
         
    """
    def stop_thread(self):
        if self.worker:
            print("Arrêt demandé...")
            self.worker.stop()

    def keyPressEvent(self, event):
        #Intercepte la touche ESC pour arrêter le thread
        if event.key() == Qt.Key_Escape:
            self.stop_thread()
        else:
            super().keyPressEvent(event)
    """
    
    def inti_go(self) :
        
        self.t0=time.time()
        if len(self.serfiles) == 0 :
            print(self.tr('Aucun fichier sélectionné'))
        else :
            if self.ext != 'ser' :
                print(self.tr('Fichiers png'))
            else :
                #
                # la routine de reconstruction
                #-----------------------------
                QApplication.setOverrideCursor(Qt.WaitCursor)
                print('____________________________')
                # init param
                self.ii=1
                self.geom=[0,1] # !! à verifier
                self.magnet_racines=[]
                self.image4=[]
               
                # recupere tous les parametres pour lancer inti_recon sur chaque fichier
                self.collect_data()
                
                # conditions pour ne pas appliquer la rotation d'angle P 
                if  (self.Flags["WEAK"] and self.Flags["FREE_TRANS"]) or self.Flags['POL'] or not(self.ui.section_angP_auto_chk.isChecked()):
                    self.ang_P = 0
                    self.ui.section_angP_auto_chk.setChecked(False)
                
                # demarre la boucle avec le lancement du premier fichier
                # recupère les resultats sur un signal emis par le thread
                # si plusieurs fichiers, gère à la fin du traitement des résultats
                # évite thread en parallèle
            

                if self.flag_ser == True :
                    # on commence par le premier fichier
                    serfile = self.serfiles[0]
                    self.serfile=serfile
                    WorkDir=os.path.dirname(serfile)+"/"
                    os.chdir(WorkDir)
                    # Creation des trois sous-repertoires 
                    subrep=WorkDir+'BASS2000'
                    if not os.path.isdir(subrep):
                       os.makedirs(subrep)
                    subrep=WorkDir+'Clahe'
                    if not os.path.isdir(subrep):
                       os.makedirs(subrep)
                    subrep=WorkDir+'Complements'
                    if not os.path.isdir(subrep):
                       os.makedirs(subrep)
                    
                    # extrait nom racine
                    base=os.path.basename(serfile)
                    if base=='':
                        print(self.tr('Erreur de nom de fichier : ')+serfile)
                        QApplication.restoreOverrideCursor()
                        return
                        
                    # basefich: nom du fichier ser sans extension et sans repertoire
                    self.basefich='_'+os.path.splitext(base)[0]
                    self.basefich_comple="Complements"+os.path.sep+self.basefich
                    self.basefich_bass="BASS2000"+os.path.sep+self.basefich
                    self.basefich_clahe="Clahe"+os.path.sep+self.basefich
                    #print(self.basefich)
        
                    # ajout d'un ordre de sequence pour mode magnetogramme
                    if self.Flags['POL'] and self.ii>1:
                        self.ratio_fixe=self.geom[0]
                        self.ang_tilt=self.geom[1]
                    
                    
                    debug_param = False
                    
                    
                    
                    if debug_param  :
                        # Shift
                        print(self.Shift)
                        # Flags
                        print(self.Flags)
                        # ratio_fixe
                        #print(self.ratio_fixe)
                        # ang_tilt
                        #print(self.ang_tilt)
                        # data_entete
                        #print(self.data_entete)
                        # ang_P
                        #print(self.ang_P)
                        # poly
                        #print(self.polynome)
                        # solar_dict
                        #print(self.solar_dict)
                    
                    
                    # ------------------------------------------------------------------------------------                     
                    # appel au module d'extraction, reconstruction et correction
                    # ------------------------------------------------------------------------------------     
                    debug_nothread = False  # mode debug, non thread pour accepter les points d'arret
                    
                    # change le fond du button pour indiquer qu'on est en traitement
                    self.ui.inti_go_btn.setStyleSheet ('background-color: rgb(40,0,0);')
                    
                    if debug_nothread :
                        image_queue = None
                        frames, header, cercle, range_dec, geom, self.polynome=sol.solex_proc(serfile,self.Shift,self.Flags,self.ratio_fixe,self.ang_tilt, 
                                                                                     self.polynome, self.data_entete,self.ang_P, self.solar_dict, self.param,image_queue)
                        result=[frames, header, cercle, range_dec, geom, self.polynome]
                        self.inti_result(result)
                    
                    else :
                        self.worker= Worker()
                        self.worker.result_ready.connect(self.inti_result)
                        image_queue = queue.Queue()
                        self.ongoing_thread = threading.Thread(target=self.worker.run, args=(serfile,self.Shift,self.Flags,self.ratio_fixe,self.ang_tilt, 
                                                                                     self.polynome, self.data_entete,self.ang_P, self.solar_dict, self.param,image_queue,),daemon=True)
                        
                        self.ongoing_thread.start()
                        #self.ongoing_thread.join()
    
                        # gere affichage dans une queue 
                        # affiche images
                        cv2.namedWindow('Ser', cv2.WINDOW_NORMAL)
                        cv2.moveWindow('Ser', 10, 0)
    
                        if self.Flags["RTDISP"] :
                            cv2.namedWindow('image', cv2.WINDOW_NORMAL)  
                            cv2.moveWindow('image', 0, 0)
                            cv2.namedWindow('disk', cv2.WINDOW_NORMAL)
                            cv2.moveWindow('disk', 100, 0)
    
            
                        img_mean = None
                        img_disk = None
                        img_image = None
                        frmax =  None
    
                        while True :
                            try :
                                kind,img,frmax = image_queue.get(timeout=0.1)
                                #print(kind)
                                if kind == 'mean' :
                                    img_mean = img
                                    img_disk = None
                                    img_image = None
                                    
                                if kind == 'stop' :
                                    #print("process terminé")
                                    break
                                
                                if kind == 'disk' :
                                    img_disk = img
                                    img_mean = None
                                    img_image = None
                                if kind == 'image' :
                                    img_mean = None
                                    img_disk = None
                                    img_image = img
                            except queue.Empty :
                                pass
    
                            if img_mean is not None : 
                                #print('ser')
                                iw,ih = img_mean.shape
                                cv2.imshow ('Ser', img_mean)
                                cv2.resizeWindow('Ser', int(ih), int(iw))
    
                            if img_disk is not None :
                                #print('disk '+str(frmax))
                                cv2.imshow ('disk', img_disk)
                                cv2.resizeWindow('disk', int(frmax/self.dpr), int(iw/self.dpr))
                                
                            if img_image is not None :
                                #print(('image'))
                                cv2.imshow('image',img_image)
                                cv2.resizeWindow('image', int(ih/self.dpr), int(iw/self.dpr))
                                
                                
                            
                            if cv2.waitKey(1) == 27:  # exit if Escape is hit otherwise wait 1 secondes
                            #if cv2.waitKey(1) & 0xFF == ord('q'):  # exit if Escape is hit otherwise wait 1 secondes
                                    break
    
                        #print("destroyall")
                        cv2.destroyAllWindows()

    
    def inti_result(self, result):
        
       
        #frames, header, cercle, range_dec, geom, self.polynome
        frames =result[0]
        header = result[1]
        cercle=result[2]
        range_dec=result[3]
        self.geom = result[4]
        self.polynome=result[5]
        #print(self.polynome)
            
        print(' ')
        self.img_allclose()
        self.display_images(self.serfile,frames,header,cercle,range_dec,self.ii)
        
        t2=time.time()
        dt=t2-self.t0
        #print('Traitement : '+"{:.2f}".format(dt)+' s')
        
        # on met à jour l'interface
        # format 5e si on veut en scientifique
        if self.Flags['WEAK'] or self.Flags['SAVEPOLY']:
            self.ui.free_a_text.setText("{:5f}".format(self.polynome[0]))
            self.ui.free_b_text.setText("{:5f}".format(self.polynome[1]))
            self.ui.free_c_text.setText("{:5f}".format(self.polynome[2]))
            
        if self.Flags['POL'] or self.Flags['SAVEPOLY']:
            self.ui.magnet_a_text.setText("{:5f}".format(self.polynome[0]))
            self.ui.magnet_b_text.setText("{:5f}".format(self.polynome[1]))
            self.ui.magnet_c_text.setText("{:5f}".format(self.polynome[2]))
        
        # Sauvegarde
        print(self.serfile)
        print(' ')
        
        self.inti_loop()
        
        
    def inti_loop (self) :
        
        # on check si on doit boucler
        index_ser = self.ii
        if   index_ser <= len(self.serfiles)-1:
            
            self.ii = self.ii + 1
            serfile = self.serfiles[index_ser]
            self.serfile = serfile
            #print('Loop : '+self.serfile)
            # extrait nom racine
            base=os.path.basename(serfile)
            # basefich: nom du fichier ser sans extension et sans repertoire
            #basefich=os.path.splitext(base)[0]
            if base=='':
                print(self.tr('Erreur de nom de fichier : ')+serfile)
                QApplication.restoreOverrideCursor()
                return
            
            print('_______________________________')
            self.basefich='_'+os.path.splitext(base)[0]
            self.basefich_comple="Complements"+os.path.sep+self.basefich
            self.basefich_bass="BASS2000"+os.path.sep+self.basefich
            self.basefich_clahe="Clahe"+os.path.sep+self.basefich
            #print(self.basefich)
            #
            # ------------------------------------------------------------------------------------                     
            # appel au module d'extraction, reconstruction et correction
            # ------------------------------------------------------------------------------------     
            debug_nothread = False
            
            if debug_nothread :
                image_queue = None
                frames, header, cercle, range_dec, geom, self.polynome=sol.solex_proc(serfile,self.Shift,self.Flags,self.ratio_fixe,self.ang_tilt, 
                                                                             self.polynome, self.data_entete,self.ang_P, self.solar_dict, self.param, image_queue)
                result=[frames, header, cercle, range_dec, geom, self.polynome]
                self.inti_result(result)
            
            else :
                new_title = self.img_list[-1].ui.windowTitle()+"  ...next : "+self.short_name(serfile) + ' (' +str(self.ii-1)+'/'+str(len(self.serfiles))+')'
                self.img_list[-1].ui.setWindowTitle(new_title)
                self.worker= Worker()
                self.worker.result_ready.connect(self.inti_result)
                image_queue = queue.Queue()
                self.ongoing_thread = threading.Thread(target=self.worker.run, args=(serfile,self.Shift,self.Flags,self.ratio_fixe,self.ang_tilt, 
                                                                             self.polynome, self.data_entete,self.ang_P, self.solar_dict, self.param,image_queue),daemon=True)
                self.ongoing_thread.start()
                
        
        else :
            # change le fond du button pour indiquer qu'on a fini le traitement
            self.ui.inti_go_btn.setStyleSheet ('background-color: rgb(100,0,0);')
            
            QApplication.restoreOverrideCursor()
            
            flag_galerie = False
            # affiche une galerie à la fin avec une visu image de clahe ou des images helium, conti...
            if   index_ser == len(self.serfiles) and index_ser !=1 and flag_galerie:
                self.mygalerie=galerie_wnd(self.image4, self.serfiles_short)
                self.mygalerie.ui.show()                                                    
        

    
    def display_images(self,serfile,frames,header,cercle,range_dec,ii) :
        # frames 0:raw, 1: recon aka disk
        suff = ['raw','disk','protus', 'clahe', 'cont', 'doppler']
        disks=[]
        
        
        # disk raw
        #-------------------------------------------------------------------
        if range_dec[0]==0 :
            if self.Shift[0] ==0 :
                ImgFile=self.basefich_comple+'_raw.fits'   
            else :
                dp_str= str(self.Shift[0])
                dp_str="_dp"+str.replace(dp_str,".","_")
                ImgFile=self.basefich_comple+dp_str+'_raw.fits'
        else:
            dp_str= str(self.Shift[0])
            dp_str="_dp"+str.replace(dp_str,".","_")
            ImgFile=self.basefich_comple+dp_str+'_raw.fits'
        
        hdulist = fits.open(ImgFile, memmap=False)
        hdu=hdulist[0]
        myspectrum=hdu.data
        hdulist.close()
        rih=hdu.header['NAXIS2']
        riw=hdu.header['NAXIS1']
        Disk=np.reshape(myspectrum, (rih,riw))
        # facteur de brightness sur image raw uniquement
        Ratio_lum=(65535/np.max(Disk))
        disk_raw=np.array((np.copy(Disk)*Ratio_lum),dtype='uint16')
        disk_raw=cv2.flip(disk_raw,0)
        disks.append(disk_raw)
        #sauvegarde en png de l'image raw
        if range_dec[0]==0:
            img_suffix=""
            cv2.imwrite(self.basefich_comple+img_suffix+'_raw.png',disk_raw)
        else:
            img_suffix="_dp"+str(range_dec[0])
        
        # disk recon
        #---------------------------------------------------------------------
        # Lecture  image recon
        if range_dec[0]==0:
            if self.Shift[0] ==0 :
                ImgFile=self.basefich+'_recon.fits'   
            else :
                v_str=str.replace(str(self.Shift[0]),".","_")
                ImgFile=self.basefich+'_dp'+v_str+'_cont.fits'
                suff[1]='_dp'+v_str+'_cont'
        else :
            dp_str="_dp"+str.replace(str(range_dec[0]),".","_")
            ImgFile=self.basefich+'_dp'+dp_str+'_recon.fits'
           
            
        hdulist = fits.open(ImgFile, memmap=False)
        hdu=hdulist[0]
        self.base_filename=hdu.header['FILENAME'].split('.')[0]
        
        
        # si disque partiel alors on annule rotation angle P auto sur l'interface
        ih, iw = frames[0].shape
        if cercle[2]*2 > ih :
            self.ui.section_angP_auto_chk.setChecked(False)
           
        
        frame1=np.copy(frames[0])
        frame_contrasted = seuil_image_dyn (frame1, 99.999, 0.05)
        frame_contrasted = np.array(frame_contrasted, dtype='uint16')
        frame_contrasted=cv2.flip(frame_contrasted,0)
        # ajout calcul intensité moyenne sur ROI centrée
        lum_roi= get_lum_moyenne(frame1)
        print()
        print(self.tr("Intensité moyenne : ")+ str(int(lum_roi)))
        disks.append(frame_contrasted)
        #sauvegarde en png disk
        cv2.imwrite(self.basefich+img_suffix+'_'+suff[1]+'.png',frame_contrasted)
        
        
        # disk protus # modif pour Qt après test 585
        #---------------------------------------------------------------------
        frame2=np.copy(frames[0])
        disk_limit_percent=0.0015 # black disk radius inferior by 3.5% to disk edge (was 2%) -june25
        if cercle[0]!=0:
            x0=cercle[0]
            y0=cercle[1]
            #wi=round(cercle[2])
            #he=round(cercle[3])
            wi=int(cercle[2])
            he=int(cercle[3])
            r=(min(wi,he))
            r=int(r- round(r*disk_limit_percent))-1 # retrait de 1 pixel modif de juin 2023
            # test seuils protus
            #fc2=seuil_image_dyn (frame2, 99.999, 0.05)
            flou=r*disk_limit_percent
            fc3=disk_gauss (frame2, x0, y0, int(r-(2*flou)), flou).astype(np.uint16)

            fc3=seuil_image_dyn (fc3, 99.99, 0) # 0 was 99.999 et 0.05
            #Seuil bas à zéro... fait remonter les défauts, mais conserve la dynamique
           
            frame_contrasted3 = np.array(fc3, dtype='uint16')
            frame_contrasted3=cv2.flip(fc3,0)
            disks.append(frame_contrasted3)
            #sauvegarde en png seuils protus
            cv2.imwrite(self.basefich+img_suffix+'_protus.png',frame_contrasted3)
            
            
        # disk clahe
        #---------------------------------------------------------------------
        #clahe = cv2.createCLAHE(clipLimit=0.8, tileGridSize=(2,2))
        clahe = cv2.createCLAHE(clipLimit=1.5, tileGridSize=(2,2))
        cl1 = clahe.apply(frames[0])
        #cc = seuil_image_percent (cl1, 99.9999, 25, 1.05)
        cc = seuil_image_percent (cl1, 99.9999, 25, 1)
        cc=cv2.flip(cc,0)
        disks.append(cc)
        #sauvegarde en png de clahe
        cv2.imwrite(self.basefich_clahe+img_suffix+'_clahe.png',cc)
        
        # images doppler, continuum ou free, ou magnet, ou volume
        if len(frames) >1:
            if len(frames)>3 or self.Flags['DOPCONT']: 
                # seuils image continuum 
                frameC=np.copy(frames[len(frames)-1])
                frame_continuum = seuil_image_dyn (frameC, 99.999, 0.25)
                # display image
                frame_continuum=cv2.flip(frame_continuum,0)
                disks.append(frame_continuum)
                
                # sauvegarde en png de continuum
                v_str = str.replace(str(range_dec[len(range_dec)-1]),".","_")
                cv2.imwrite(self.basefich+'_dp'+v_str+'_cont.png',frame_continuum)
            
            if self.Flags["DOPCONT"] and self.ui.dop_groupbox.isChecked():
                #on se tente une image couleur, en 8 bits
                flag_erreur_doppler=False
                try :
                    img_doppler=np.zeros([frames[1].shape[0], frames[1].shape[1], 3],dtype='uint8')
                    img_doppler_ec=np.zeros([frames[1].shape[0], frames[1].shape[1], 3],dtype='uint8')
                    if self.Flags["VOL"] :
                        f1=np.array(frames[1], dtype="float64")
                        fn=np.array(frames[len(frames)-1], dtype="float64")
                        moy=np.array(((f1+fn)/2), dtype='float64')
                    else:
                        f1=np.array(frames[1], dtype="float64")
                        f2=np.array(frames[2], dtype="float64")
                        moy=np.array(((f1+f2)/2), dtype='float64') 
                        # on equilibre les plans
                        lum_roi1 = get_lum_moyenne(f1) # calcul moyenne au centre sur zone de 200x200
                        lum_roi2 = get_lum_moyenne(f2)
                        lum_roimoy = get_lum_moyenne(moy)
                        ratio_l1 = lum_roimoy/lum_roi1
                        ratio_l2 = lum_roimoy/lum_roi2
                        frames[2] = np.clip(f2*ratio_l2, 0, 65535).astype( np.uint16)
                        frames[1] = np.clip(f1*ratio_l1, 0, 65535).astype( np.uint16) 
                        moy=np.array(moy, dtype='uint16')
                            
                        #frames[1]=np.array((moy-f1), dtype='uint16')
                        #frames[2]=np.array((f2-moy), dtype='uint16')
                        
                 
                    #i2,Seuil_haut, Seuil_bas=seuil_image(moy) # seuil bas = 0
                    Seuil_haut = np.percentile(moy, 99.99) # was 99.999
                    Seuil_bas=np.percentile(moy,5) # was zéra
                    i2= seuil_image_force (moy,Seuil_haut, Seuil_bas)
                    i1=seuil_image_force (frames[1],Seuil_haut, Seuil_bas)
                    i3=seuil_image_force(frames[2],Seuil_haut, Seuil_bas)
                    
                    # tableau pour sauvegarde des fits moy, dop1, dop2
                    moy_doppler = np.copy(i2)
                    dop1= np.copy(i1)
                    dop2= np.copy(i3)
                    
                    #conversion en 8 bits
                    i1=np.clip(i1/256,0,255).astype(np.uint8)
                    i2=np.clip(i2/256,0,255).astype(np.uint8)
                    i3=np.clip(i3/256,0,255).astype(np.uint8)
                    
                    
                    if self.Flags['DOPFLIP'] !=0 : 
                        img_doppler[:,:,0] = i3
                        img_doppler[:,:,1] = i2
                        img_doppler[:,:,2] = i1
                    else :
                        img_doppler[:,:,0] = i1 # blue
                        img_doppler[:,:,1] = i2 # green
                        img_doppler[:,:,2] = i3 # red
                    
                    flag_dop_color_force= self.ui.dop_color_force_chk.isChecked()
                    
                    if  flag_dop_color_force :
                        # remap orange en rouge 
                        img_doppler = img_force_color(img_doppler)
                        
                    img_doppler =cv2.flip(img_doppler,0)
                    disks.append(img_doppler)
                    
                    if not self.Flags["VOL"] :
                        # ajout d'une image protus doppler
                        moy=cv2.circle(moy, (x0,y0),r,80,-1,lineType=cv2.LINE_AA)
                        frames[1]=cv2.circle(frames[1], (x0,y0),r,80,-1,lineType=cv2.LINE_AA)
                        frames[2]=cv2.circle(frames[2], (x0,y0),r,80,-1,lineType=cv2.LINE_AA)
                        Th_Upper=np.percentile(moy,99.9999)*1 #preference for high contrast was 0.6
                        Th_low= 0
                        i2=seuil_image_force(moy, Th_Upper, Th_low)
                        i1=seuil_image_force (frames[1],Th_Upper, Th_low)
                        i3=seuil_image_force(frames[2],Th_Upper, Th_low)
                        
                        # passage en 8 bits
                        i3 = np.clip(i3 /256, 0, 255).astype(np.uint8)
                        i2 = np.clip(i2/256, 0, 255).astype(np.uint8)
                        i1 = np.clip(i1 /256, 0, 255).astype(np.uint8)

                        if self.Flags['DOPFLIP'] !=0 : 
                            img_doppler_ec[:,:,0] = i3
                            img_doppler_ec[:,:,1] = i2
                            img_doppler_ec[:,:,2] = i1
                        else :
                            img_doppler_ec[:,:,0] = i1 # blue
                            img_doppler_ec[:,:,1] = i2 # green
                            img_doppler_ec[:,:,2] = i3 # red
                            
                        flag_dop_color_force= self.ui.dop_color_force_chk.isChecked()
                        
                        if  flag_dop_color_force :
                            # remap orange en rouge 
                            img_doppler = img_force_color(img_doppler)
                            
                        img_doppler_ec=cv2.flip(img_doppler_ec,0)
                        # on remplace l'image protus noir et blanc
                        disks[2]=img_doppler_ec
               
                except:
                    print (self.tr('Erreur image doppler.'))
                    flag_erreur_doppler=True
                    pass
    
                #sauvegarde en png de doppler
                if flag_erreur_doppler==False:
                    dp_str= str(abs(range_dec[1]))
                    dp_str='_'+str.replace(dp_str,".","_")
                    img_doppler_cv2 =  cv2.cvtColor(img_doppler, cv2.COLOR_RGB2BGR)
                    img_doppler_ec_cv2 =  cv2.cvtColor(img_doppler_ec, cv2.COLOR_RGB2BGR)
                    cv2.imwrite(self.basefich+'_doppler'+dp_str+'.png',img_doppler_cv2)
                    cv2.imwrite(self.basefich+'_doppler_protus'+dp_str+'.png',img_doppler_ec_cv2)
    
            
            if self.Flags["POL"] :

                try :
                    fr1=np.copy(frames[1])
                    fr2=np.copy(frames[2])
                    fr0=np.copy(frames[0])
                    img_diff=np.array(np.array(fr2, dtype='float32')-np.array(fr1, dtype='float32'),dtype='int32')
                except:
                    print (self.tr('Erreur image difference '))
                
                # et on en fait quoi ?? !!!

            if self.Flags["WEAK"] :
                flag_erreur_weak=False
                
                try :
                    fr0 = frames[0].astype(np.float64)
                    fr1 = frames[1].astype(np.float64)
                    fr2 = frames[2].astype(np.float64)
                    
                    moy = (fr1 + fr2) * 0.5
                    img_diff = fr2 - fr1
                    
                    d = fr0 - moy
                    
                    offset = max(0.0, -d.min())
                    print("offset: "+str(offset))

                    if offset > 0:
                        offset += 100
                        if self.Flags['Couronne']:
                            offset = min(offset, 1000)
                            #offset = 0
                    else :
                        offset=0
                   
                    #print("Offset "+str(offset))
                    img_weak_array = np.clip(d + offset, 0, 65535)
                    img_weak_uint = img_weak_array.astype(np.uint16)
                    
                    if self.Flags['Couronne'] and 2==1:
                        #disk_limit_percent=0.0015 # black disk radius inferior by 3.5% to disk edge (was 2%) -june25
                        disk_limit_percent=0.003
                        if cercle[0]!=0:
                            x0=cercle[0]
                            y0=cercle[1]
                            #wi=round(cercle[2])
                            #he=round(cercle[3])
                            wi=int(cercle[2])
                            he=int(cercle[3])
                            r=(min(wi,he))
                            #r=int(r- round(r*disk_limit_percent))-1 # retrait de 1 pixel modif de juin 2023
                            r=int(r- round(r*disk_limit_percent))
                            
                            flou=r*disk_limit_percent
                            img_weak_uint=disk_gauss (img_weak_uint, x0, y0, int(r-(2*flou)), flou).astype(np.uint16)
                    
                    #Seuil_bas=int(offset//2)
                    Seuil_bas=0
                    Seuil_haut=int(np.percentile(img_weak_uint,99.99))
                    
                    
                    #print("free seuils : ", Seuil_bas, Seuil_haut)
                    if (Seuil_haut-Seuil_bas) != 0 :
                        img_weak=seuil_image_force(img_weak_uint, Seuil_haut, Seuil_bas)
                    else:
                        img_weak=np.array((img_weak_array), dtype='uint16')
                    
                    #DiskHDU=fits.PrimaryHDU(img_weak,header)
                    #DiskHDU.writeto(self.basefich+'_qt_test0.fits', overwrite='True')
                    
                    if self.Flags["FREE_TRANS"] :
                        R=int(cercle[2])
                        #print(R, img_weak.shape[0])
                        if R*2 >= img_weak.shape[0] :
                            print("Pas de correction de transversallium sur un disque partiel")
                        else :
                            result_image = corrige_trans_helium(img_weak, R)
                            
                            #DiskHDU=fits.PrimaryHDU(result_image,header)
                            #DiskHDU.writeto(self.basefich+'_qt_test1.fits', overwrite='True')
                            
                            # On ajoute le continuum méthode sunscan
                            image1 = np.array(result_image, dtype=np.int32) 
                            image2 = np.array(moy, dtype=np.int32)
    
                            constant = 32767
    
                            coef = 0.8 # 0.6 sur sunscan app
                            image2_transformed = np.where(image2 > 0, image2 - constant, 0)
                            result_image = image1 + coef * image2_transformed
                            result_image = np.clip(result_image, 0, 65535).astype(np.uint16)
                            max_value = np.max(result_image)
                            result_image = (result_image / max_value) * 65535.0
                            result_image = result_image.astype(np.uint16)
                            
                            # etape supplémentaire de merge avec fond
                            height, width = img_weak.shape
                            center = (cercle[0], cercle[1])
                            feather_width=15
                            radius=R-feather_width-1
        
                            # Create the circular mask
                            mask = create_circular_mask((height, width), center, radius, feather_width)
                            # Blend the images
                            img_weak = blend_images(img_weak, result_image, mask)
                    #
                    img_weak=np.array(img_weak, dtype='uint16')
                    img_weak=cv2.flip(img_weak,0)
                    disks.append(img_weak)
                    img_weak_array=d
                    suff[4] = 'free'
                
                except:
                    print (self.tr('Erreur image '))
                    flag_erreur_weak=True
                    pass
    
                #sauvegarde en png de raie libre 
                if flag_erreur_weak==False:
                    cv2.imwrite(self.basefich+'_free'+'.png',img_weak)
                    
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Affiche les images   
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
        self.img_list=[]
        self.myzoom_wnd=zoom_wnd()
        
        if len(disks)<=4 :
            self.image4.append(disks[3])
        else :
            self.image4.append(disks[4])
        
        if self.Flags['zoom'] :
            
            self.myzoom_wnd.show()
        posx0=50
        dx=self.myscreen_w//20
        dy=20
        
        for i in range(0, len(disks)) :
            myimg=img_wnd(self.myzoom_wnd, self.Flags['zoom'], self)
            
            self.img_list.append(myimg)
            self.img_list[-1].show()
            self.img_list[-1].set_img(disks[i])
            if suff[i] == 'doppler' :
                suff[i] = suff[i]+dp_str
            if suff[i] == 'cont' :
                if len(range_dec)==1 :
                    v_str=str.replace(str(self.Shift[0]),'.','_')
                    suff[i] = 'dp'++'_cont'
                else :
                    v_str=str.replace(str(range_dec[len(range_dec)-1]),'.','_')
                    suff[i] = 'dp'+v_str+'_cont'
            self.img_list[-1].set_title(self.racine_name(self.short_name(serfile))+'_'+suff[i])
            self.img_list[-1].set_file_name(self.working_dir+os.sep+'_'+self.racine_name(self.short_name(serfile))+'_'+suff[i]+'.png')
            self.img_list[-1].on_ferme.connect(self.img_allclose)
            #
            dposx=dx*i
            dposy=dy*i
            posx=posx0+dposx
            self.img_list[-1].set_pos(posx, dposy)
        #self.myzoom_wnd.show()


# ------------------------------------------------------------------------------------             
# génère images additionnelles non affichées
# ------------------------------------------------------------------------------------
        
        # image contraste inversé
        # ------------------------
        #frame_sub=65535-frame_contrasted
        frame = frame_contrasted - frame_contrasted.min()
        frame_sub = frame.max() - frame        
        if self.files_to_save['inv'] : cv2.imwrite(self.basefich+img_suffix+'_inv.png',frame_sub)
        
        # test mix mais c'est moche...        
        frame_sub[frame_sub==65535]=1000
        frame_combo=np.maximum(frame_contrasted3, frame_contrasted)
        frame_combo[frame_combo==65535]=0        
        if self.files_to_save['mix'] : cv2.imwrite(self.basefich+img_suffix+'_mix.png',frame_combo)
        
        
        
        couleur = self.bass_entete['wave_label']
        
        # image disk colorisée
        # --------------------
        if self.files_to_save['color'] :
            img_color = Colorise_Image(couleur, frame_contrasted, self.basefich, img_suffix)
            #img_color = Colorise_Image(couleur, cce, basefich, img_suffix)
        
        # si pas de flag mode magneto et ou raie libre
        if self.Flags["POL"] != True  and  self.Flags["WEAK"] != True : 
            
            # sauve en format jpg pour BASS2000 si pas en manuel
            if couleur != 'Manual' :
                cv2.imwrite("BASS2000"+os.path.sep+self.base_filename+'.jpg',frame_contrasted*0.0038)
                cv2.imwrite("BASS2000"+os.path.sep+self.base_filename+'_protus.jpg',frame_contrasted3*0.0038)

            # image avec stonyhurtdisk
            # ------------------------
            nomfich=self.basefich+img_suffix+'_disk.png'
            if not file_exist(nomfich) :
                nomfich=self.basefich+img_suffix+'_cont.png'
            
            nomrep1=self.working_dir+os.path.sep
            nomrep2=nomrep1+"BASS2000"+os.path.sep
            fich_param={}
            graph_param={}
            fich_param['date'] = header['DATE-OBS']
            #hdr['SEP_LAT']=float(solar_dict['B0'])
            #hdr['SEP_LON']=float(solar_dict['L0'])
            #hdr['SOLAR_P']=float(ang_P)
            #hdr['CAR_ROT']=float(solar_dict['Carr'])

            fich_param['P']=0
            fich_param['PDisp']=0
            #grid_on=True
            if self.Flags ['Grid']==1 :
                try :
                    fich_param['P'] = header['SOLAR_P']  # Only for display, INTI puts heliographic pole on top
                    fich_param['PDisp'],a,b,c=angle_P_B0 (fich_param['date'] ) # modif 14 janvier, il faut l'afficher
                except:
                    fich_param['PDisp'],a,b,c=angle_P_B0 (fich_param['date'] )
                try :
                    fich_param['B0'] = float(header['SEP_LAT'])
                    fich_param['L0'] = float(header['SEP_LON'])
                    if fich_param['B0'] == 0 :  # modif du 14 janvier 24
                        a, fich_param['B0'],fich_param['L0'],c=angle_P_B0 (fich_param['date'] )
                        fich_param['B0']=float(fich_param['B0'])
                        fich_param['L0']=float(fich_param['L0'])
                        
                except:
                    a, fich_param['B0'],fich_param['L0'],c=angle_P_B0 (fich_param['date'] )
                    fich_param['B0']=float(fich_param['B0'])
                    fich_param['L0']=float(fich_param['L0'])
                    
                fich_param['xcc'] = x0
                
                ih = frames[0].shape[0]
                if y0<ih or y0>ih :
                    fich_param['ycc'] = ih-y0
                else :
                    fich_param['ycc'] = y0
                
                fich_param['radius'] = r
                graph_param['gradu'] = self.flag_grid_disp
                graph_param['opacity'] = 0.5
                graph_param['lwidth'] = 0.2
                graph_param['color'] = self.grid_color
                graph_param['color_inv'] = 'black'
                graph_param['disp']=False
                
                # creation et sauvegarde du fichier avec grille coor helio
                sth.draw_stonyhurst(nomrep1, nomrep2,nomfich, fich_param, graph_param)
 
# ------------------------------------------------------------------------------------             
# sauve les png multispectraux et cree une animation gif
       
        if (self.Flags["VOL"] or self.Flags["POL"]) and len(range_dec)!=1:
            
            k=1
            
            #ROi de 200 pixels autour du centre du disque
    
            dim_roi=200
            rox1= cercle[0]-dim_roi
            rox2=cercle[0]+dim_roi
            roy1=cercle[1]-dim_roi
            roy2=cercle[1]+dim_roi
            
            # calcul moyenne intensité au centre du disque
            frame0=np.copy(frames[0])
            sub_frame=frame0[5:,:-5]
            Seuil_haut= np.percentile(sub_frame,99.999)
            
            lum_roi=np.mean(frame0[roy1:roy2,rox1:rox2])
            lum_roi_ref=lum_roi*(65535/Seuil_haut)
            #lum_roi_ref=lum_roi
            lum_roi_max=np.mean(frames[1][roy1:roy2,rox1:rox2])
           
            sub_frame=frame0[5:,:-5]*(lum_roi_ref/lum_roi)
            Seuil_haut=np.percentile(sub_frame,99.999)
            if Seuil_haut>65535: Seuil_haut=65535
            #Seuil_bas=(Seuil_haut*0.15)
            Seuil_bas=0
    
    
            a=1-(lum_roi_ref/lum_roi_max)*0.6       
            mlum=a/((len(range_dec)-1)/2)
            #print('mlum : ', mlum)
            
            if self.Flags["POL"] :
                # zeeman
                flag_lum_auto=False
                flag_seuil_auto=True
            else:
                # volume
                flag_lum_auto=True
                flag_seuil_auto=False
            
            # creation volume animation
            ih,iw=frames[0].shape
            self.vol_images=np.zeros((len(range_dec),iw, ih),dtype='uint16')
            
            # on recrer une liste de range_dec en mettant le dec à zéro au centre
            ra1=range_dec[1:((len(range_dec)-1)//2)+1]
            ra2=range_dec[((len(range_dec)-1)//2)+1:]
            range_dec2 = np.concatenate([ra1,[0],ra2])
            #print(range_dec2)
            mid_size=len(frames)//2
            frames[:mid_size]=frames[1:mid_size+1]
            frames[mid_size]=frame0
            
            # on boucle was 1
            for i in range (0,len(range_dec2)) :

                # image seuils moyen 
                framedec=np.copy(frames[i])
                sub_frame=framedec[5:,:-5]
                
                # calcul moyenne intensité au centre du disque
                lum_roi=np.mean(sub_frame[roy1:roy2,rox1:rox2])
                lum_coef=(lum_roi_ref/lum_roi)+(abs(range_dec2[i])*mlum)
                #print(lum_roi_ref/lum_roi)
                
                if flag_lum_auto==True : 
                    framedec=framedec*lum_coef
                
                if flag_seuil_auto==True:
                    Seuil_haut=np.percentile(framedec[5:-5],99.999)
                    if Seuil_haut>65535: Seuil_haut=65535
                    #Seuil_bas=(Seuil_haut*0.15)
                    Seuil_bas=0
              
                framedec[framedec>Seuil_haut]=Seuil_haut
                fdec=(framedec-Seuil_bas)* (65535/(Seuil_haut-Seuil_bas))
                fdec[fdec<0]=0
                frame_dec=np.array(fdec, dtype='uint16')
                #print (k, lum_coef)
                
                frame_dec=cv2.flip(frame_dec,0)
                cv2.imwrite(self.basefich+'_dp'+str(k)+'.png',frame_dec)
                num=np.flipud(np.rot90(frame_dec))
                self.vol_images[k-1]=num
                k=k+1
            
                   
            debug= False
            if debug :
                for j in range(0, len(self.vol_images)) :
                    plt.imshow(self.vol_images[j])
                    plt.title(str(j))
                    plt.show()
            
            if  self.Flags["POL"]==False:
                self.img_list[4].set_vol(self.vol_images)
                self.img_list[4].set_title('Volume')
                filename_gif=self.basefich+'.gif'
                
                self.frps = 1
                self.gif_images=[]
                # mp4
                height,width = self.vol_images[0].shape
                filename_mp4 = self.basefich+'.mp4'
                fourcc = cv2.VideoWriter_fourcc(*'avc1')
                out = cv2.VideoWriter(filename_mp4, fourcc,self.frps, (width, height),0)
                
                for i in range(1,len(range_dec)+1) :
                    filename=self.basefich+'_dp'+str(i)+'.png'
                    img=cv2.imread(filename)
                    #img2=cv2.resize(img,(width, height),interpolation = cv2.INTER_AREA)
                    self.gif_images.append(Image.fromarray(img))
                    out.write(img)
                
                out.release()
                im=self.gif_images[0]
                im.save(filename_gif, save_all=True, append_images=self.gif_images[1:], duration=int(1/self.frps)*1000, loop=0)
        
                    
 # ------------------------------------------------------------------------------------      
 # sauvegarde les fits
 
        frame2=np.copy(frames[0])
        if self.Flags['WEAK']:
            frame2=np.copy(img_weak_array)
            DiskHDU=fits.PrimaryHDU(frame2,header)
            DiskHDU.writeto(self.basefich+'_free.fits', overwrite='True')
            if len(self.serfiles)>1 :
                DiskHDU.writeto(self.basefich+'_free-'+str(ii)+'.fits', overwrite='True')
            frame3=np.copy(img_diff)
            # sauve la différence
            if self.files_to_save['diff'] : 
                DiskHDU=fits.PrimaryHDU(frame3,header)
                DiskHDU.writeto(self.basefich_comple+'_diff.fits', overwrite='True')
            # sauve les images individuelles pour debug
            if self.files_to_save['x0x1x2']:
                DiskHDU=fits.PrimaryHDU(fr0,header)
                DiskHDU.writeto(self.basefich_comple+'_x0.fits', overwrite='True')
                DiskHDU=fits.PrimaryHDU(fr1,header)
                DiskHDU.writeto(self.basefich_comple+'_x1.fits', overwrite='True')
                DiskHDU=fits.PrimaryHDU(fr2,header)
                DiskHDU.writeto(self.basefich_comple+'_x2.fits', overwrite='True')
            # sauve la somme
            if self.files_to_save['sum'] :
                DiskHDU=fits.PrimaryHDU(moy,header)
                DiskHDU.writeto(self.basefich_comple+'_sum.fits', overwrite='True')
            
        if self.Flags['DOPCONT'] and not self.Flags['Contonly'] :
            if self.files_to_save['dop']:
                DiskHDU=fits.PrimaryHDU(moy_doppler,header)
                DiskHDU.writeto(self.basefich_comple+'_moy.fits', overwrite='True')
                
                DiskHDU=fits.PrimaryHDU(dop1,header)
                DiskHDU.writeto(self.basefich_comple+'_dop1.fits', overwrite='True')
    
                DiskHDU=fits.PrimaryHDU(dop2,header)
                DiskHDU.writeto(self.basefich_comple+'_dop2.fits', overwrite='True')

        
        if self.Flags['POL']:
            # image difference en fits
            if self.files_to_save['diff'] :
                frame3=np.copy(img_diff)
                DiskHDU=fits.PrimaryHDU(frame3,header)
                DiskHDU.writeto(self.basefich_comple+'_diff.fits', overwrite='True')
            # renomme les dp fits en r et b
            ImgFileb=self.basefich+'_dp'+str(range_dec[1])+'_recon.fits'
            ImgFiler=self.basefich+'_dp'+str(range_dec[2])+'_recon.fits'
            os.replace(ImgFileb,self.magnet_racines[0]+str(ii)+".fits")
            os.replace(ImgFiler,self.magnet_racines[1]+str(ii)+".fits")
            print('Sorties :')
            print(self.magnet_racines[0]+str(ii)+".fits")
            print(self.magnet_racines[1]+str(ii)+".fits")
            # renomme les dp png en r et b
            ImgFileb=self.basefich+'_dp1.png'
            ImgFiler=self.basefich+'_dp3.png'
            os.replace(ImgFileb,self.magnet_racines[0]+str(ii)+".png")
            os.replace(ImgFiler,self.magnet_racines[1]+str(ii)+".png")
    
        if self.Flags['VOL'] and self.Flags["FITS3D"]:
            
            hdr_vol=header
            hdr_vol['NAXIS']=3
            hdr_vol['NAXIS3']=len(frames)
            hdr_vol['CTYPE3']='wavelength'
            hdr_vol['CUNIT3']='Angstrom'
            # insert fullprofile dans le nom
            a=header['FILENAME'].split('SOLEX')
            try :
                filename_3D=a[0]+'SOLEX_fullprofile'+a[1]
            except:
                filename_3D=a[0]+'SOLEX_fullprofile.fits'
            hdr_vol['FILENAME']=filename_3D
            #print(hdr_vol['FILENAME'])
            hdr_vol['CDELT3']=0.155
            """
            n=len(frames)-1
            demi_n=int(n/2)
            # remet la premiere frame au centre car correspond a dec=0
            t=frames[0]
            frames[0:demi_n-1]=frames[1:demi_n]
            frames[demi_n+1:]=frames[demi_n:]
            frames[demi_n]=t
            """
            DiskHDU=fits.PrimaryHDU(frames,hdr_vol)
            DiskHDU.writeto(filename_3D, overwrite='True')

    
    def collect_data (self ):
        
        # Flags & Shifts
        self.Shift=[]
        self.get_Flags_tab()
        
        self.Flags["ALLFITS"] = True #historic      
        self.Flags["RTDISP"] = self.ui.section_realtime_chk.isChecked()
        self.Flags['NOISEREDUC']=self.ui.section_reduc_chk.isChecked()
        self.Flags['Autocrop']=self.ui.section_autocrop_chk.isChecked()
        self.Flags["SAVEPOLY"]=self.ui.section_savepoly_chk.isChecked()
        self.Flags['Grid']=self.ui.section_grid_chk.isChecked()
        self.Flags["FORCE"] = self.ui.section_tilt_chk.isChecked()
        self.Flags['FLIPRA']=self.ui.ori_inv_EW_chk.isChecked()
        self.Flags['FLIPNS']=self.ui.ori_inv_NS_chk.isChecked()
        self.Flags["DOPFLIP"] = self.ui.dop_inv_chk.isChecked()
        self.Flags["FITS3D"] = self.ui.seq_fits3D_chk.isChecked()
        self.Flags["ZEE_AUTOPOLY"] = self.ui.magnet_poly_auto_chk.isChecked()
        self.Flags['FREE_AUTOPOLY']=self.ui.free_poly_auto_chk.isChecked()
        self.Flags["FREE_TRANS"] = self.ui.free_trans_helium_chk.isChecked()
        self.Flags['zoom'] = self.ui.section_zoom_chk.isChecked()
        self.Flags['h20percent'] = self.ui.inti_h20_radio.isChecked()
        self.Flags['Couronne'] = self.ui.free_corona_chk.isChecked()
        self.Flags['Contonly']= not self.ui.dop_groupbox.isChecked()
        
        
        # ratio_fixe, ang_tilt
        self.ratio_fixe= to_float(self.ui.inti_ratio_text.text())
        self.ang_tilt = self.ui.inti_tilt_text.text()
        
        # les shifts
        if self.Flags["WEAK"] :
            poly=[]
            poly.append(float(self.ui.free_a_text.text()))
            poly.append(float(self.ui.free_b_text.text()))
            poly.append(float(self.ui.free_c_text.text()))
            self.Shift.append(0)
            self.Shift.append(float(self.ui.free_shift1_text.text())) # decalage free bleu ou decalage 1
            self.Shift.append(float(self.ui.free_shift2_text.text()))
            self.Shift.append(to_float(conv(self.ui.free_shift_text.text())))
            self.Flags["FORCE"]=self.ui.free_tilt_chk.isChecked()
            
            if self.ui.free_tilt_chk.isChecked() != True:
                self.ratio_fixe=0
                self.ang_tilt=0
            else :
                self.ratio_fixe=to_float(conv(self.ui.free_ratio_text.text()))
                self.ang_tilt=to_float(conv(self.ui.free_tilt_text.text()))
            
            self.ui.magnet_tilt_chk.setChecked( self.ui.free_tilt_chk.isChecked())
            self.polynome=poly
        
        if self.Flags["POL"]:
            poly=[]
            poly.append(float(self.ui.magnet_a_text.text()))
            poly.append(float(self.ui.magnet_b_text.text()))
            poly.append(float(self.ui.magnet_c_text.text()))
            self.Flags["FORCE"]=self.ui.magnet_tilt_chk.isChecked()
            
            if self.ui.magnet_tilt_chk.isChecked() != True: 
                # mode auto
                self.ratio_fixe=0
                self.ang_tilt=0
            else :
                # mode force
                self.ratio_fixe=to_float(conv(self.ui.magnet_ratio_text.text()))
                self.ang_tilt=to_float(conv(self.ui.magnet_tilt_text.text()))

            self.ui.free_tilt_chk.setChecked( self.ui.magnet_tilt_chk.isChecked())
            self.polynome=poly
        
        if self.ui.section_tilt_chk.isChecked() != True and self.Flags["WEAK"] != True :
            if  self.Flags["POL"] != True :
                self.ratio_fixe=0
                self.ang_tilt=0
        
        self.Flags["FORCE_FREE_MAGN"]=self.ui.magnet_tilt_chk.isChecked() or self.ui.free_tilt_chk.isChecked()
        
        if self.Flags["DOPCONT"] :
            self.Shift.append(0)   
        else :
            self.Shift.append(float(self.ui.inti_shift_text.text()))
        
        self.Shift.append(float(self.ui.dop_shift_text.text()))
        self.Shift.append(float(self.ui.dop_conti_shift_text.text()))

        
        if self.Flags["POL"]:
            self.Shift.append(int(self.ui.magnet_wide_text.text()))
            self.Shift.append(float(self.ui.magnet_shift_text.text()))
        else:
            self.Shift.append(int(self.ui.seq_range_text.text()))
            self.Shift.append(0.0)
        
        self.Shift.append(float(self.ui.dop_offset_text.text()))
        #print('shift dop ', Shift[5])
        
        self.magnet_racines.append(self.ui.magnet_racineB_text.text())
        self.magnet_racines.append(self.ui.magnet_racineR_text.text())
        
        # data_entete
        """
        self.my_dictini['observer'] = self.ui.db_observer_text.text()
        self.my_dictini['instru'] = self.ui.db_instru_text.text()
        self.my_dictini['site_lat'] = self.ui.db_lat_text.text()
        self.my_dictini['site_long'] = self.ui.db_long_text.text()
        self.my_dictini['contact'] = self.ui.db_contact_text.text()
        self.my_dictini['wavelength'] = self.list_wave[1][self.ui.db_wave_text.currentIndex()]
        self.my_dictini['wave_label'] = self.ui.db_wave_text.currentText()
        """
                
        self.bass_entete['observer'] = self.ui.db_observer_text.text()
        self.bass_entete['instru'] = self.ui.db_instru_text.text()
        self.bass_entete['lat'] = self.ui.db_lat_text.text()
        self.bass_entete['long'] = self.ui.db_long_text.text()
        self.bass_entete['wavelength'] = self.list_wave[1][self.ui.db_wave_text.currentIndex()]
        self.bass_entete['wave_label'] = self.ui.db_wave_text.currentText()
        self.bass_entete['version']= "INTI "+self.version
        self.bass_entete['binning']= self.ui.db_binning_combo.currentText()
        self.bass_entete['binID']= self.ui.db_binning_combo.currentIndex()
        
        self.data_entete=[self.bass_entete['observer'], 
                          self.bass_entete['instru'],
                          self.bass_entete['long'],
                          self.bass_entete['lat'],
                          self.bass_entete['contact'],
                          self.bass_entete['wavelength'],
                          self.bass_entete['wave_label'],
                          self.bass_entete['camera'],
                          self.bass_entete['binning'],
                          self.bass_entete['pixel'],
                          self.bass_entete['version'],
                          self.bass_entete['focinstru'],
                          self.bass_entete['diaminstru'],
                          self.bass_entete['spectro'],
                          self.bass_entete['objcollim'],
                          self.bass_entete['objcam'],
                          self.bass_entete['reseau'],
                          self.bass_entete['ordre'],
                          self.bass_entete['angle'],
                          self.bass_entete['fentelong'],
                          self.bass_entete['fentelarg'],
                          self.bass_entete['diaph'],
                          self.bass_entete['dfilter'],
                          
                          ]
        
        # ang_P mis a jour à l'ouverture fichier avec la date
        if self.ui.ori_angP_text.text() == '' :
            self.ang_P = 0 
        else :
            try :
                self.ang_P=float(self.ui.ori_angP_text.text())
            except :
                print(self.tr('erreur format angle P'))
                self.ang_P=0
    
        # solar_dict > à l'ouverture fichier avec la date
        
        # param > uniquement dans read_ini
        
        
    
    # tab general
    #-------------------------------------------------------------------------
    
    def inti_section_show(self) :
        if self.ui.inti_section_btn.isChecked() :
            self.ui.dock_main.show()
            self.ui.panel_stack.setCurrentIndex(0)
            self.ui.inti_section_btn.setArrowType(Qt.LeftArrow)
            self.ui.inti_database_btn.setChecked(False)
            self.ui.inti_database_btn.setArrowType(Qt.RightArrow)
            self.ui.inti_orientation_btn.setChecked(False)
            self.ui.inti_orientation_btn.setArrowType(Qt.RightArrow)
        else :
            self.ui.dock_main.hide()
            self.ui.inti_section_btn.setArrowType(Qt.RightArrow)

    def inti_database_show (self) :
        if self.ui.inti_database_btn.isChecked() :
            self.ui.dock_main.show()
            self.ui.panel_stack.setCurrentIndex(1)
            self.ui.inti_database_btn.setArrowType(Qt.LeftArrow)    
            self.ui.inti_section_btn.setChecked(False)
            self.ui.inti_section_btn.setArrowType(Qt.RightArrow)
            self.ui.inti_orientation_btn.setChecked(False)
            self.ui.inti_orientation_btn.setArrowType(Qt.RightArrow)
        else :
            self.ui.dock_main.hide()
            self.ui.inti_database_btn.setArrowType(Qt.RightArrow)
            
    def inti_orientation_show (self) :
        if self.ui.inti_orientation_btn.isChecked() :
            self.ui.dock_main.show()
            self.ui.panel_stack.setCurrentIndex(2)
            self.ui.inti_orientation_btn.setArrowType(Qt.LeftArrow)      
            self.ui.inti_database_btn.setChecked(False)
            self.ui.inti_database_btn.setArrowType(Qt.RightArrow)
            self.ui.inti_section_btn.setChecked(False)
            self.ui.inti_section_btn.setArrowType(Qt.RightArrow)
        else :
            self.ui.dock_main.hide()
            self.ui.inti_orientation_btn.setArrowType(Qt.RightArrow)

    def crop_manu (self) :
        pos = QtGui.QCursor.pos()
        self.mycrop_box = crop_popup(self)
        self.mycrop_box.valeurs_onclose.connect(self.valeurs_received)
        self.mycrop_box.ui.move(pos)
        self.mycrop_box.ui.show()

    def valeurs_received (self, H_val, L_val) :
        print(f"Crop H : {H_val}, L : {L_val}")
        self.param[2] = H_val
        self.param[3] = L_val
        if H_val != '0' or L_val != '0' :
            self.ui.section_autocrop_chk.setStyleSheet('color: orange;')
            self.ui.inti_h20_radio.setChecked(False)
        else :
            self.ui.section_autocrop_chk.setStyleSheet('color: ivory;')
        
    def db_edit_entete (self):
        
        self.bass_entete['observer'] = self.ui.db_observer_text.text()
        self.bass_entete['instru'] = self.ui.db_instru_text.text()
        self.bass_entete['lat'] = self.ui.db_lat_text.text()
        self.bass_entete['long'] = self.ui.db_long_text.text()
        self.bass_entete['waveID'] =self.ui.db_wave_text.currentIndex()
        self.bass_entete['binID']=self.ui.db_binning_combo.currentIndex()
        
        self.myentete_dlg=entete_dialog(self.bass_entete)
        self.myentete_dlg.ui.finished.connect(self.get_entete)
        
    def db_print_wave (self) :
        index= self.ui.db_wave_text.currentIndex()
        print(self.list_wave[1][index])
    
    def get_entete (self) :
        self.bass_entete = self.myentete_dlg.get_parameters()
        #print("Nouveaux paramètres :", self.bass_entete)
        
        self.ui.db_observer_text.setText(self.bass_entete['observer'])
        self.ui.db_instru_text.setText(self.bass_entete['instru'])
        self.ui.db_lat_text.setText(self.bass_entete['lat'])
        self.ui.db_long_text.setText(self.bass_entete['long'])
        self.ui.db_wave_text.setCurrentIndex(self.bass_entete['waveID'])
        self.ui.db_binning_combo.setCurrentIndex(self.bass_entete['binID'])

    # tab doppler
    #-------------------------------------------------------------------------
    
    # tab free line
    #-------------------------------------------------------------------------
    
            
    def inti_calc (self) :
        self.my_calc=calc_dialog(self.shg_foc_camera_other, self.shg_angle_other)
        self.my_calc.show()
    
    
        
    #--------------------------------------------------------------------------
    # fonctions utilitaires
    #--------------------------------------------------------------------------
    
    
    
    
        
    def get_Flags_tab(self) :
        # check dernier tab
        self.Flags['VOL'] = False
        self.Flags['DOPCONT']=False
        self.Flags['WEAK']= False
        self.Flags['POL'] = False
        
        app_tab= self.ui.tab_main.tabText(self.ui.tab_main.currentIndex())
        #print('tab  :'+str(app_tab))
        
        if app_tab =='Helium-Couronne-Libre' or app_tab == 'Helium-Coronna Free':
            self.Flags['WEAK'] =True
        if app_tab == 'Magnétogramme' or app_tab=='Magnetogram':
            self.Flags['POL'] =True
        if app_tab == 'Doppler - Continuum' :
            self.Flags['DOPCONT'] =True
        if app_tab == 'Séquence Doppler' or app_tab=='Doppler Sequence':
            self.Flags['VOL'] =True

    def img_allclose (self):
        #print("all close")
        
        try :
            if self.Flags['zoom'] :
                self.myzoom_wnd.ui.close()
            for f in self.img_list :
                f.write_settings()    
                f.ui.close()                         
        except :
            pass
    
    def ori_gong(self):
        # cela n'a du sens que si on a deja une image _recon
        fits_dateobs = self.ui.ori_date_text.text()
        if fits_dateobs!='' :
            fname='_'+self.get_basename(self.short_name(self.serfiles[-1]))+'_disk.png'
            filename= self.working_dir+os.sep+fname
        
            self.mygong=gong_wnd()
            
            datemonth=fits_dateobs.split('T')[0].replace('-','')[:6]
            dateday=fits_dateobs.split('T')[0].replace('-','')
            r1="https://gong2.nso.edu/HA/hag/"+datemonth+"/"+dateday+"/"
            Vo_req=r1

            reponse_web=rq.get(Vo_req)
            sun_meudon=reponse_web.text.split("\n")
            t=sun_meudon[11].split('href=')[1].split(">")[0]
            t=t.replace('"','')
            #print(r1+t)
            url = r1+t 
            img_data = requests.get(url).content
            nparr =np.frombuffer(img_data,np.uint8)
            img_gong=cv2.imdecode(nparr,cv2.IMREAD_GRAYSCALE)
            hg,wg=img_gong.shape
            img_gong_orig=np.copy(img_gong)
            
            if os.path.exists(filename):
                #web.open(filename)
                #pixmap2 = QtGui.QPixmap(filename)
                img_disk=cv2.imread(filename,cv2.IMREAD_UNCHANGED)
                ih, iw = img_disk.shape
                if ih != iw :
                    # a priori image toujours plus large que haute...
                    pad_h= (iw-ih)//2
                    pad_zone = np.zeros((pad_h,iw), dtype= 'uint16')
                    img_disk= np.concatenate((pad_zone,img_disk,pad_zone))
                    
            
            
            # traitement images pour trouver les deux diamètres img_gong et img_disk
            diam_gong = img_getdiam (img_gong)
            diam_disk = img_getdiam (img_disk)
            
            # mettre les deux images à la même dimension centrée et avec marge de 100 pixels
            ratio =diam_disk/diam_gong
            img_gong=cv2.resize(img_gong,(int(hg*ratio),int(wg*ratio)), cv2.INTER_AREA)            
            diam_gong = int(diam_gong *ratio)
            img_gong = img_crop_diam(img_gong, diam_gong, marge=100)
            img_disk = img_crop_diam(img_disk, diam_disk, marge=100)
            gh, gw = img_gong.shape
            ih, iw = img_disk.shape
            
            # affiche iamge gong
            pixmap = QtGui.QPixmap()
            #pixmap.loadFromData(img_gong)
            img_gong = np.ascontiguousarray(img_gong)
            myqgong = QtGui.QImage(img_gong, gw, gh, img_gong.strides[0],QtGui.QImage.Format_Grayscale8 )
            pixmap= QtGui.QPixmap.fromImage(myqgong)
            
            lbl_w, lbl_h= (self.mygong.ui.gong_gongimg_lbl.width(),self.mygong.ui.gong_gongimg_lbl.height())
            pixmap.scaled(lbl_w, lbl_h,Qt.IgnoreAspectRatio)
            self.mygong.ui.gong_gongimg_lbl.setPixmap(pixmap)
            self.mygong.ui.gong_gongimg_lbl.adjustSize()

            #affiche image disk
            img_disk_8bits=(img_disk/256).astype(np.uint8)
            img_disk_8bits = np.ascontiguousarray(img_disk_8bits)
            myqimage = QtGui.QImage(img_disk_8bits, iw, iw ,iw, QtGui.QImage.Format_Grayscale8)
            pixmap2= QtGui.QPixmap.fromImage(myqimage)
            
            lbl_w, lbl_h= (self.mygong.ui.gong_myimg_lbl.width(),self.mygong.ui.gong_myimg_lbl.height()) 
            pixmap2.scaled(lbl_w, lbl_w,Qt.IgnoreAspectRatio)
            self.mygong.ui.gong_myimg_lbl.setPixmap(pixmap2)
            self.mygong.ui.gong_myimg_lbl.adjustSize()
          
            #print(self.mygong.ui.gong_myimg_lbl.size())
            
            try :
                angle_P=float(self.ui.ori_angP_text.text())
                if not self.ui.section_angP_auto_chk.isChecked() :
                    angle_P = 0
                # detection des inversions - uniquement sur disque entier et sur H-alpha
                inversion = gong_orientation_auto(img_gong, img_disk, diam_disk, angle_P)
                #print("Inversions : "+inversion)
                self.mygong.update_inversions(inversion)
                self.mygong.ui.finished.connect(self.ori_get_inversions)
            except:
                print('Erreur détection inversions')
                print('Affiche gong original')
                # affiche iamge gong
                pixmap = QtGui.QPixmap()
                #pixmap.loadFromData(img_gong)
                img_gong = np.ascontiguousarray(img_gong_orig)
                myqgong = QtGui.QImage(img_gong, gw, gh, img_gong.strides[0],QtGui.QImage.Format_Grayscale8 )
                pixmap= QtGui.QPixmap.fromImage(myqgong)

                
            # on affiche la fenetre
            self.mygong.show()
            
            
            
        else :
            print("Pas de fichier sélectionné")                                          
    
    def ori_get_inversions (self) :
        inv_list = self.mygong.get_inversions()
        #print(inv_list)
        
        self.ui.ori_inv_NS_chk.setChecked((inv_list[0]+self.ui.ori_inv_NS_chk.isChecked())%2)
        self.ui.ori_inv_EW_chk.setChecked((inv_list[1]+self.ui.ori_inv_EW_chk.isChecked())%2)
        
        
    def ori_angP(self):
        fits_dateobs = self.ui.ori_date_text.text()
        angP, paramB0, longL0, RotCarr = angle_P_B0(fits_dateobs)
        self.ui.ori_angP_text.setText(angP)
    
    def ori_grid_format(self):
        self.mygrid_dlg=grid_dialog()
        self.mygrid_dlg.ui.finished.connect(self.get_grid_format)

        
    def get_grid_format (self) :
        self.grid_color, self.flag_grid_disp = self.mygrid_dlg.get_format()
    
    def trame_sel_file (self, f1,f2) :
        msg = QMessageBox()
        msg.setWindowTitle("Choix de fichier")
        msg.setText("Quel fichier _mean voulez-vous utiliser ?")
        btn_a = msg.addButton("Mean", QMessageBox.AcceptRole)
        btn_b = msg.addButton("Mean_Start", QMessageBox.AcceptRole)
        msg.addButton(QMessageBox.Cancel)
        
        msg.exec()
        
        if msg.clickedButton() == btn_a:
            print(f1)
            f = f1
            ext=''
        elif msg.clickedButton() == btn_b:
            print(f2)
            f = f2
            ext='start'
        else:
            print("Action annulée")
            f=''
            ext=''
        
        return f,ext

    
    def trame_mean_img (self):
        # version image redressee
        # Lecture et affiche image disque brut
        if  self.serfiles :
            self.get_Flags_tab()
            poly=self.polynome
            trame_dir=self.working_dir+os.sep+"Complements"
            serfile = self.serfiles[0]
            self.Flags['Couronne'] = self.ui.free_corona_chk.isChecked()
            # extrait nom racine
            base=os.path.basename(serfile)
            if base=='':
                print(self.tr('Erreur de nom de fichier : ')+serfile)
            else :
                
                # basefich: nom du fichier ser sans extension et sans repertoire
                self.basefich='_'+os.path.splitext(base)[0]
                f1 = self.basefich+"_mean.fits"
                if self.Flags['WEAK'] and not self.Flags['Couronne']:
                    f2 = self.basefich+"_mean_start.fits"
                    f, m_ext = self.trame_sel_file (f1,f2)
                else :
                    m_ext=''
                    f= f1
                ImgFile = trame_dir+os.sep + f
    
                
                if ImgFile != None :
                    try :
                        
                        mean_trame, header= self.read_fits_image(ImgFile)
                        mean_trame=np.flipud(mean_trame)
                                    
                        # redresse raies pour faciliter la prise de position de la raie
                        try :
                            image_mean_redresse = translate_img(mean_trame, poly)
                        except :
                            image_mean_redresse = np.copy(mean_trame)
                            print("pb...")
                        
                        # affiche l'image redressée
                        self.my_trame = trame_img(self)
                        self.my_trame.const=poly[2]
                        self.my_trame.show()
                        self.my_trame.on_click.connect(self.update_constante)
                        if m_ext == 'start' : 
                            self.my_trame.auto.connect(self.update_constante)
                        self.my_trame.set_img(image_mean_redresse)
                        
                        
                        flag_save_mean=False
                        
                        if flag_save_mean :
                            #sauvegarde du fits en _check
                            mih, miw = image_mean_redresse.shape()
                            WorkDir=os.path.dirname(ImgFile)+"/"
                            base=os.path.basename(ImgFile)
                            basefich=os.path.splitext(base)[0].split('_mean')[0]
                            header['NAXIS1']=miw
                            fname=WorkDir+basefich+'_check.fits'
                            self.save_fits_image (fname,image_mean_redresse,header, 16)
                            print(self.working_dir+basefich+'_check.fits')
                        
                        
                        #print("image affichée")
                    except :
                        if not file_exist(ImgFile) :
                            print(self.tr('Fichier non trouvé : ' + ImgFile))
                        else :
                            pass
        else :
            print(self.tr('Aucun fichier sélectionné'))
    
    def console_clear (self):
        self.ui.log_edit.setText('INTI Version : ' + str(self.version))
        
    def update_constante(self):
        self.Flags["POL"]=False
        self.Flags["WEAK"]=False
        app_tab= self.ui.tab_main.tabText(self.ui.tab_main.currentIndex())
        if app_tab =='Helium-Couronne-Libre' or app_tab == 'Helium-Coronna Free':
            self.Flags['WEAK'] =True
        if app_tab == 'Magnétogramme' or app_tab == "Magnetogram":
            self.Flags['POL'] =True
        
        const=self.my_trame.coord
        
        if self.Flags['WEAK'] :
            self.ui.free_shift_text.setText(str(const))
        if self.Flags['POL']:
            self.ui.magnet_shift_text.setText(str(const))
    
    def dop_profil (self) :
        try :
            base=os.path.basename(self.serfiles[0])
            basefich=os.path.splitext(base)[0]
            pro_file=self.working_dir+os.path.sep+"Complements"+os.path.sep+'_'+basefich+'.dat'
            print(self.tr('Fichier profil : '), pro_file)
            if basefich != '' :
                try :
                    pro_data=np.loadtxt(pro_file)
                    self.mypro=profil_wnd()
                    self.mypro.set_profil(pro_data)
                    self.mypro.show()
                    self.mypro.set_title(basefich)
                except :
                    pass
        except :
            print(self.tr('Aucun fichier sélectionné'))
            
    def cfg_files_to_save (self) :
        self.myconfig_dlg=config_dialog()
        self.myconfig_dlg.set_files_to_save(self.files_to_save)
        self.myconfig_dlg.ui.finished.connect(self.get_files_to_save)

    def crop_clicked (self):
        if self.ui.section_autocrop_chk.isChecked():
            self.ui.inti_h20_radio.setEnabled(True)
        else :
            self.ui.inti_h20_radio.setChecked(False)
            self.ui.inti_h20_radio.setEnabled(False)


    def corona_clicked(self) :
        # si on passe en mode couronne cela n'a de sens que si le mode poly auto est actif et desactive helium
        if self.ui.free_corona_chk.isChecked() :
            self.ui.free_poly_auto_chk.setChecked(True)
            self.ui.free_trans_helium_chk.setChecked(False)
            
    
    def get_files_to_save (self) :
        self.files_to_save = self.myconfig_dlg.get_files_to_save()
        #print(self.files_to_save)
    
    
    def read_fits_image(self, file_name):
        with fits.open(file_name,memmap=False) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            
            if np.min(data) < 0 :
                offset=32767
                data=data+offset
            
            data=np.flipud(data)
            
        return data, header
    
    def save_fits_image (self, name, data, hdr, nb_bytes):
        if nb_bytes == 16 :
            data = np.array( data, dtype='uint16')   # conversion en flottant 32 bits
        if nb_bytes == 32 :
            data = data.astype(np.float32)   # conversion en flottant 32 bits
        
        DiskHDU = fits.PrimaryHDU(data, hdr)
        DiskHDU.writeto(name, overwrite='True')
        
    def save_png_image (self, name, data):
        data = np.array( data, dtype='uint16')   # conversion en flottant 32 bits
        cv2.imwrite(name, data)
    
    
    
    
    def name_to_png (self,filename):
        # retire extension fits pour ajouter extension png
        name_png=os.path.splitext(filename)[0]+'.png'
        #print(name_png)
        return name_png

    def file_exist(name):    
        if not os.path.exists(name):
            print('ERROR: File ' + name + ' not found.')
            sys.exit(1)
            
    def short_name(self,file_name):
        #f_short=file_name.split('/')[-1]
        f_short=os.path.split(file_name)[1]
        return f_short
    
    def racine_name(self,file_name):
        f_short=os.path.split(file_name)[1]
        f_racine=f_short.split('.')[0]
        return f_racine
    
    def get_dirpath(self,file_name):
        dir_path=os.path.split(file_name)[0]
        #print('dir : ', dir_path)
        return dir_path
    
    def get_extension(self, file_name) :
        pos=file_name.rfind('.')
        #ext=os.path.split(file_name)[1].split('.')[1]
        ext = file_name[pos+1:]
        return ext

    def get_basename(self, file_name) :
        pos=file_name.rfind('.')
        #ext=os.path.split(file_name)[1].split('.')[1]
        basename = file_name[:pos]
        return basename
    


# ----------------------------------------------------------------------------
# class new window to display image floating
#-----------------------------------------------------------------------------

class img_wnd(QMainWindow) :
    
    on_ferme = Signal()
    
    def __init__(self, myzoom_wnd, flag_zoom_on, mainwindow):
       
        #super().__init__(parent)
        super(img_wnd, self).__init__()
        self.mainwindow = mainwindow
        
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('img_qt.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        #self.ui.closeAll_btn.clicked.connect(self.mon_closeEvent)
        self.ui.closeAll_btn.clicked.connect(self.emit_on_ferme)
        
        self.flag_zoom_on=flag_zoom_on
        if self.flag_zoom_on : 
            self.myzoom_wnd=myzoom_wnd
        

        # set icon application
        self.ui.setWindowIcon(QtGui.QIcon(resource_path("inti_logo.png")))
    
        
        # recupere la position
        self.read_settings()
        
        self.ui.inti_view.ui.roiBtn.hide()
        self.ui.inti_view.ui.menuBtn.hide()
        self.ui.inti_view.sigTimeChanged.connect(self.frame_changed)
        self.ui.inti_view.scene.sigMouseMoved.connect(self.on_mouse_move)
        self.ui.inti_view.view.setBackgroundColor((20,20,20))
        
        # Create a keyboard shortcut Ctrl+M        #shortcut = QtGui.QKeySequence(Qt.CTRL + Qt.Key_M)
        shortcut = QtGui.QKeySequence("Escape")
        self.shortcut = QtGui.QShortcut(shortcut, self.ui)
        self.shortcut.activated.connect(self.emit_on_ferme)
        self.ui.save_btn.clicked.connect(self.save_file)

        
    def emit_on_ferme(self) :
        #print("emit")
        self.on_ferme.emit()
    
    def set_img(self, img_data) :
        img_proc = np.fliplr(np.rot90(img_data, 3))
        self.ui.inti_view.setImage(img_proc)
        #self.ui.inti_view.setLevels(0,65535)
    
    def set_vol(self, img_data) :
        #img_proc = np.fliplr(np.rot90(img_data, 3))
        self.ui.inti_view.setImage(img_data)

    def save_file(self):
        new_name,_=QFileDialog.getSaveFileName(self, self.tr("Sauver fichier png"), self.file_name, self.tr("Fichiers png (*.png);;Tous les fichiers (*)"))
        if new_name :
            myimage=self.ui.inti_view.image
            myimage=np.flipud(np.rot90(myimage))
            if len(myimage.shape)==3 :
                #ajuste les seuils
                levels = self.ui.inti_view.getLevels()
                sbas,shaut = levels
                if shaut != sbas :
                    myimage = np.clip((myimage.astype(np.float32)-sbas)/(shaut-sbas),0,1)
                    myimage = (myimage*255).astype(np.uint8)
                myimage=cv2.cvtColor(myimage, cv2.COLOR_BGR2RGB)
                cv2.imwrite(new_name, myimage)
            else :
                # ajustment des seuils
                levels = self.ui.inti_view.getLevels()
                sbas,shaut = levels
                if shaut != sbas :
                    myimage = np.clip((myimage-sbas)/(shaut-sbas)*65535,0,65535).astype(np.uint16)
                cv2.imwrite(new_name, myimage)
        # on duplique le fichier log avec le nom du nouveau fichier
        file_log_name = get_baseline(os.path.splitext(self.file_name)[0])+'_log.txt'
        new_log_name = os.path.splitext(new_name)[0]+'_log.txt'
        if file_exist(file_log_name) :
            shutil.copy(file_log_name, new_log_name)
        
        
    def on_mouse_move (self, pos):
        if self.ui.inti_view.imageItem.sceneBoundingRect().contains(pos) :
            mouse_point= self.ui.inti_view.view.mapSceneToView(pos)
            x,y =int(mouse_point.x()), int(mouse_point.y())
            
            if 150 <= x < self.ui.inti_view.image.shape[0]-150 and 150<= y < self.ui.inti_view.image.shape[1]-150:
                
                img = self.ui.inti_view.image
                img=np.fliplr(np.rot90(img,3))
                #y=img.shape[0]-y
                img=img[y-150:y+150,x-150:x+150]
                try :
                    if self.flag_zoom_on :
                        self.myzoom_wnd.update_img(img)
                except :
                    pass
                
            else:
                pass
                #self.myzoom_wnd.hide()
            
            # on visualise l'intensité
            if 0 <= x < self.ui.inti_view.image.shape[0] and 0<= y < self.ui.inti_view.image.shape[1] :
                
                # pour affichage premiere ligne n'est pas 0 mais 1
                msg="x : "+str(x+1)+' , y : '+str(y+1)
                #print(msg)
                pix_value = self.ui.inti_view.image[x,y]
                try :
                    if len(pix_value) != 1:
                        msg=msg+' , R : '+str(int(pix_value[0]))+ ' , G : '+str(int(pix_value[1]))+' , B : '+str(int(pix_value[2]))
                except :
                    msg=msg+' , I : '+str(int(pix_value))
                self.ui.statusbar.showMessage(self.nomfich+'  '+msg)
            else:
                #print ("mouse out of bounds")
                self.ui.statusbar.showMessage(self.nomfich)
                    
                    
        
    def frame_changed(self) :
        test=False
        if test :
            trame_index=self.ui.inti_view.currentIndex
            vol_img=self.ui.inti_view.image
            #my_trame = vol_img[trame_index]
        else :
            pass
    
    def show(self) :
        self.ui.show()

    def mon_closeEvent (self,event):
        
        self.write_settings()
        self.on_ferme.emit()
        #self.ui.close()
        #self.mainwindow.close_all_subwindows()
        #self.ui.close()

    def read_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        self.ui.restoreGeometry(settings.value("ImgWindow/geometry"))
        self.ui.restoreState(settings.value("ImgWindow/windowState"))
    
    def write_settings(self) :
        # sauve settings
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("ImgWindow/geometry", self.ui.saveGeometry())
        settings.setValue("ImgWindow/windowState", self.ui.saveState())
        
    def set_title (self, title ):
        self.nomfich = title
        self.ui.setWindowTitle(title)
        self.ui.statusbar.showMessage(title)
        
    def set_file_name (self, file_name ):
        self.file_name = file_name
        
        
    def set_pos (self, posx,posy) :
        self.ui.move(posx,posy)


# ----------------------------------------------------------------------------
# class new window to display zoom
#-----------------------------------------------------------------------------

class zoom_wnd(QWidget) :
    
    def __init__(self):
        super(zoom_wnd, self).__init__()
       
        
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('zoom.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        self.ui.setWindowFlags(Qt.Tool | Qt.WindowStaysOnTopHint)
        self.ui.setWindowTitle('Zoom 1:1')
        self.ui.setFixedSize(300,300)
        self.ui.move(200,100)
        
    def update_img (self, image):
       
        try :
            ih,iw,nbplan= np.array(image).shape
        except :
            nbplan=1
        if nbplan == 1 :
            # transforme en 8 bits
            img_8bits= (image/256).astype(np.uint8)
            h,w =img_8bits.shape
            img_8bits = np.ascontiguousarray(img_8bits)
            image=QtGui.QImage(img_8bits.data, w, h ,w ,QtGui.QImage.Format_Grayscale8).copy()
            pixmap=QtGui.QPixmap.fromImage(image)
            self.ui.zoom_lbl.setPixmap(pixmap)
        else :
            icolor=np.zeros([300,300,4], dtype='uint16')
            icolor[:,:,0]=image[:,:,0]
            icolor[:,:,1]=image[:,:,1]
            icolor[:,:,2]=image[:,:,2]
            icolor[:,:,3]= 65535 # couche alpha
            im=np.require (icolor, np.uint16, 'C') # contiguous
            h, w, nbplan =im.shape
            qicolor=QtGui.QImage(im.data,w,h,im.strides[0] ,QtGui.QImage.Format_RGBX64).copy()
            pixmap=QtGui.QPixmap.fromImage(qicolor)
            self.ui.zoom_lbl.setPixmap(pixmap)
        
    def show(self):
        self.ui.show()


class gong_wnd(QDialog) :

    
    def __init__(self):
        super(gong_wnd, self).__init__()
       
        
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        ui_file_name=resource_path('gong.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        # geometry écrans
        screen_geom = QtGui.QGuiApplication.primaryScreen().availableGeometry()
        self.dpr = QtGui.QGuiApplication.primaryScreen().devicePixelRatio()
        self.myscreen_w = int(screen_geom.right()*self.dpr)
        self.myscreen_h = int(screen_geom.bottom()*self.dpr)
        
        # settings
        flag_reset=False
        settings=QSettings("Desnoux Buil", "inti_qt")
        if settings.value("GongWindow/geometry") and not flag_reset:

            self.ui.restoreGeometry(settings.value("GongWindow/geometry"))
        else :
            w_w=self.myscreen_w*0.4
            #r=self.myscreen_h/self.myscreen_w
            w_h = (w_w*0.5)+9
            if w_h > self.myscreen_h :
                w_w=w_h*2
            self.ui.resize(int(w_w), int(w_h))
        
        
        # connect signaux boutons
        self.ui.gong_GD_btn.clicked.connect(self.gong_gd)
        self.ui.gong_HB_btn.clicked.connect(self.gong_hb)
        self.ui.gong_blink_btn.clicked.connect(self.gong_blink)
        #self.ui.gong_apply_btn.clicked.connect(self.update_ori_image)
        
        # Sauvegarder quand la boîte se termine, quelle que soit la cause
        self.ui.finished.connect(self.save_settings)
        

        
        # init nb inversions
        self.nb_ns=0
        self.nb_ew = 0 
        
        self.ui.move(10,10)
        
    def update_inversions(self, inv) :
        self.inv = inv
        self.ui.gong_inversions_lbl.setText('Inversions proposées : '+ inv)
    
    def update_ori_image (self):
        inv=self.inv
        if self.inv != 'None' :
            mypixmap= self.ui.gong_myimg_lbl.pixmap()
            if inv == 'NS' :
               flipped = mypixmap.transformed(QtGui.QTransform().scale(1,-1))
               self.nb_ns=self.nb_ns+1
            if inv == 'EW' :
               flipped = mypixmap.transformed(QtGui.QTransform().scale(-1,1))
               self.nb_ew=self.nb_ew+1
            if inv == 'NS-EW' :
               flipped = mypixmap.transformed(QtGui.QTransform().scale(-1,-1))
               self.nb_ns=self.nb_ns+1
               self.nb_ew=self.nb_ew+1
            
            self.ui.gong_myimg_lbl.setPixmap(flipped)
    
    def gong_gd (self):
        mypixmap= self.ui.gong_myimg_lbl.pixmap()
        flipped = mypixmap.transformed(QtGui.QTransform().scale(-1,1))
        self.nb_ew=self.nb_ew+1
        self.ui.gong_myimg_lbl.setPixmap(flipped)
    
    def gong_hb (self):
        mypixmap= self.ui.gong_myimg_lbl.pixmap()
        flipped = mypixmap.transformed(QtGui.QTransform().scale(1,-1))
        self.nb_ns=self.nb_ns+1
        self.ui.gong_myimg_lbl.setPixmap(flipped)
        
    def get_inversions (self) :
        inv_to_return = [self.nb_ns%2, self.nb_ew%2]
        return inv_to_return
    
    def gong_blink (self) :
        self.pixmap_disk = self.ui.gong_myimg_lbl.pixmap()
        if self.ui.gong_blink_btn.isChecked() :
            self.ui.gong_gongimg_lbl.setPixmap (self.pixmap_disk)
        else :
            self.ui.gong_gongimg_lbl.setPixmap (self.pixmap_gong)
    
    def show(self):
        self.ui.show()
        self.pixmap_gong = self.ui.gong_gongimg_lbl.pixmap()
        
    def save_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("GongWindow/geometry", self.ui.saveGeometry())

        

# ----------------------------------------------------------------------------
# class new window to display profil
#-----------------------------------------------------------------------------

class profil_wnd(QMainWindow) :
    
    def __init__(self, parent=None):
       
        #super().__init__(parent)
        super(profil_wnd, self).__init__()
                
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('profil_qt.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()

        # set icon application
        self.ui.setWindowIcon(QtGui.QIcon(resource_path("inti_logo.png")))
        
        #app.aboutToQuit.connect(self.close)
        
        # init
        self.label = pg.TextItem(text='',color='black', anchor=(1,1))
        #const=0
        
        # recupere la position
        self.read_settings()
    
    def set_profil(self, pro_data) :
        pro=pro_data[:,1]
        lamb=pro_data[:,0]  
        # calcul centre de la raie
        centre_raie = np.argmin(pro)
        
        lamb=lamb-centre_raie
        # affiche
        self.ui.profil_view.setBackground('w')
        pen=pg.mkPen(color='blue',width=1.5)
        self.myplot = self.ui.profil_view.plot(lamb, pro, pen=pen, symbol = 'o',symbolSize=5,name='profil', maxTickLength=-100)
        self.myplot.setCurveClickable(True)
        self.myplot.sigClicked.connect(self.handle_curve_click)
        self.ui.profil_view.getPlotItem().addItem(self.label)

        
        self.curr_data = self.myplot.curve.getData()
        
    def handle_curve_click(self, curve,event):
        #print('curve scene : ',event.pos())
        #self.curr_data = curve.getData()
        mx=event.pos().x()
        my=event.pos().y()
        x=round(mx)
        #msg='X : '+ str(x)
        #self.ui.statusbar.showMessage(msg)
        self.label.setText(str(x)+' pixels')
        self.label.setPos(mx,my)
    
    
    def show(self) :
        self.ui.show()

    def closeEvent (self,event):
        self.write_settings()
        self.ui.close()

    def read_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        self.ui.restoreGeometry(settings.value("ProWindow/geometry"))
        self.ui.restoreState(settings.value("ProWindow/windowState"))
    
    def write_settings(self) :
        # sauve settings
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("ProWindow/geometry", self.ui.saveGeometry())
        settings.setValue("ProWindow/windowState", self.ui.saveState())
        
    def set_title (self, title ):
        self.ui.setWindowTitle(title)
        self.ui.statusbar.showMessage(title)
        
    def set_pos (self, posx,posy) :
        self.ui.move(posx,posy)
        
    def on_mouse_move (self) :
        pass
        
# ----------------------------------------------------------------------------
# class Widget trame avec click souris
#-----------------------------------------------------------------------------
class trame_img(QWidget) :
    
    on_click = Signal()
    auto =  Signal()
    
    def __init__(self, parent=None):
       
        #super().__init__(parent)
        super(trame_img, self).__init__()
                
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('trame.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()

        # set icon application
        self.ui.setWindowIcon(QtGui.QIcon(resource_path("inti_logo.png")))
        
        # initialisation 
        self.ui.trame_view.ui.roiBtn.hide()
        self.ui.trame_view.ui.menuBtn.hide()
        #const=0
        
        # gestion des signaux
        self.ui.trame_view.scene.sigMouseMoved.connect(self.on_mouse_move)
        self.ui.trame_view.scene.sigMouseClicked.connect(self.on_mouse_click)
        self.ui.trame_ok_btn.clicked.connect(self.on_ok_click)
        

        # recupere la position
        #self.read_settings()
        
    def set_img(self, img_data) :
        # modif 11 juillet 2025 - ne plus faire de zoom
        img_data_pro=np.copy(img_data)
        img_data = np.fliplr(np.rot90(img_data, 3))
        # en zoom 1:1
        #curr_scale=self.ui.trame_view.view.viewPixelSize()[0] # 1 pixel screen en pixel image
        #zo=1/curr_scale
        #zoom=(zo,zo)
        
        # affiche image 
        self.ui.trame_view.setImage(img_data)
        self.ui.trame_view.setLevels (0, np.percentile(img_data,99.5))
        #self.ui.trame_view.view.scaleBy(zoom)
        try :
            # mais on zoom pour meilleure precision 
            ih, iw = img_data_pro.shape
            mih = ih//2
            zone =mih//4
            img_data2=img_data_pro[mih-zone:mih+zone, :]
            img_data2 = np.fliplr(np.rot90(img_data2, 3))
            
            pro,peak = self.compute_profil(img_data2)
            #plt.plot(pro)
            #plt.show()
            if self.const <0 : 
                self.coord = int(peak+self.const)
            else :
                self.coord = int(peak-self.const)
            self.auto.emit()
        except :
            pass
    
        
    def compute_profil (self, img_data):  

        profil=np.mean(img_data, axis=1)
        peak=np.argmax(profil)

        return profil, peak
    
    def show(self) :
        self.ui.show()

    def closeEvent (self,event):
        self.write_settings()
        self.ui.close()

    def read_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        self.ui.restoreGeometry(settings.value("ImgTrame/geometry"))
        self.ui.restoreState(settings.value("ImgTrame/windowState"))
    
    def write_settings(self) :
        # sauve settings
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("ImgTrame/geometry", self.ui.saveGeometry())
        settings.setValue("ImgTrame/windowState", self.ui.saveState())
        
    def set_title (self, title ):
        self.ui.setWindowTitle(title)

        
    def set_pos (self, posx,posy) :
        self.ui.move(posx,posy)
        
    def on_mouse_move(self,pos) :
        if self.ui.trame_view.imageItem.sceneBoundingRect().contains(pos) :
            mouse_point= self.ui.trame_view.view.mapSceneToView(pos)
            x,y =int(mouse_point.x()), int(mouse_point.y())
            
            if 0 <= x < self.ui.trame_view.image.shape[0] and 0<= y < self.ui.trame_view.image.shape[1] :
                # pour affichage premiere ligne n'est pas 0 mais 1
                msg="x : "+str(x+1)+' , y : '+str(y+1)
                #print(msg)
                pix_value = self.ui.trame_view.image[x,y]
                msg=msg+' , I : '+str(int(pix_value))
                #print(pix_value)
                self.ui.trame_move_lbl.setText(msg)
            else:
                #print ("mouse out of bounds")
                self.ui.trame_move_lbl.setText('')
                
    def on_mouse_click (self, ev) :

        pos=ev.pos()
        if self.ui.trame_view.imageItem.sceneBoundingRect().contains(pos) :
            
            mouse_point= self.ui.trame_view.view.mapSceneToView(pos)
            x,y =int(mouse_point.x()), int(mouse_point.y())
            if self.const <0 : 
                p=int(x+self.const)
            else :
                p=int(x-self.const)
            if cfg.LG == 1 :
                msg='X : '+str(x+1)+' ; Décalage : '+str(p)
            else :
                msg='X : '+str(x+1)+' ; Shift : '+str(p)
            self.ui.trame_click_lbl.setText(msg)
            self.coord=p
            #self.on_click.emit()
        else:
            #print ("mouse out of bounds")
            self.ui.trame_click_lbl.setText('')
            
    def on_ok_click(self) :
        self.on_click.emit()


class calc_dialog(QDialog):
    def __init__(self, foc, ang):
        super().__init__()
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        #loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('calc.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()
        self.wav_name='Ha'
        self.wav_dict={'Ha':6562.762,'Ca':3968.469,'He':5877.3, 'Fe':5302.86}
        
        
        validator_4dec = DotOnlyDoubleValidator()
        validator_4dec.setLocale(QLocale(QLocale.English))
        validator_4dec.setDecimals(4)           # Nombre max de décimales autorisées
        validator_4dec.setNotation(QtGui.QDoubleValidator.StandardNotation)
        
              
        self.ui.calc_ang_btn.clicked.connect(self.calc_valid_angtopix)
        self.ui.calc_pix_btn.clicked.connect(self.calc_valid_pixtoang)
        self.ui.calc_lamb_text.setText(str(self.wav_dict[self.wav_name]))
        self.ui.buttonGroup.idClicked.connect(self.wave_radio_clicked)
        self.ui.buttonGroup_2.idClicked.connect(self.shg_radio_clicked)
        self.ui.shg_popup_btn.clicked.connect(self.shg_custom)
        self.shg_foc_camera_other = foc
        self.shg_angle_other = ang
        if foc !='0' and ang!='0' :
            self.ui.shg_autre.setEnabled(True)
        else :
            self.ui.shg_autre.setEnabled(False)
                 
        self.ui.calc_lamb_text.setValidator(validator_4dec)
        self.ui.calc_size_text.setValidator(validator_4dec)
        self.ui.calc_ang_text.setValidator(validator_4dec)
        self.ui.calc_pix_text.setValidator(validator_4dec)
        self.ui.calc_bin_text.setValidator(QtGui.QIntValidator(1, 10))
        
        # update
        self.shg_radio_clicked()
     
        
        
    def show(self) :
        self.ui.show()
        
    def shg_radio_clicked(self ):
        if self.ui.shg_solex.isChecked() :
            self.shg_angle = 34
            self.shg_foc_camera = 125
        if self.ui.shg_autre.isChecked() :
            self.shg_foc_camera = self.shg_foc_camera_other
            self.shg_angle = self.shg_angle_other
            
    def shg_custom (self) :
        pos = QtGui.QCursor.pos()
        self.theshg = shg_popup(self)
        self.theshg.valeurs_onclose.connect(self.valeurs_received)
        self.theshg.ui.move(pos)
        self.theshg.ui.show()
    
    def valeurs_received (self, foc_cam, angle_res) :
        self.shg_foc_camera_other = foc_cam
        self.shg_angle_other = angle_res
        if foc_cam != '0' or angle_res != '0':
            self.ui.shg_autre.setEnabled(True)
        else :
            self.ui.shg_autre.setEnabled(False)
        
       
  
    def wave_radio_clicked(self ):
        if self.ui.wave_btn_ha.isChecked() :
            self.wav_name='Ha'
        if self.ui.wave_btn_he.isChecked() :
            self.wav_name='He'
        if self.ui.wave_btn_ca.isChecked() :
            self.wav_name='Ca'
        if self.ui.wave_btn_fe.isChecked() :
            self.wav_name='Fe'
        
        self.ui.calc_lamb_text.setText(str(self.wav_dict[self.wav_name]))
        
    def calc_valid_angtopix (self):
        self.wave=float(self.ui.calc_lamb_text.text())*1e-7
        alpha=np.degrees(np.arcsin((2400*self.wave)/(2*np.cos(np.radians(17)))))+17
        beta=alpha-int(self.shg_angle)
        size_pix_cam=float(self.ui.calc_size_text.text())
        bin_cam=float(self.ui.calc_bin_text.text())
        val_toconvert=float(self.ui.calc_ang_text.text())
        disp = 1e7 * size_pix_cam* bin_cam*np.cos(np.radians(beta)) / 2400 / int(self.shg_foc_camera)
        val_ang_onepix= disp 
        #print(val_ang_onepix)
        resultat=  "{:.3f}".format(val_toconvert / val_ang_onepix)
        self.ui.calc_pix_text.setText(str(resultat))
        self.ui.calc_disp_lbl.setText("{:.3f}".format(disp))
        
    def calc_valid_pixtoang (self):           
         self.wave=float(self.ui.calc_lamb_text.text())*1e-7
         alpha=np.degrees(np.arcsin((2400*self.wave)/(2*np.cos(np.radians(17)))))+17
         beta=alpha-int(self.shg_angle)
         size_pix_cam=float(self.ui.calc_size_text.text())
         bin_cam=float(self.ui.calc_bin_text.text())
         val_toconvert=float(self.ui.calc_pix_text.text())
         disp = 1e7 * size_pix_cam* bin_cam*np.cos(np.radians(beta)) / 2400 /  int(self.shg_foc_camera)
         val_ang_onepix= disp 
         #print(val_ang_onepix)
         resultat= "{:.3f}".format(val_toconvert * val_ang_onepix)
         self.ui.calc_ang_text.setText(str(resultat))
         self.ui.calc_disp_lbl.setText("{:.3f}".format(disp))
    

class grid_dialog(QDialog):
    def __init__(self, parent=None):
        super().__init__()
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        #loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('grid.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()
                
        self.ui.grid_ok_btn.clicked.connect(self.ui.accept)
        self.ui.show()
    
    def get_format(self) :
        flag_gradu = self.ui.grid_grad_disp_chk.isChecked()
        grid_color = self.ui.grid_color_combo.currentText()
        if grid_color =='noir' :
            grid_color='black'
        if grid_color =='jaune' :
            grid_color='yellow'
        return grid_color, flag_gradu
        
    def show(self) :
        self.ui.show()   
    
    def set_pos (self, posx,posy) :
        self.ui.move(posx,posy)
        
        
# ----------------------------------------------------------------------------
# class dialog entete BASS2000
#-----------------------------------------------------------------------------

class entete_dialog(QDialog):
    def __init__(self, entete, parent=None):
        super().__init__()
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        #loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('param.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()

        # validation des entrées
        #validator_2dec = QtGui.QDoubleValidator()
        validator_2dec = DotOnlyDoubleValidator()
        validator_2dec.setLocale(QLocale(QLocale.English))
        validator_2dec.setDecimals(2)           # Nombre max de décimales autorisées
        validator_2dec.setNotation(DotOnlyDoubleValidator.StandardNotation)
        validator_2dec.setRange(-180.0, 180.0)  # Bornes min et max
        
        self.ui.db_lat_text.setValidator(validator_2dec)
        self.ui.db_long_text.setValidator(validator_2dec)
        self.ui.BASS_pixel_text.setValidator(validator_2dec)
        self.ui.BASS_focinstru_text.setValidator(QtGui.QIntValidator(0, 9999))
        self.ui.BASS_diaminstru_text.setValidator(QtGui.QIntValidator(0, 9999))
        self.ui.BASS_objcollim_text.setValidator(QtGui.QIntValidator(0, 1000))
        self.ui.BASS_objcam_text.setValidator(QtGui.QIntValidator(0, 1000))
        self.ui.BASS_reseau_text.setValidator(QtGui.QIntValidator(0, 9999))
        self.ui.BASS_ordre_text.setValidator(QtGui.QIntValidator(1, 999))
        self.ui.BASS_ang_text.setValidator(QtGui.QIntValidator(0, 360))
        self.ui.BASS_fentelong_text.setValidator(validator_2dec)
        self.ui.BASS_fentelarg_text.setValidator(QtGui.QIntValidator(0, 50))
        self.ui.BASS_diaph_text.setValidator(QtGui.QIntValidator(0, 9999))
        
        # signaux
        self.ui.BASS_ok_btn.clicked.connect(self.ok)
        self.ui.db_wave_text.currentIndexChanged.connect(self.db_print_wave)
        
        # initialisations
        self.list_wave=[['Manual','Ha','Ha2cb','Cah','Cah1v','Cak','Cak1v','HeID3'],[0,6562.762,6561.232,3968.469,3966.968,3933.663,3932.163,5877.3]]
        self.ui.db_wave_text.addItems(self.list_wave[0])
        
        # met à jour les champs collectés dans l'entete
        self.ui.db_observer_text.setText (entete['observer'])
        self.ui.db_lat_text.setText (entete['lat'])
        self.ui.db_long_text.setText(entete['long'])
        self.ui.db_instru_text.setText(entete['instru'])
        self.ui.db_contact_text.setText(entete['contact'])
        self.ui.db_wave_text.setCurrentIndex(entete['waveID'])
        self.ui.db_binning_combo.setCurrentIndex(entete['binID'])
        self.ui.BASS_camera_text.setText(entete["camera"])     
        self.ui.BASS_pixel_text.setText(entete["pixel"] )
        self.ui.BASS_focinstru_text.setText(entete["focinstru"] )
        self.ui.BASS_diaminstru_text.setText(entete["diaminstru"] )
        self.ui.BASS_spectro_text.setText(entete["spectro"] )
        self.ui.BASS_objcollim_text.setText(entete["objcollim"] )
        self.ui.BASS_objcam_text.setText(entete["objcam"] )
        self.ui.BASS_reseau_text.setText(entete["reseau"] )
        self.ui.BASS_ordre_text.setText(entete["ordre"] )
        self.ui.BASS_ang_text.setText(entete["angle"] )
        self.ui.BASS_fentelong_text.setText(entete["fentelong"] )
        self.ui.BASS_fentelarg_text.setText(entete["fentelarg"] )
        self.ui.BASS_diaph_text.setText(entete["diaph"] )
        self.ui.BASS_nd_text.setText(entete["dfilter"] )


        #self.ui.grid_ok_btn.clicked.connect(self.ui.accept)
        self.ui.show()
        
    def db_print_wave (self) :
        index= self.ui.db_wave_text.currentIndex()
        print(self.list_wave[1][index])

    def get_parameters(self):
        
        myparam={"observer" : self.ui.db_observer_text.text(),
                "contact" : self.ui.db_contact_text.text(),
                "lat": self.ui.db_lat_text.text(),
                "long": self.ui.db_long_text.text(),
                "instru" : self.ui.db_instru_text.text(),               
                "waveID" : self.ui.db_wave_text.currentIndex(),
                "binID" : self.ui.db_binning_combo.currentIndex(),
                "camera" : self.ui.BASS_camera_text.text(),             
                "pixel" : self.ui.BASS_pixel_text.text(),
                "focinstru" : self.ui.BASS_focinstru_text.text(),
                "diaminstru" : self.ui.BASS_diaminstru_text.text(),
                "spectro" : self.ui.BASS_spectro_text.text(),
                "objcollim" : self.ui.BASS_objcollim_text.text(),
                "objcam" : self.ui.BASS_objcam_text.text(),
                "reseau" : self.ui.BASS_reseau_text.text(),
                "ordre" : self.ui.BASS_ordre_text.text(),
                "angle" : self.ui.BASS_ang_text.text(),
                "fentelong" :self.ui.BASS_fentelong_text.text(),
                "fentelarg" :self.ui.BASS_fentelarg_text.text(),
                "diaph" : self.ui.BASS_diaph_text.text(),
                "dfilter" : self.ui.BASS_nd_text.text()
                }
        
        return myparam
    
    def ok (self):
        self.ui.close()

   
            

# ----------------------------------------------------------------------------
# class dialog fichiers à enregistrer
#-----------------------------------------------------------------------------

class config_dialog(QDialog):
    def __init__(self, parent=None):
        super().__init__()
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        #loader.registerCustomWidget(ImageView)
        ui_file_name=resource_path('config_save.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        ui_file.close()
                
        #self.ui.grid_ok_btn.clicked.connect(self.ui.accept)
        self.ui.show()

        
    def set_files_to_save(self, flags) :
        #self.ui.cfg_disk_chk.setChecked(flags['disk']) 
        self.ui.cfg_raw_chk.setChecked(flags['raw'])
        self.ui.cfg_color_chk.setChecked(flags['color']) 
        self.ui.cfg_inv_chk.setChecked(flags['inv'])
        self.ui.cfg_mix_chk.setChecked(flags['mix']) 
        #flags['recon'] = self.ui.cfg_recon_chk.isChecked()
        #flags['rawfits'] = self.ui.cfg_rawfits_chk.isChecked()
        self.ui.cfg_diff_chk.setChecked(flags['diff']) 
        self.ui.cfg_sum_chk.setChecked(flags['sum']) 
        self.ui.cfg_x0x1x2_chk.setChecked(flags['x0x1x2']) 
        self.ui.cfg_dop_chk.setChecked(flags['dop']) 
        
    def get_files_to_save(self) :
        flags={}
        flags['raw'] = self.ui.cfg_raw_chk.isChecked()
        flags['color'] = self.ui.cfg_color_chk.isChecked()
        flags['inv'] = self.ui.cfg_inv_chk.isChecked()
        flags['mix'] = self.ui.cfg_mix_chk.isChecked()

        flags['diff'] = self.ui.cfg_diff_chk.isChecked()
        flags['sum'] = self.ui.cfg_sum_chk.isChecked()
        flags['x0x1x2'] = self.ui.cfg_x0x1x2_chk.isChecked()
        flags['dop'] = self.ui.cfg_dop_chk.isChecked()
        return flags
        
    def show(self) :
        self.ui.show()   
    
    def set_pos (self, posx,posy) :
        self.ui.move(posx,posy)


# ----------------------------------------------------------------------------
# class popup crop_box
#-----------------------------------------------------------------------------
class crop_popup(QWidget):
    
    valeurs_onclose = Signal(str,str)
    
    def __init__(self, parent=None):
        super().__init__(parent, Qt.Tool | Qt.FramelessWindowHint)
        self.setWindowModality(Qt.NonModal)

        loader = QUiLoader()
        ui_file_name=resource_path('crop_box.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        #self.ui.adjustSize()
        self.ui.setFocusPolicy(Qt.StrongFocus)
        self.ui.setWindowFlags(Qt.Tool | Qt.FramelessWindowHint)
        
        
        self.ui.ok_btn.clicked.connect(self.ok)
        self.ui.cancel_btn.clicked.connect(self.cancel)
        self.ui.H_crop_text.setValidator(QtGui.QIntValidator(0, 999999))
        self.ui.L_crop_text.setValidator(QtGui.QIntValidator(0, 999999))

    def ok(self):
        self.valeurs_onclose.emit(self.ui.H_crop_text.text(), self.ui.L_crop_text.text())
        self.ui.close()
    
    def cancel(self):
        self.ui.H_crop_text.setText('0')
        self.ui.L_crop_text.setText('0')
        self.valeurs_onclose.emit(self.ui.H_crop_text.text(), self.ui.L_crop_text.text())
        self.ui.close()
        
# ----------------------------------------------------------------------------
# class popup shg_box
#-----------------------------------------------------------------------------
class shg_popup(QWidget):
    #  return focale camera et angle reseau
    valeurs_onclose = Signal(str, str)
    
    def __init__(self, parent=None):
        super().__init__(parent, Qt.Tool | Qt.FramelessWindowHint)
        self.setWindowModality(Qt.NonModal)

        loader = QUiLoader()
        ui_file_name=resource_path('shg_box.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        #self.ui.adjustSize()
        self.ui.setFocusPolicy(Qt.StrongFocus)
        self.ui.setWindowFlags(Qt.Tool | Qt.FramelessWindowHint)
        
        
        self.ui.ok_btn.clicked.connect(self.ok)
        self.ui.cancel_btn.clicked.connect(self.cancel)
        
        # validatators
        self.ui.shg_focale_camera_text.setValidator(QtGui.QIntValidator(0, 999999))
        self.ui.shg_angle_reseau_text.setValidator(QtGui.QIntValidator(0, 999999))
        
        # read settings
        if settings.value("Shg/foc") is not None :
            self.shg_foc_camera_other=settings.value("Shg/foc")
        else :
            self.shg_foc_camera_other='0'
        
        if settings.value("Shg/angle") is not None :
            self.shg_angle_other=settings.value("Shg/angle")
        else :
            self.shg_angle_other='0'
            
        self.ui.shg_focale_camera_text.setText( self.shg_foc_camera_other)
        self.ui.shg_angle_reseau_text.setText(self.shg_angle_other)

    def ok(self):
        # save settings
        settings.setValue("Shg/foc", self.ui.shg_focale_camera_text.text())
        settings.setValue("Shg/angle", self.ui.shg_angle_reseau_text.text())
        
        self.valeurs_onclose.emit( self.ui.shg_focale_camera_text.text(),self.ui.shg_angle_reseau_text.text())
        self.ui.close()
    
    def cancel(self):
       
        self.ui.shg_focale_camera_text.setText('0')
        self.ui.shg_angle_reseau_text.setText('0')
        self.valeurs_onclose.emit( self.ui.shg_focale_camera_text.text(),self.ui.shg_angle_reseau_text.text())
        self.ui.close()

# ----------------------------------------------------------------------------
# class galerie images à la fin d'un batch
#-----------------------------------------------------------------------------
class galerie_wnd(QDialog) :

    
    def __init__(self, img_list, name_list):
        super(galerie_wnd, self).__init__()
       
        
        #fichier GUI par Qt Designer
        loader = QUiLoader()
        ui_file_name=resource_path('galerie.ui')
        ui_file = QFile(ui_file_name)
        
        if not ui_file.open(QIODevice.ReadOnly):
            print(f"Cannot open {ui_file_name}: {ui_file.errorString()}")
            sys.exit(-1)
        
        self.ui = loader.load(ui_file)
        
        ui_file.close()
        
        self.img_list = img_list
        self.name_list=name_list
        
        self.read_settings()
        self.display_img(img_list, name_list)
        self.ui.img_list_view.itemClicked.connect(self.img_click)
        
    def display_img (self, img_list, name_list) :
        self.ui.img_list_view.clear()
        galx=300
        galy=300
        i=0
        
        for d in img_list :
            pro_item=QListWidgetItem()
     
            # transforme en 8 bits
            img_8bits= (d/256).astype(np.uint8)
            h,w =img_8bits.shape
            img_8bits = np.ascontiguousarray(img_8bits)
            image=QtGui.QImage(img_8bits.data, w, h ,w ,QtGui.QImage.Format_Grayscale8).copy()
            pix=QtGui.QPixmap.fromImage(image)
            self.pixmap = pix.scaled(galx, galy, Qt.AspectRatioMode.KeepAspectRatio,
                                     Qt.TransformationMode.SmoothTransformation)    
            pro_icon=QtGui.QIcon()
            pro_icon.addPixmap(self.pixmap)
            pro_item.setIcon(pro_icon)
            pro_item.setText(name_list[i])
            self.ui.img_list_view.addItem(pro_item)
            i=i+1
        
        screen_geom = QtGui.QGuiApplication.primaryScreen().availableGeometry()
        #dpr = QtGui.QGuiApplication.primaryScreen().devicePixelRatio()
        self.myscreen_w = int(screen_geom.right())
        h = galy
        w = len(img_list) * (galx+25)
        nh= int( w / self.myscreen_w )
        nw = int(self.myscreen_w/(galx+25))
        nw=min(nw,len(img_list))
        w= nw * (galx+25)+16
        h=(galy+50)*(nh+1)
        self.ui.resize(w,h)
        self.myimg_solo = img_wnd(None, False, self)
    
    def img_click (self,item) :
        index= self.ui.img_list_view.row(item)
        self.myimg_solo.show()
        #self.myimg_solo.on_ferme.connect(self.img_allclose)
        img_proc = np.fliplr(np.rot90(self.img_list[index], 3))
        self.myimg_solo.ui.inti_view.setImage(img_proc, autoRange=False)
        self.myimg_solo.set_title(self.name_list[index])
    
    
    def closeEvent (self,event):
        self.write_settings()
        #self.ui.on_ferme.emit()
        self.ui.close()
    
    def read_settings(self):
        settings=QSettings("Desnoux Buil", "inti_qt")
        self.ui.restoreGeometry(settings.value("GalWindow/geometry"))

    
    def write_settings(self) :
        # sauve settings
        settings=QSettings("Desnoux Buil", "inti_qt")
        settings.setValue("GalWindow/geometry", self.ui.saveGeometry())
    
    
# ----------------------------------------------------------------------------
# class redirection console to textEdit
#-----------------------------------------------------------------------------
class StdoutRedirector(QObject):
    new_text = Signal(str)

    def write(self, text):
        if text==' ' :
            self.new_text.emit(text)
        else :
            if text.strip() :
                self.new_text.emit(text)
        
        
    def flush(self) :
        pass

# ----------------------------------------------------------------------------
# class worker de thread
#-----------------------------------------------------------------------------
class Worker(QObject) :
    result_ready = Signal(object)
    #finished = Signal()
    
    """
    def __init__(self):
        super().__init__()
        self._running = True

    
    def stop(self):
        # Demande l'arrêt du traitement
        self._running = False
    """
    
    def run (self, *args, **kwargs) :
        #kwargs["stop_flag"] = lambda: not self._running
        result =sol.solex_proc (*args, **kwargs)
        self.result_ready.emit(result)
        #self.finished.emit()
    
    

# ----------------------------------------------------------------------------
# class validateur de nombre decimaux
#-----------------------------------------------------------------------------
class DotOnlyDoubleValidator(QtGui.QDoubleValidator):
    def validate(self, input_str, pos):
        # remplace la virgule par un point à la volée
        if ',' in input_str:
            input_str = input_str.replace(',', '.')
        return super().validate(input_str, pos)

    def fixup(self, input_str):
        # idem pour la correction finale
        return input_str.replace(',', '.')




#-----------------------------------------------------------------------------
# Utilitaires
#-----------------------------------------------------------------------------  

def to_float(s: str) -> float:
    return float(s.replace(',', '.'))

def file_exist(name):    
   
    if not os.path.exists(name):

        #print('ERROR: File ' + name + ' not found.')
        #print('End.') 
        return  False
    else:
        return True

def conv (mystr):
    try :
        a= to_float(mystr)
    except:
        mystr='0'
    return mystr

def get_baseline(f) :
    # f est le nom complet du fichier sans l'extension
    
    # f_short est le nom sans le chemin
    #f_short=os.path.split(f)[1]
    # f_path est le chemin
    #f_path=os.path.split(f)[0]
    #f=f_path+os.sep+f_short
    
    baseline = f
    
    suffixes = ["_disk", "_protus","_clahe", "_dp*_cont","_dp*", "_cont",
               "_doppler*", "_free", "_raw", "_color*", "_mix", "_inv",
               "_recon","_dp*recon"]
    
    for suffixe in suffixes :
        # On transforme le motif wildcard en regex : "*" → ".*"
        motif = re.escape(suffixe).replace("\\*", ".*") + r"$"
        m = re.match(f"^(.*){motif}", f)
        if m:
            baseline =m.group(1)
            #print("baseline "+str(baseline))
            return baseline
            
    
    return baseline

def get_data_ser (serfile) :
    try:
        scan = Serfile(serfile, False)
        dateSerUTC = scan.getHeader()['DateTimeUTC']
        f_dateSerUTC=datetime.fromtimestamp(SER_time_seconds(dateSerUTC), UTC)
        fits_dateobs=f_dateSerUTC.strftime('%Y-%m-%dT%H:%M:%S.%f7%z')[:19]
        FrameCount = scan.getLength()    #      return number of frame in SER file.
        Width = int(scan.getWidth())     #      return width of a frame, int est uint64
        Height = int(scan.getHeight())   #      return height of a frame, int car est uint64
        scan_size =(Height, Width)
    except:
        print(('Erreur ouverture fichier : ')+serfile)
    
    return fits_dateobs,FrameCount, scan_size

# from J.Meeus
def angle_P_B0 (date_utc):
    time = astropy.time.Time(date_utc)
    myJD=time.jd
    #myJD= 2460711.9700
    #date_1853JD2=2398167.2763889 # ref 9 nov 1853 18:37 
    #theta = ((myJD - date_1853JD2) /27.2743) +1
    #a=360*(theta-int(theta))
    #L0=360-a
    theta0 = (myJD - 2398220) * 360/25.38
    Rot_Carrington= (myJD-2398140.2270)/27.2752316 #JMeus 2nd edition p191
    


    I = 7.25
    K = 73.6667 + 1.3958333*(myJD - 2396758)/36525
    T = (myJD - 2451545)/36525
    Lo = (0.0003032*T + 36000.76983)*T + 280.46645
    M = ((-0.00000048*T - 0.0001559)*T + 35999.05030)*T + 357.52910
    C = ((-0.000014*T - 0.004817)*T + 1.914600)*math.sin(math.radians(M))
    #C = C +(-0.000101*T - 0.019993)*math.sin(math.radians(2*M)) + 0.000290*math.sin(math.radians(3*M))
    C = C +(-0.000101*T + 0.019993)*math.sin(math.radians(2*M)) + 0.000290*math.sin(math.radians(3*M)) #p164
    S_true_long = Lo + C
    Lambda = S_true_long - 0.00569 - 0.00478*math.sin(math.radians(125.04 - 1934.136*T))
    Lambda_cor = Lambda + 0.004419
    x = math.degrees(math.atan(-math.cos(math.radians(Lambda_cor)) * math.tan(math.radians(23.440144))))
    y = math.degrees(math.atan(-math.cos(math.radians(Lambda - K)) * math.tan(math.radians(I))))
    P = x + y
    Bo = math.degrees(math.asin(math.sin(math.radians(Lambda - K)) * math.sin(math.radians(I))))
    eta= math.degrees(math.atan(math.tan(math.radians(Lambda - K)) * math.cos(math.radians(I))))
    a = theta0 /360
    theta=(a-int(a))*360
    L0= eta - theta
    if L0 <0 :
        L0=L0+360
    
    return(str(round(P,2)),str(round(Bo,2)), str(round(L0,2)), str(int(Rot_Carrington)))

def img_force_color (img_doppler) :
    # BGR → HSV
    hsv = cv2.cvtColor(img_doppler, cv2.COLOR_RGB2HSV).astype(np.float32)
    H, S, V = cv2.split(hsv)
    
    # Correction orange → plus rouge
    mask_orange = (H > 5) & (H < 25)   # plage d'orange en degrés OpenCV (0-179) was 25
    H[mask_orange] -= 20               # décale la teinte vers le rouge
    S[mask_orange] *= 1.5              # booste la saturation
                       
    
    # Correction bleu → plus bleu
    mask_blue = (H > 90) & (H < 150)    # plage de bleu was 90
    H[mask_blue] += 20 
    S[mask_blue] *= 1.5 
    #V[mask_blue] *= 1.1
    #V[mask_blue] += 20
    
    # Après ton traitement en float32
    H = np.clip(H, 0, 179)
    S = np.clip(S, 0, 255) 
    V = np.clip(V, 0, 255) 
    
    # Reconstruction
    hsv_mod = cv2.merge([H, S, V]).astype(np.uint8)
    img_doppler = cv2.cvtColor(hsv_mod, cv2.COLOR_HSV2RGB)
    
    return img_doppler

def img_getdiam (img) :
    # detection des rayons
    axis=0 
    offset=0
    flag_disk=True
    ih,iw=img.shape
    m= ih//2
    img_c = img[m-ih//10:m+ih//10,50:-50]
    img_c = cv2.GaussianBlur(img_c,(101,101), sigmaX=50)
    a1,a2 = detect_bord (img_c, axis, offset, flag_disk)
    diam = a2-a1
    return diam

def img_crop_diam (img, diam, marge) :
    # image du disque est centré
    # si image disque etait un disque partiel il aura au préalable été paddé
    
    ih, iw = img.shape
    c=ih//2
    half = (diam + marge) // 2
    
    img_crop = img[c-half : c+half, c-half : c+half]
    
    return img_crop


def gong_orientation_auto(img1, img2, diam,ang_P) :
    # img1 image jpg de gong en 8 bits
    # img2 image _disk de inti en 16 bits
    
    debug = False
    
    # on seuil image gong
    sb=80
    img1[img1<sb]=0
    img1= (((img1-sb)/(255-sb))*255).astype(np.uint8)
    
    # on passe image en 8 bits
    img2=(img2 / 256).astype(np.uint8)
    
    gh,gw=img1.shape
    ih, iw = img2.shape
    
    if debug :
        plt.imshow(img1, cmap="grey")
        plt.show()
        plt.imshow(img2, cmap="grey")
        plt.show()

    # genere les 4 images possibles avec rotation angle P

    #rayon=(b2-b1)//2
    rayon=diam//2
    cote=int(rayon*(2**0.5))
    mi_cote=(cote//2)
    c= ih//2
    cote= mi_cote*2
    
    # carré inscrit est c-mi_cote
    img2=img2[c-mi_cote:c+mi_cote,c-mi_cote:c+mi_cote]
    img_g=np.array(img1[c-mi_cote:c+mi_cote,c-mi_cote:c+mi_cote], dtype='uint8')
    
    # Ajustement d'histogramme - transforme img_gong en fonction de l'histogramme de img2
    img_g = exposure.match_histograms(img_g, img2)
    img_g = np.clip(img_g, 0, 255).astype(np.uint8)
    # image Gong a beaucoup de détails
    #img_g=cv2.GaussianBlur(img_g, (11,11), 5)
    
    
    # test augmentation de contrast avec clahe
    # augmentation contrast sur image gong et disk
    clahe = cv2.createCLAHE(clipLimit=1, tileGridSize=(2,2)) # cliplimit was 0.8
    img_g = clahe.apply(img_g)
    img2=clahe.apply(img2)
    
    
    
    r_lum=np.mean(img2)/np.mean(img_g)
    #print("r_lum : ",r_lum)

    if r_lum<=1 : 
        img2=np.array(img2, dtype='uint8')    
        img_g=np.array(img_g*r_lum, dtype='uint8')
    else :
        img2=np.array(img2/r_lum, dtype='uint8')    
        img_g=np.array(img_g, dtype='uint8')

    
    
    if debug :
        plt.imshow(img_g, cmap="grey")
        plt.show()
        plt.imshow(img2, cmap="grey")
        plt.show()


    inv=['None','EW','NS','NS-EW']
    
    """
    # zones noires dues à la rotation perturbe la correlation
    
    # calcul de la matrice de rotation, angle en degre
    rotation_mat=cv2.getRotationMatrix2D((cote//2,cote//2),float(-ang_P),1.0)
                
    # application de la matrice de rotation
    #fr_rot=cv2.warpAffine(fr_avant_rot,rotation_mat,(w,h),flags=cv2.INTER_LINEAR)
    #frame=np.array(fr_rot, dtype='uint16')
    
    # calcul score SSD
    # Transformations possibles
    transformations = {
        "None": cv2.warpAffine(img2,rotation_mat,(cote,cote),flags=cv2.INTER_LINEAR),
        "NS": cv2.warpAffine(np.flipud(img2),rotation_mat,(cote,cote),flags=cv2.INTER_LINEAR),
        "EW":  cv2.warpAffine(np.fliplr(img2),rotation_mat,(cote,cote),flags=cv2.INTER_LINEAR),
        "NS-EW": cv2.warpAffine(np.flipud(np.fliplr(img2)),rotation_mat,(cote,cote),flags=cv2.INTER_LINEAR),
    }
    
    # Calcul du score pour chaque transformation
    scores = {}
    for nom, img in transformations.items():
        diff = img.astype(np.int16) - img_g.astype(np.int16)
        scores[nom] = np.sum(diff * diff)  # SSD

    # Trouver la transformation avec le score minimal
    myinv = min(scores, key=scores.get)
    #print(scores)
    #print(myinv)
   
    scores_psnr = {}
    scores_corr ={}
    for nom, img in transformations.items():
        scores_psnr[nom] = cv2.PSNR(img_g, img) 
        scores_corr[nom] = cv2.matchTemplate(img_g, img, cv2.TM_CCORR_NORMED)
        
    myinv_psnr = max(scores_psnr, key=scores_psnr.get)
    myinv_corr = max(scores_corr, key=scores_corr.get)
    
    if myinv_psnr != myinv_corr :
        print("Low confidence : psnr "+myinv_psnr +' corr '+myinv_corr)     
        
    """
    
    # calcul de deux scores
    score_PSNR =[cv2.PSNR(img_g, img2),cv2.PSNR(img_g, np.fliplr(img2)),
                cv2.PSNR(img_g, np.flipud(img2)),cv2.PSNR(img_g, np.flipud(np.fliplr(img2)))]
    score_corr =[cv2.matchTemplate(img_g, img2, cv2.TM_CCORR_NORMED),cv2.matchTemplate(img_g, np.fliplr(img2), cv2.TM_CCORR_NORMED),
                 cv2.matchTemplate(img_g, np.flipud(img2), cv2.TM_CCORR_NORMED),cv2.matchTemplate(img_g, np.flipud(np.fliplr(img2)), cv2.TM_CCORR_NORMED)]
    
    icorr=np.argmax(score_corr)
    ipsnr=np.argmax(score_PSNR)
    
    #print('Score PSNR : '+ str(score_PSNR))
    #print('Score correlation : ' + str(score_corr))
    #print(inv[ipsnr])
    myinv_psnr = inv[ipsnr]
    myinv_corr = inv[icorr]
    
    if icorr != ipsnr :
        print("Low confidence")
        print("Low confidence : psnr "+myinv_psnr +' corr '+myinv_corr)
    

    #return inv[icorr]
    return myinv_psnr



def gong (fits_dateobs, filename) :
    #fmt = '%Y%m%d'
    if fits_dateobs!='' :
        datemonth=fits_dateobs.split('T')[0].replace('-','')[:6]
        dateday=fits_dateobs.split('T')[0].replace('-','')
        r1="https://gong2.nso.edu/HA/hag/"+datemonth+"/"+dateday+"/"
        Vo_req=r1

        reponse_web=rq.get(Vo_req)
        sun_meudon=reponse_web.text.split("\n")
        t=sun_meudon[11].split('href=')[1].split(">")[0]
        t=t.replace('"','')
        web.open(r1+t)
        if os.path.exists(filename):
            web.open(filename)


def Colorise_Image (couleur, frame_contrasted, basefich, suff):

    # gestion couleur auto ou sur dropdown database compatibility
    # 'Manual','Ha','Ha2cb','Cah','Cah1v','Cak','Cak1v','HeID3'
    couleur_lbl=couleur
    if couleur_lbl == 'Manual' :
        couleur = 'on' # mode detection auto basé sur histogramme simple
    else :
        if couleur_lbl[:2] == 'Ha' :
            couleur='H-alpha'
        if couleur_lbl[:3] == 'Ha2' :
            couleur='Pale'
        if couleur_lbl[:2] == 'Ca' :
            couleur='Calcium'
        if couleur_lbl[:2] == 'He' :
            couleur='Pale'
    
    f=frame_contrasted/256
    f_8=f.astype('uint8')
    
    #hist = cv2.calcHist([f_8],[0],None,[256],[10,256])
    # separe les 2 pics fond et soleil
    th_otsu,img_binarized=cv2.threshold(f_8, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    hist = cv2.calcHist([f_8],[0],None,[256],[0,256])
    hist2=np.copy(hist)
    hist2[0:int(th_otsu)]=0
    pos_max=np.argmax(hist2)
    """
    p=np.percentile(hist2,90)
    print('p: ',p)
    h3=int(p)
    hist2[hist2<h3]=0
    hc=np.argwhere(hist2)
    e=hc[-1][0]-hc[0][0]
    print('e: ',hc[0][0],hc[-1][0],e)
    """
    debug_hist=False
    if debug_hist :
        plt.plot(hist2)
        plt.show()
        print('couleur : ',pos_max)

    
    # test ombres >> provoque des applats 
    ombres=False
    if ombres :
        
        i_low=[]
        i_hi=[]
        fr=np.copy(frame_contrasted)
        i_low=np.array((fr<(pos_max*256))*fr*1.01, dtype='uint16')
        i_hi=(fr>=pos_max)*fr
        fr=i_low+i_hi
        f=fr/256
        f_8=f.astype('uint8')
    
    
    if couleur =='on' :  
        if pos_max<200 and pos_max>=98 : #was 70
            couleur="H-alpha"
        if pos_max<98 :
            couleur="Calcium"
        if pos_max>=200 :
            couleur="Pale"

    
    # test ombres >> provoque des applats 
    ombres=False
    if ombres :
        f8_low=[]
        f8_hi=[]
        f8_low=np.array((f_8<pos_max)*f_8*1.05, dtype='uint8')
        f8_hi=(f_8>=pos_max)*f_8
        f_8=f8_low+f8_hi
    
    
    #couleur="H-alpha"
    
    if couleur != '' :
        # image couleur en h-alpha
        if couleur == 'H-alpha' :
            # build a lookup table mapping the pixel values [0, 255] to
            # their adjusted gamma values
            gamma=0.3   # was gam 1.3 > 0.3 ok un peu plus clair et 0.1 plus sombre sombre
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f1_gam= cv2.LUT(f_8, table)
            
            gamma=0.55 # was gam 0.5 - 0.3 trop rouge, 0.6 un peu jaune - 0.55 ok
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f2_gam= cv2.LUT(f_8, table)
            
            gamma=1 # gam is 1.0
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f3_gam= cv2.LUT(f_8, table)
            
            i1=(f1_gam*0.1).astype('uint8')     # was 0.05 - 1 trop pale - 0.1 ok
            i2=(f2_gam*1).astype('uint8')       # is 1
            i3=(f3_gam*1).astype('uint8')       # is 1
            
            gamma=1.5 # gam total image 2 est trop fade, 1.2 pas assez, 1.5 pas mal
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            i1= cv2.LUT(i1, table)
            i2= cv2.LUT(i2, table)
            i3= cv2.LUT(i3, table)
            
            img_color=np.zeros([frame_contrasted.shape[0], frame_contrasted.shape[1], 3],dtype='uint8')
            img_color[:,:,0] = np.array(i1, dtype='uint8') # blue
            img_color[:,:,1] = np.array(i2, dtype='uint8') # blue
            img_color[:,:,2] = np.array(i3, dtype='uint8') # blue
            
            # gestion gain alpha et luminosité beta
            #alpha=(255//2+10)/pos_max
            #print('alpha ', alpha)
            #img_color=cv2.convertScaleAbs(img_color, alpha=alpha, beta=0) # was 1.3 - 1.1 plus sombre - 1.2 ok
            
            # affiche dans clahe window for test
            #cv2.imshow('clahe',img_color)
            #cv2.setWindowTitle("clahe", "color")

            
        # image couleur en calcium
        if couleur == 'Calcium' :
            # build a lookup table mapping the pixel values [0, 255] to
            # their adjusted gamma values
            gamma=1.2  # was 1
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f1_gam= cv2.LUT(f_8, table)
            
            gamma=1 # was 0.8
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f2_gam= cv2.LUT(f_8, table)
            
            gamma=1 # was 0.8
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f3_gam= cv2.LUT(f_8, table)
            
            # i1: bleu, i2: vert, i3:rouge
            i1=(f1_gam*1).astype('uint8')     # was 0.05 - 1 trop pale - 0.1 ok
            i2=(f2_gam*0.7).astype('uint8')       # is 1
            i3=(f3_gam*0.7).astype('uint8')       # was 0.8 un peu trop violet
            
            gamma=1 # gam total image finalement aucun, 1.2 un peu fade
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            i1= cv2.LUT(i1, table)
            i2= cv2.LUT(i2, table)
            i3= cv2.LUT(i3, table)
            
            img_color=np.zeros([frame_contrasted.shape[0], frame_contrasted.shape[1], 3],dtype='uint8')
            img_color[:,:,0] = np.array(i1, dtype='uint8') # blue
            img_color[:,:,1] = np.array(i2, dtype='uint8') # green
            img_color[:,:,2] = np.array(i3, dtype='uint8') # red
            
            vp=np.percentile(f_8, 99.7)
            alpha=(255//2)/(vp*0.5)
            #print('alpha ', alpha)
            
            img_color=cv2.convertScaleAbs(img_color, alpha=alpha) # was 1.5 ok
            
            # affiche dans clahe window for test
            #cv2.imshow('clahe',img_color)
            #cv2.setWindowTitle("clahe", "color")
        
        # image couleur en jaune-orange (helium, sodium, continuum)
        if couleur == 'Pale' :
            # build a lookup table mapping the pixel values [0, 255] to
            # their adjusted gamma values
            gamma=1  # 
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f1_gam= cv2.LUT(f_8, table)
            
            gamma=1 # was 0.7
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f2_gam= cv2.LUT(f_8, table)
            
            gamma=1 # 
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            f3_gam= cv2.LUT(f_8, table)
            
            # i1: bleu, i2: vert, i3:rouge
            i1=(f1_gam*0.92).astype('uint8')     # was 0.5 
            i2=(f2_gam*0.98).astype('uint8')       # was 0.9
            i3=(f3_gam*1).astype('uint8')       # is 1
            
            gamma=0.5 # gam total image 1 trop fade, 0.7 pas mal
            invGamma = 1.0 / gamma
            table = np.array([((i / 255.0) ** invGamma) * 255
            for i in np.arange(0, 256)]).astype("uint8")
                # apply gamma correction using the lookup table
            i1= cv2.LUT(i1, table)
            i2= cv2.LUT(i2, table)
            i3= cv2.LUT(i3, table)
            
                
            img_color=np.zeros([frame_contrasted.shape[0], frame_contrasted.shape[1], 3],dtype='uint8')
            img_color[:,:,0] = np.array(i1, dtype='uint8') # blue
            img_color[:,:,1] = np.array(i2, dtype='uint8') # green
            img_color[:,:,2] = np.array(i3, dtype='uint8') # red
            
            #alpha=(255//2+50)/pos_max
            #alpha=1
            #print('alpha ', alpha)
            #img_color=cv2.convertScaleAbs(img_color, alpha=alpha) # was 1
            
            # affiche dans clahe window for test
            #cv2.imshow('clahe',img_color)
            #cv2.setWindowTitle("clahe", "color")
        
        #img_color=cv2.flip(img_color,0)
        
        #cv2.imshow('clahe',img_color)
        #cv2.setWindowTitle("clahe", "color")
            
        cv2.imwrite(basefich+suff+'_color_'+str(couleur)+'.png',img_color)
        
        
        print("Couleur : "+ str(couleur))

        return img_color     

def seuil_image_dyn (frameC, percent_h, dyn_bas):
    sub_frame=frameC[5:,:-5]
    Seuil_haut=np.percentile(sub_frame,percent_h)
    Seuil_bas=(Seuil_haut*dyn_bas)
    frameC[frameC>Seuil_haut]=Seuil_haut
    fcc=(frameC-Seuil_bas)* (65535/(Seuil_haut-Seuil_bas))
    fcc[fcc<0]=0
    #print('dyn : ' + str(Seuil_haut)+' '+str( Seuil_bas))
    frame_seuil=np.array(fcc, dtype='uint16')
    return frame_seuil

def seuil_image_percent (frameC, percent_h, percent_b, fmult_h):
    sub_frame=frameC[5:,:-5]
    Seuil_haut=np.percentile(sub_frame,percent_h)*fmult_h
    if Seuil_haut>65535 :
        Seuil_haut=65535
    Seuil_bas=np.percentile(sub_frame,percent_b)
    #frameC[frameC>Seuil_haut]=Seuil_haut
    fcc=(frameC-Seuil_bas)* (65535/(Seuil_haut-Seuil_bas))
    fcc[fcc<0]=0
    fcc[fcc>65535]=65535
    frame_seuil=np.array(fcc, dtype='uint16')
    return frame_seuil

def seuil_image (img):
    Seuil_haut=np.percentile(img,99.999)
    Seuil_bas=(Seuil_haut*0) # was 0.25 test dopppler protu
    img[img>Seuil_haut]=Seuil_haut
    img_seuil=(img-Seuil_bas)* (65535/(Seuil_haut-Seuil_bas))
    img_seuil[img_seuil<0]=0
    
    return img_seuil, Seuil_haut, Seuil_bas

def seuil_image_force (img, Seuil_haut, Seuil_bas):
    img[img>Seuil_haut]=Seuil_haut
    img_seuil=(img-Seuil_bas)* (65535/(Seuil_haut-Seuil_bas)) # was 65500
    img_seuil[img_seuil<0]=0
    
    return img_seuil

def get_lum_moyenne(img) :
    # ajout calcul intensité moyenne sur ROI centrée
    ih, iw =img.shape
    dim_roi = 200 # was 100 avant 6.7
    rox1 = iw//2 - dim_roi
    rox2 = iw//2 + dim_roi
    roy1 = ih//2 - dim_roi
    roy2 = ih//2 + dim_roi
    #print('roi ', rox1,rox2,roy1,roy2)
    try :
        lum_roi=np.mean(img[roy1:roy2,rox1:rox2])
    except:
        lum_roi=0
    return lum_roi

def disk_gauss (img, cx, cy, radius, flou) :
    # Créer un masque avec un disque blanc sur fond noir
    mask = np.zeros_like(img, dtype=np.uint16)
    # Paramètres du disque
    centre = (cx, cy)  # (x, y)
    # Dessiner un disque blanc
    cv2.circle(mask, centre, radius, 65535, -1)
    # Flouter le disque pour créer un dégradé doux
    mask_flou = cv2.GaussianBlur(mask, (0, 0), flou)
    # Mélange pondéré : image * (1 - alpha) + noir * alpha
    alpha = mask_flou.astype(np.float32) / 65535.0
    img_result = (img.astype(np.float32) * (1.0 - alpha)).astype(np.uint16)
    return img_result
    
    
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def data_path(relative_path):
    """ Get path to exe, works for dev and for PyInstaller """
    data_path = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), relative_path))
    return data_path


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# main app Qt
#-----------------------------------------------------------------------------  
#-----------------------------------------------------------------------------  

if __name__ == "__main__":
    
    # ajout d'_internal dans sys.path pour appel avec subprocess
    if getattr(sys, 'frozen', False):
        base_path = sys._MEIPASS if hasattr(sys, '_MEIPASS') else os.path.dirname(sys.executable)
    else:
        base_path = os.path.dirname(__file__)

    internal_path = os.path.join(base_path, "_internal")
    if os.path.isdir(internal_path) and internal_path not in sys.path:
        sys.path.insert(0, internal_path)
    
    # test du type d'execution pour trouver l'emplacement
    # des fichiers une fois compilé avec pyinstaller
   
    ui_file= True
    
    
    
    ## doit etre initialiser ici
    if ui_file==True :
        loader = QUiLoader()    

    # recherche de la langue dans les Qsettings
    try : # au cas ou je ne gere pas bien la valeur par defaut au premier lancement
        settings=QSettings("Desnoux Buil", "inti_qt")
        #LG = settings.value("App/lang",'Fr') # LG est 'Fr' (par defaut) ou 'En'
        LG = settings.value("App/lang")
        #print(LG)
    except :
        LG='FR'
        
    # pour eviter de devoir tuer app qt sous spyder
    app = QApplication.instance() 
    if not app:

        app = QApplication(sys.argv)
    else:
        app = QApplication.instance()
    
    app.setStyle('fusion') # pour forcer un beau look sur Mac
    
    #print(resource_path('lang_EN.qm'))
    
    fichier = sys.argv[1] if len(sys.argv) > 1 else None
    flag = sys.argv[2] if len(sys.argv) > 2 else None
    #print("sys arg"+str(fichier))
    
    if LG =='EN':
        translator=QTranslator(app)
        #translator.load("lang_EN")
        translator.load(resource_path('lang_EN.qm'))
        app.installTranslator(translator)
    
    my_wnd_class=main_wnd_UI(fichier, flag) 

    my_wnd_class.show()
    
    
    sys.exit(app.exec())
    #app.exec()
    
    # retour a la redirection des prints sur la console systeme
    sys.stdout=sys.__stdout__
