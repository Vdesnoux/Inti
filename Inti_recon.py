# -*- coding: utf-8 -*-
"""
Created on Thu Dec 31 11:42:32 2020

@author: valerie desnoux
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
import scipy.ndimage as ndimage
import os
#import time
from scipy.signal import savgol_filter
import cv2
import sys
import math
import config as cfg
from datetime import datetime
import time


from Inti_functions import *

from serfilesreader_vhd import Serfile

    


def solex_proc(serfile,Shift, Flags, ratio_fixe,ang_tilt, poly, data_entete,ang_P, solar_dict,param,image_queue):
    """
    ----------------------------------------------------------------------------
    Reconstuit l'image du disque a partir de l'image moyenne des trames et 
    des trames extraite du fichier ser
    avec un fit polynomial
    Corrige de mauvaises lignes et transversallium
 
    basefich: nom du fichier de base de la video sans extension, sans repertoire
    shift: ecart en pixel par rapport au centre de la raie pour explorer 
    longueur d'onde decalée
    ----------------------------------------------------------------------------
    """
    #plt.gray()              #palette de gris si utilise matplotlib pour visu debug
    
    # flag pour logger les temps d'execution
    debug_time = False # si True ne marche que sur un décalage à la fois, pas en doppler 
    t0=time.time()
    # flag antialisaing fort ou normal
    flag_aliasing_fort = True
    
    clearlog()
    #print (Shift)
    shift = Shift[0] # shift géré dans la constante du polynome

    shift_offdop=Shift[5]
    shift_dop1 = shift_offdop-Shift[1]
    shift_dop2 = shift_offdop+Shift[1]
    shift_cont = Shift[2]
    shift_offdop = Shift[5]

    shift_vol=Shift[3] # decacalage zeeman +/- 
    if Flags["WEAK"]:
        free_shift=Shift[3]
    if Flags["POL"]:
        shift_vol=Shift[3] # decacalage zeeman +/- 
        shift_zeeman=Shift[4] # decalage zeeman additionel
  
    flag_dopcont = Flags["DOPCONT"]
    flag_display = Flags["RTDISP"]
    sfit_onlyfinal = Flags["ALLFITS"]
    flag_pol = Flags['POL']
    flag_volume = Flags["VOL"]
    flag_weak = Flags["WEAK"]
    flag_flipRA = Flags["FLIPRA"]
    flag_flipNS = Flags["FLIPNS"]
    auto_crop = Flags["Autocrop"]
    flag_h20=Flags["h20percent"]
    # pour forcer les valeurs tilt et ratio
    flag_force = Flags["FORCE"]
    flag_corona = Flags['Couronne']
    flag_contonly = Flags['Contonly']
    deb1=int(param[4]) # en relatif par rapport à l'apparition du disque, première trame
    deb2=int(param[5]) # en relatif par rapport à l'apparition du disque, dernière trame

    
    # cas difficile fente 
    try :
        pos_fente_min=int(param[0])
        pos_fente_max=int(param[1])
    except :
        pos_fente_min = 0
        pos_fente_max = 0

    
    geom=[]

    # gere ici le nombre de longueur d'onde aka decalages à gérer
    # decalage nul pour l'image de la raie spectrale
    # decalage de shift_dop pour le doppler, decalage de shift_cont pour celui du continuum sur flag_dopcont
    # serie de décalage de 1 pixel entre les bornes +/- shift_vol sur flag_volume
    
    #img_suff=(['.fits','_d1_'+str(shift_dop)+'.fits', '_d2_'+str(shift_dop)+'.fits','_cont_'+str(shift_cont)+'.fits'])
    
    if flag_dopcont :
        if flag_contonly :
            range_dec=[0,shift_cont]
            logme ('shift cont : '+ str(shift_cont))
        else :
            range_dec=[0,shift_dop1,shift_dop2,shift_cont]
            logme ('shift cont : '+ str(shift_cont))
            logme ('shift doppler 1 : '+ str(shift_dop1))
            logme ('shift doppler 2 : '+ str(shift_dop2))
    else:
        range_dec=[0]
        """
        if shift==0:
            range_dec=[0]
        else:
            range_dec=[int(shift)]
        """
    
    if flag_volume :
        range_dec1=np.arange(-shift_vol,0)
        range_dec2=np.arange(1,shift_vol+1)
        range_dec=np.concatenate(([0],range_dec1,range_dec2), axis=0)
        
    
    if flag_pol :
        range_dec=[0,-shift_vol,shift_vol]
        
       
        
    if flag_weak :
        range_dec=[0, Shift[1], Shift[2]]
        flag_dopcont=True
        
        
    kend=len(range_dec)
        
    #print ('nombre de traitements ',kend)
    WorkDir=os.path.dirname(serfile)+"/"
    os.chdir(WorkDir)
    base=os.path.basename(serfile)
    basefich='_'+os.path.splitext(base)[0]
    
    # log timing execution pour analyse
    if debug_time :
        tim = open("timing"+basefich+".txt", "w", encoding="utf-8")
        tim.write(serfile+'\n')
    
    if shift != 0 :
        #add shift value in filename to not erase previous file
        dp_str= str(shift)
        dp_str="_dp"+str.replace(dp_str,".","_")
        basefich=basefich+dp_str # ajout '_' pour fichier en tete d'explorer
    
    # ouverture du fichier ser

    try:
        scan = Serfile(serfile, False)
    except:
        logme('Erreur ouverture fichier : '+serfile)
        
    FrameCount = scan.getLength()    #      return number of frame in SER file.
    Width = int(scan.getWidth())    #      return width of a frame, int est uint64
    Height = int(scan.getHeight())   #      return height of a frame, int car est uint64
    dateSerUTC = scan.getHeader()['DateTimeUTC']  
    dateSer=scan.getHeader()['DateTime']
    bitdepth=scan.getHeader()['PixelDepthPerPlane']
    #ser_header=scan._readExistingHeader()
    #print(ser_header)
    dtype=np.uint16
    if bitdepth==8:
        dtype = np.uint8
        if cfg.LG == 1:
            logme('Acquisition 8 bits.')
        else:
            logme('Acquisition 8 bits.') 
        #sys.exit()
    logme (serfile)
    
    if cfg.LG == 1:
        logme ('Largeur et hauteur des trames SER : ' + str(Width)+','+str(Height))
        logme ('Nombre de trames : '+str( FrameCount))
    else:
        logme ('Width and Height of SER frames : ' + str(Width)+','+str(Height))
        logme ('Number of frames : '+str( FrameCount))
 
    try:
        # Date UTC
        # Correction bug fromtimestamp convertie en prenant en compte heure locale PC
        # Fonction utc ne tient pas compte du reglage pc
        f_dateSerUTC=datetime.utcfromtimestamp(SER_time_seconds(dateSerUTC))
        logme('SER date UTC :' + f_dateSerUTC.strftime('"%Y-%m-%dT%H:%M:%S.%f7%z"'))
        fits_dateobs=f_dateSerUTC.strftime('%Y-%m-%dT%H:%M:%S.%f7%z')
        # date ser du PC
        f_dateSer=datetime.utcfromtimestamp(SER_time_seconds(dateSer))
        logme('SER date local :' + f_dateSer.strftime('"%Y-%m-%dT%H:%M:%S.%f7%z"'))
        

    except:
        if cfg.LG == 1:
            logme('Erreur lecture date-time fichier SER')
        else :
            logme('Error reading date-time SER file')
    
    #ok_flag=True              # Flag pour sortir de la boucle de lecture avec exit
    #FrameIndex=1              # Index de trame    

    # fichier ser avec spectre raies verticales ou horizontales (flag true)
    if Width>Height:
        flag_rotate=True
        cam_height=Width
    else:
        flag_rotate=False
        cam_height=Height
    
    # initialisation d'une entete fits (etait utilisé pour sauver les trames individuelles
    hdr= fits.Header()
    hdr['SIMPLE']='T'
    hdr['BITPIX']=32
    hdr['NAXIS']=2
    if flag_rotate:
        hdr['NAXIS1']=Height
        hdr['NAXIS2']=Width
    else:
        hdr['NAXIS1']=Width
        hdr['NAXIS2']=Height
    hdr['BZERO']=0
    hdr['BSCALE']=1
   
    #hdr['EXPTIME']=0
    try :
        hdr['DATE-OBS']=fits_dateobs
    except:
        fits_dateobs="2000-01-01T00:00:00"
        hdr['DATE-OBS']=fits_dateobs
        
    
    hdr['OBSERVER']=data_entete[0]
    hdr['INSTRUME']=data_entete[1]
    hdr['SITELONG']=data_entete[2] #Longitude of the telescope in decimal degrees, could East as positive
    hdr['SITELAT']=data_entete[3]  #Latitude of the telescope in decimal degrees 
    hdr['OBJNAME']='Sun'
    hdr['CONTACT']=data_entete[4]
    hdr['WAVELNTH']=data_entete[5]
    hdr['PHYSPARA']= 'Intensity'
    hdr['WAVEUNIT']= -10
    filename_suffixe=data_entete[1]+'_'+ data_entete[6]+'_'+fits_dateobs.split('T')[0].replace('-','')+'_'+fits_dateobs.split('T')[1].split('.')[0].replace(':','')
    hdr['CAMERA']  = data_entete [7] #String, Camera identification
    # Extract bin de data_entete data_entete[8]
    hdr['BIN1']=int(data_entete[8][0:1])
    hdr['BIN2']=int(data_entete[8][2:3])
    hdr['CAMPIX']  = to_float(data_entete [9])*1e-3 #float, pixel size of the CMOS or CCD sensor in mm
    hdr['CREATOR'] = data_entete [10] # String, name of the software which created the file and version number
    if data_entete[11]=='' :
        data_entete[11]=420
    if data_entete[12]=='':
        data_entete[12]=72
    hdr['FEQUIV']  = to_float(data_entete [11]) #float, equivalent focal length of the telescope in mm              
    hdr['APERTURE']= to_float(data_entete [12]) # float, aperture of the telescope in mm                 
    hdr['DIAPH'] = int(data_entete[21]) # integer, diametre après diaphragme
    hdr['DN'] = data_entete[22] # string, type de filtre genre ND8, Helioscope Lacerta
    hdr['SPECTRO'] = data_entete [13] #String, Spectrograph identification (SOLEX, SHG700...)
    hdr['FCOL']    = to_float(data_entete [14]) #float, focal length of the collimator in mm (80 SOLEX, 72 SHG700)
    hdr['FCAM']    = to_float(data_entete [15]) #float, focal length of the camera lens in mm (125 SOLEX, 72 SHG700)
    hdr['GROOVES'] = int(data_entete [16]) #integer, number of grooves/mm of the grating
    hdr['ORDER']   = int(data_entete [17]) #integer, interference order on the grating
    hdr['SHGANGLE']   = to_float(data_entete [18]) #float, angle between the central incident and diffracted rays in degrees
    hdr['SLWIDTH'] = round(to_float(data_entete [20])/1000, 3) #float, slit width in mm (0.010 or 0.007)
    hdr['SLHEIGHT'] = to_float(data_entete [19]) #float, slit height in mm (4.5 or 6)
    hdr['WAVEBAND'] = to_float(str(0.0)) #float, peut-etre obtenu avec param d'avant par calcul
    

    """
    ---------------------------------------------------------------------------
    #---- lecture fichier ser et accumulation pour moyenne
    ---------------------------------------------------------------------------
    """
    data_offset = 178
    frame_size = Width * Height
    block_size = 1000  # frames par bloc


    
    # Accumulateurs
    sum_all  = np.zeros((Height, Width), dtype=np.uint64)
    sum_good = np.zeros((Height, Width), dtype=np.uint64)
    
    count_all = 0
    count_good = 0
    
    mean_list=[]
    
    with open(serfile, "rb") as f:
        f.seek(data_offset)
        frame_index = 0
        
        while frame_index < FrameCount:
            n_to_read = min(block_size, FrameCount - frame_index)
            count = n_to_read * frame_size
            
            block = np.fromfile(f, dtype=dtype, count=count).reshape((n_to_read, Height, Width))
            
            if bitdepth == 8:
                block = block.astype(np.uint16) * 256
            else:
                block = block.astype(np.uint16)
            
            
            
            # moyennes des frames du bloc
            means = np.mean(block, axis=(1, 2))
            # On ajoute au tableau global
            mean_list.extend(means)  
            
            
            # accumulateur total
            sum_all += np.sum(block, axis=0, dtype=np.uint64)
            count_all += n_to_read
            
            # accumulateur bonnes trames
            good_mask = means > 3000
            if np.any(good_mask):
                sum_good += np.sum(block[good_mask], axis=0, dtype=np.uint64)
                count_good += np.count_nonzero(good_mask)
            
            frame_index += n_to_read
    
    # Décision finale
    if count_good > 500:
        mydata = sum_good
        kept_frame = count_good
    else:
        mydata = sum_all
        kept_frame = count_all
        
    if flag_rotate:
        # Rotation de 90° dans le sens horaire
        mydata = np.flip(mydata.swapaxes(0,1), axis=0)


    """
    ---------------------------------------------------------------------------
    #---- Calcul image moyenne
    ---------------------------------------------------------------------------
    """
    
    # calcul de l'image moyenne
    myimg=mydata/(kept_frame-1)             # Moyenne sur les kept frame
    myimg=np.array(myimg, dtype='uint16')   # Passe en entier 16 bits
    ih= int(hdr['NAXIS2'])                    # Hauteur de l'image
    iw= int(hdr['NAXIS1'])                    # Largeur de l'image
    myimg=np.reshape(myimg, (ih, iw))   # Forme tableau X,Y de l'image moyenne
    
        
    # sauve en fits l'image moyenne avec suffixe _mean
    savefich="Complements"+os.path.sep+basefich+'_mean'              
    SaveHdu=fits.PrimaryHDU(myimg,header=hdr)
    SaveHdu.writeto(savefich+'.fits',overwrite=True)
    
    if flag_weak :
        mytrame=np.zeros((hdr['NAXIS2'],hdr['NAXIS1']),dtype='uint64')
        if not flag_corona :
            # was -10,10
            deb1=-15
            deb2=5
        # calcul du début d'apparition du spectre, disparition et milieu
        mean_list=np.array(mean_list, dtype='uint16')
        mean_filtre = savgol_filter(mean_list,5,3)
        mean_gradient= np.gradient(mean_filtre)
        tram1=np.argmax(mean_gradient)
        tram2=np.argmin(mean_gradient)
        #print("Trame apparition : "+ str(tram1))
        #print("Trame disparition : "+ str(tram2))
        #print("Trame milieu : "+ str(tram1+(tram2-tram1)//2))
        
        
        with open(serfile, "rb") as f:
            frame1 = tram1+deb1
            frame2 = tram1+deb2
            
            if flag_corona :

                print("Selection : "+ str(frame1+deb1)+" , "+str(frame2+deb2))

            
            n_to_read = frame2 - frame1
            count = n_to_read * frame_size

            # déplacement du curseur en tenant compte de l'offset
            f.seek(data_offset + frame1 * frame_size * np.dtype(dtype).itemsize)
            frames = np.fromfile(f, dtype=dtype, count=count).reshape((n_to_read, Height, Width))
            
            # rotation éventuelle
            if flag_rotate:
                frames = np.flip(frames.swapaxes(1, 2), axis=1)
                
            
        for i in range (0, n_to_read) :          
            mytrame += frames[i]
        
        nbtrame = deb2- (deb1)
        # et on la sauve avec mode couronne ou pas
        mytrame = mytrame/(nbtrame)  
        mytrame = np.array(mytrame, dtype='uint16')
        mytrame = np.reshape(mytrame, (ih, iw))   # Forme tableau X,Y de l'image moyenne
        # mode couronne pour calcul polynome
        if flag_corona :
            savefich="Complements"+os.path.sep+basefich+'_mean'  
            myimg=np.copy(mytrame)
        else :
            savefich="Complements"+os.path.sep+basefich+'_mean_start'              
        SaveHdu=fits.PrimaryHDU(mytrame,header=hdr)
        SaveHdu.writeto(savefich+'.fits',overwrite=True)
    
    #gestion taille des images Ser et Disk
    # plus petit for speed up
    screensizeH = (864-50)*0.8
    screensizeW = (1536)*0.8
    
    # gere reduction image png
    nw = screensizeW/iw
    nh = screensizeH/ih
    sc=min(nw,nh)

    if sc >= 1 :
        sc = 1
    
    if image_queue == None : # modif macOS 
        #affiche image moyenne
        cv2.namedWindow('Ser', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('Ser', int(iw*sc), int(ih*sc))
        cv2.imshow ('Ser', myimg)
        cv2.moveWindow('Ser', 100, 0)
        
        if cv2.waitKey(2000) == 27:  # exit if Escape is hit otherwise wait 2 secondes
            cv2.destroyAllWindows()
            sys.exit()
        
        cv2.destroyAllWindows()
    else :
        image_queue.put(('mean', myimg,0))
        time.sleep(2)
        image_queue.put(('mean', None,0))
        

    t1=time.time()
    dt=t1-t0
    if debug_time : tim.write('fin lecture : '+"{:.2f}".format(dt)+'\n')
    
  
    """
    ----------------------------------------------------------------------------
    #---- Calcul polynome ecart sur l'image moyenne
    ----------------------------------------------------------------------------
    """
    # detect up and down limit of the spectrum of the mean image
    y1,y2=detect_bord(myimg, axis=1, offset=5, flag_disk=False)
    
    #if flag_corona :
       #y1= 100
       #y2=ih-100
        
    
    # calcul profile spectral en bonus
    pro = bin_to_spectre(myimg,y1,y2)
    lamb = np.arange(0,iw)
    pro_lamb=np.column_stack((lamb,pro))

    
    if cfg.LG == 1:
        logme('Image moyenne - Limites verticales y1,y2 : '+str(y1)+' '+str(y2))
    else:
        logme('Mean Image - Vertical limits y1,y2 : '+str(y1)+' '+str(y2))
    
    PosRaieHaut=y1
    PosRaieBas=y2
    
    # test deg3 ou deg2 poly
    deg=2
    
    if (flag_weak and Flags["FREE_AUTOPOLY"]!=1)   or (flag_pol and Flags["ZEE_AUTOPOLY"]!=1) : 
        # on force le polynome
        a=poly[0]
        b=poly[1]
        c=poly[2]
        LineRecal=1
        if Flags["WEAK"]==1 and free_shift !=0 :
            c=poly[2]+free_shift
    else :
        # on calcul le polynome à partir des min de la raie du haut jusqu'en bas
        c_offset=0 # offset de la constante polynome si l'image est croppée pour detection fente
        try :
            if pos_fente_min != 0 :
                myimg = myimg[:,pos_fente_min:]
                c_offset=pos_fente_min
                if cfg.LG == 1 :
                    logme ("zone detection fente xmin : "+str(pos_fente_min))
                else :
                    logme ("slit detection zone xmin : "+ str(pos_fente_min))
            if pos_fente_max != 0:
                myimg = myimg[:,0:-int(iw)+pos_fente_max]
                if cfg.LG == 1 :
                    logme ("zone detection fente xmax : "+ str(pos_fente_max))
                else :
                    logme ("slit detection zone xmax : "+ str(pos_fente_max))
        except:
            pass
        
        nogauss = True
        
        if nogauss :
            # THE calcul pour obtenir la position de la fente sur la raie la plus sombre
            MinX=np.argmin(myimg, axis=1)
            
            # on reduit à la zone du spectre
            marge=30
            MinX=MinX[PosRaieHaut+marge:PosRaieBas-marge]
            IndY=np.arange(PosRaieHaut+marge, PosRaieBas-marge,1)
            LineRecal=1
            #best fit d'un polynome degre 2, les lignes y sont les x et les colonnes x sont les y
        else :
            # test ajustement pour centre raie >> pas  de difference, abandon
            MinX_img=np.argmin(myimg, axis=1)
            # on reduit à la zone du spectre
            MinX=MinX_img[PosRaieHaut+30:PosRaieBas-30]
            posx=[]
            debug=False
            
            for j in range (PosRaieHaut+30, PosRaieBas-30) :
                posx.append (get_line_pos_absorption(myimg[j:j+1:][0],MinX_img[j],13,debug)[0])            
            #ecart=np.array((MinX-posx), dtype='float64')
            #np.savetxt("c:/codepy/Simu/polymin.txt", ecart,fmt='%.2f',delimiter=',',newline='\n')
            IndY=np.arange(PosRaieHaut+30, PosRaieBas-30,1)
            MinX=np.array(posx)
            LineRecal=1
        
        # ajustement robuste polynome
        if deg == 2 :
            sigma_clip=6.0 
            K=2
            
            mask = np.ones_like(MinX, dtype=bool)
            
            # detection outliers en K passes
            for _ in range(K):
                p=np.polyfit(IndY[mask],MinX[mask],2)
                #p=np.polyfit(IndY,MinX,2)
                fitted_x=np.polyval(p, IndY)
                residus = MinX-fitted_x
                std = np.std(residus[mask])
                mask = np.abs(residus) < sigma_clip * std
            
            # Calcul de l'erreur RMS
            #rms = np.sqrt(np.mean((MinX[mask] - fitted_x[mask])**2))
           # print(f"\nErreur RMS de l’ajustement polynome : {rms:.3f} pixels")

        else :
            # test rms -no diff, idem sous ISIS
            # on reste en deg2
            p=np.polyfit(IndY,MinX,3)
          
        
        #calcul des x colonnes pour les y lignes du polynome
        a=p[0]
        b=p[1]
        if deg == 2 :
            c=p[2] +c_offset
        else :
            c=p[2] # modif poly deg3
            d=p[3]+c_offset
        
        if Flags["WEAK"]==1 :
            if free_shift !=0 :
                c=p[2]+c_offset+free_shift
                if cfg.LG == 1 :
                    logme("Inclut decalage de : " + str(free_shift))
                else :
                    logme("Include shift of : " + str(free_shift))
            else:
                c=p[2]+c_offset
            
            
        
        if not Flags["WEAK"] and not Flags["POL"]: 
            if shift !=0:
                c= p[2]+c_offset+shift
                if cfg.LG == 1 :
                    logme("Inclut decalage de : " + str(shift))
                else :
                    logme("Include shift of : " + str(shift))
            
        
        # ajout du 11 Sept 2022
        poly=np.copy(p)
    
    
    if flag_pol :
        
        a=poly[0]
        b=poly[1]
        px=poly[2] # constante a y=0
        
        # Extraction zeeman
        # Calcul position X precise de la raie d'absortion d'où on extrait la polarisation
        # On extrait la ligne centrale dans myimg (au centre Y de l'image)
        posy = int(ih/2)
        line = myimg[posy:posy+1, 0: iw]
        # calcul constante à mi hauteur
        pos_line_mihauteur = a*posy**2 + b*posy + c
        print('pos_line_mihauteur '+ str(pos_line_mihauteur))
        # On calcul la position précise de la raie à mi hauteur
        posx = get_line_pos_absoption_old(line[0], pos_line_mihauteur, 13)
        print('pos_line_mihauteur ajustee '+ str(posx))
        # On calcule la constante actualisée du polynôme en incluant le decalage du a l'effet zeeman (zeeman wide)
        c = posx - a*posy**2 - b*posy + shift_zeeman
        #c = posx - a*posy**2 - b*posy 
        
        if shift_zeeman !=0 :
            if cfg.LG == 1 :
                logme("Inclut decalage de : " + str(shift_zeeman))
            else :
                logme("Include shift of : " + str(shift_zeeman))

              
    fit=[]
    #ecart=[]
    
    for y in range(0,ih):
        if deg==2 : 
            x=a*y**2+b*y+c 
        else : 
            x=a*y**3+b*y**2+c*y+d
        
        deci=x-int(x)
        # gere si modele polynome de la fente depasse la dimension image iw
        # sinon code vectorisation va generer un pb d'indice 
        if x<iw-3:
            fit.append([int(x)-LineRecal,deci,y])
            #fit.append([int(x),deci,y])
        else:
            fit.append([iw-3,0,y])
        #ecart.append([x-LineRecal,y])
    

        
    
    #if deg == 2 : 
    logme('Coef a*x2,b*x,c :'+"{:.4e}".format(a)+' '+str("{:.4e}".format(b))+' '+str("{:.2f}".format(c)))
    #else :
        #logme('Coef a*x3,b*x2,cx,d :'+"{:.4e}".format(a)+' '+str("{:.4e}".format(b))+' '+str("{:.2f}".format(c))+' '+str("{:.2f}".format(d)))
    
    np_fit=np.asarray(fit)
    xi, xdec,y = np_fit.T
    xdec=xi+xdec+LineRecal
    xi=xi+LineRecal
   
    posX_min_slit=int(min(xi[PosRaieHaut], xi[PosRaieBas]))
    posX_max_slit=int(max(xi[PosRaieHaut], xi[PosRaieBas]))
    
    debug_poly=False # attention plante si on recalcule pas le polynome (free not auto)
    
    if debug_poly :
        imgplot1 = plt.imshow(myimg)
        plt.scatter(xdec,y,s=0.1, marker='.', edgecolors=('red'))
        plt.show()
        #plt.scatter(vals,np.arange(ih),s=0.1, marker='.', edgecolors=('green'))
        
        #plt.scatter(MinX[mask],IndY[mask],s=0.1, marker='.', edgecolors=('blue'))
        plt.scatter(MinX,IndY,s=0.1, marker='.', edgecolors=('blue'))
        
        plt.scatter(xdec,y,s=0.1, marker='.', edgecolors=('red'))
        #plt.scatter(xi,y,s=0.1, marker='.', edgecolors=('green'))
        plt.show()
  
    
    # Test position de la fente décalée par rapport aux bords de l'image
    # modif interpolaton sur 4 points donc marge devient 3 et non plus 1
    if len(range_dec)!=1 :   # si il y des decalages à prendre en compte
        shift_max=range_dec[len(range_dec)-1]
        shift_min=range_dec[1]

        if ((posX_min_slit-shift_max)<=2) or ((posX_max_slit+shift_max)>=iw-3):
            v=min(posX_min_slit-2, iw-posX_max_slit-3)
            logme('*******************************')
            if cfg.LG == 1:
                logme ('Warning !! Valeur de decalage trop importante')
                logme ('Valeur de decalage maximum:'+str(v))
                logme ('Image au centre la raie sera calculee à la place')
            else:
                logme ('Warning!! Shift value too large')
                logme ('Maximum shift value :'+str(v))
                logme ('Image at line center will be computed instead')
            logme('*******************************')
            range_dec=[0]
            kend=len(range_dec)
            
    else:
        shift_max=abs(range_dec[0])
        if ((posX_min_slit-shift_max)<=2) or ((posX_max_slit+shift_max)>=iw-3):
            v=min(posX_min_slit-2, iw-posX_max_slit-3)
            logme('*******************************')
            if cfg.LG == 1:
                logme ('Warning !! Valeur de decalage trop importante')
                logme ('Valeur de decalage maximum:'+str(v))
                logme ('Image au centre la raie sera calculee à la place')
            else:
                logme ('Warning!! Shift value too large')
                logme ('Maximum shift value :'+str(v))
                logme ('Image at line center will be computed instead')
            logme('*******************************')
            range_dec=[0]
            kend=len(range_dec)
     

    t2=time.time()
    dt=t2-t1
    if debug_time : tim.write('calcul poly : '+"{:.2f}".format(dt)+'\n')

    """
    ----------------------------------------------------------------------------
    ----------------------------------------------------------------------------
    #---- Applique les ecarts a toute les lignes de chaque trame de la sequence
    ----------------------------------------------------------------------------
    ----------------------------------------------------------------------------
    """
    FrameIndex=0            # Index de trame

    FrameMax=FrameCount
    
    
    Disk=[]
    
    # TODO : verif reduc noise si on applique un shift
   
    # reduction de bruit - moyenne de 3 colonnes
    if len(range_dec)==1 and Flags["NOISEREDUC"] == 1 :
        range_dec_ref=range_dec
        range_dec=[0,-1,1]
    if Flags["WEAK"] and Flags["NOISEREDUC"] == 1 :
        range_dec_ref=range_dec
        range_dec=[0,-1,1,Shift[1]-1, Shift[1], Shift[1]+1,Shift[2]-1,Shift[2],Shift[2]+1]
    if Flags["POL"] and Flags["NOISEREDUC"] == 1 :
        range_dec_ref=range_dec
        range_dec=[0,-1,1,-shift_vol-1, -shift_vol, -shift_vol+1,shift_vol-1,shift_vol,shift_vol+1]

    
    for k in range_dec:
        Disk.append(np.zeros((ih,FrameMax), dtype='uint16'))
        
    
    
    if flag_display and image_queue == None :
   
        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
        cv2.moveWindow('image', 0, 0)
        cv2.resizeWindow('image', int(iw*sc), int(ih*sc))    

        cv2.namedWindow('disk', cv2.WINDOW_NORMAL)
        FrameMax=FrameCount
        cv2.resizeWindow('disk', int(FrameMax*sc), int(ih*sc))
        cv2.moveWindow('disk', int(iw*sc)+1, 0)
        #initialize le tableau qui va recevoir les intensités spectrale de chaque trame
        

    x_floors=[]
    t=[]
    
    for decalage in range_dec:
        s= int(decalage)
        decimale_decalage = decalage-s
        # indice de la partie entiere de la position avec eventuellement un decalage s
        x_floor = (np.asarray(fit)[:, 0] + np.ones(ih) * (LineRecal + s)).astype(int) # valeur entière de la position
        frac = np.asarray(fit)[:, 1] + decimale_decalage # valeur decimale de la position
        # il peut y a voir des valeurs de t qui sont supérieure à 1, il faut les ajouter à x_floor
        depassement = np.trunc(frac).astype(int)
        tt = frac - depassement
        x_floor = x_floor+depassement
        
        # teste des bornes pour indices left
        # antialiasing +1.5, cattmull +2, partie decimale max +1
        x_floor[(x_floor-4)<=0]=4
        x_floor[(x_floor+5)>=iw]=iw-5
              
        x_floors.append(x_floor)
        t.append(tt)
        
     
    
    # Lance la reconstruction du disk a partir des trames par block
    data_offset = 178
    frame_size = Width * Height
    block_size = 1000  # nombre de trames lues en une fois
    
    with open(serfile, "rb") as f:
        f.seek(data_offset)
        
        frame_index = 0
        while frame_index < FrameCount:
            # combien de frames à lire dans ce bloc
            n_to_read = min(block_size, FrameCount - frame_index)
            count = n_to_read * frame_size
    
            # lecture du bloc brut
            block = np.fromfile(f, dtype=dtype, count=count).reshape((n_to_read, Height, Width))
    
            # mise à l’échelle selon bitdepth
            if bitdepth == 8:
                block = block.astype(np.uint16) * 256
            else:
                block = block.astype(np.uint16)
    
            # rotation éventuelle
            if flag_rotate:
                block = np.flip(block.swapaxes(1, 2), axis=1)
    
            # --- ici tu fais tes calculs frame par frame ---
            for j in range(n_to_read):
                frame_abs_index = frame_index + j
                img = block[j]
    
                # si flag_display vrai montre trame en temps reel       
                if flag_display :
                    if image_queue == None :
                        cv2.imshow('image', img)
                        if cv2.waitKey(1)==27:
                            cv2.destroyAllWindows()
                            sys.exit()
                    else :
                        image_queue.put(("image", img, FrameCount))
                        
                # Boucle sur les decalages
                i=0
             
                # --- Précalcul des coefficients pour une fenêtre de 7 points et polyorder=2 ---
                # filtre polynomiale antialiasing sans reduction de la resolution spectrale
                coeffs = np.array([-2, 3, 6, 7, 6, 3, -2], dtype=np.float32) / 21.0
                offsets = np.array([-1.5,-1,-0.5,0,0.5, 1,1.5])
                    
                for i in range(0,len(range_dec)): 
                    
                    I_samples = np.zeros((ih, len(offsets)), dtype=np.float32)
                   
                    
                    for idx, dx in enumerate(offsets):
                        # Position fractionnaire
                        pos = x_floors[i] + dx +t[i]
                        pos_floor = np.floor(pos).astype(int)
                        frac = pos - pos_floor
                        try :
                            pt0=np.float64(img[np.arange(ih), pos_floor-1])
                            pt1=np.float64(img[np.arange(ih), pos_floor])
                            pt2=np.float64(img[np.arange(ih), pos_floor+1])
                            pt3=np.float64(img[np.arange(ih), pos_floor+2])
                            img4=[pt0,pt1,pt2,pt3]
                            I_samples[:, idx] = catmull_rom_vectorized(img4, frac)
                        except :
                            # a priori erreur catchée par bornes
                            I_samples[:,idx]= pt1
                        
                    
                    IntensiteRaie = I_samples @ coeffs
                    # IntensiteRaie = I_samples[:,3] si on enleve l'antialiasing
                    
                                
                    # converti en 16 bits
                    IntensiteRaie=np.clip(IntensiteRaie,0,65535).astype(np.uint16)
        
                    # Ajoute au tableau disk 
                    Disk[i][:,frame_abs_index]=IntensiteRaie
          
                
                # Display reconstruction of the disk refreshed every 30 lines
                refresh_lines=int(20)
                if flag_display and frame_abs_index %refresh_lines == 0 :
                    disk_display=[]
                    disk_display = np.copy(Disk[0]).astype(np.uint32)*2
                    disk_display=np.clip(disk_display,0,65535).astype(np.uint16)
                    if image_queue == None :
                        cv2.imshow ('disk', disk_display)
                        if cv2.waitKey(1) == 27:             # exit if Escape is hit
                                cv2.destroyAllWindows()
                                sys.exit()
                    else : 
                        image_queue.put(("disk",disk_display, FrameCount))
    
            # on passe au chunk suivant
            frame_index += n_to_read
    
           

    # on signale la fin de la transmission dans la queue
    if image_queue is not None :
        image_queue.put(('stop',None, 0))
    
    # reduction de bruit - moyenne de 3 colonnes
    if (kend == 1 and Flags["NOISEREDUC"] ==1)  or (Flags["NOISEREDUC"]== 1 and Flags["WEAK"]) or (Flags["NOISEREDUC"]== 1 and Flags["POL"]) :
        d0=np.array(Disk[0], dtype='float64')
        d1=np.array(Disk[1], dtype='float64')
        d2=np.array(Disk[2], dtype='float64')
        disk_reducnoise = (d0+d1*0.5+d2*0.5) / 2
        Disk[0]=np.array(disk_reducnoise, dtype='uint16')
        if Flags["WEAK"] or Flags["POL"] :
            d3=np.array(Disk[3], dtype='float64')
            d4=np.array(Disk[4], dtype='float64')
            d5=np.array(Disk[5], dtype='float64')
            disk_reducnoise = (d3*0.5+d4+d5*0.5) / 2
            Disk[1]=np.array(disk_reducnoise, dtype='uint16')
            d6=np.array(Disk[6], dtype='float64')
            d7=np.array(Disk[7], dtype='float64')
            d8=np.array(Disk[8], dtype='float64')
            disk_reducnoise = (d6*0.5+d7+d8*0.5) / 2
            Disk[2]=np.array(disk_reducnoise, dtype='uint16')
            del Disk[3:]


        range_dec = range_dec_ref
        kend=len(range_dec)
        if cfg.LG == 1 :
            logme("Option Reduction de Bruit")
        else :
            logme("Noise reduction option")

    # !!!! special serbannes fente photometrique
    # !!! special crop pour serbannes fente photometrique
    #hdr['NAXIS2']=hdr['NAXIS2'] -200
    #ih= int(hdr['NAXIS2'])                    
    #Disk[i] = Disk[i][0:ih,:]
    
    # Sauve fichier disque reconstruit pour affichage image raw en check final
    hdr['NAXIS1']=FrameCount-1
    img_suff=[]
    for i in range(0,kend):
        DiskHDU=fits.PrimaryHDU(Disk[i],header=hdr)
        if i>0 : 
            dp_str= str(range_dec[i])
            #dp_str='_dp_'+str.replace(dp_str,".","_")
            dp_str="_dp"+str.replace(dp_str,".","_")
            img_suff.append(dp_str)
            DiskHDU.writeto("Complements"+os.path.sep+basefich+img_suff[i]+'_raw.fits',overwrite='True')
        else:
            # le suffixe dp si shift non nul a deja été ajouté dans baseline
            img_suff.append('')
            
            DiskHDU.writeto("Complements"+os.path.sep+basefich+img_suff[i]+'_raw.fits',overwrite='True')

    if flag_display and image_queue == None:
        cv2.destroyAllWindows()
        
    t3=time.time()
    dt=t3-t2
    if debug_time : tim.write('image raw : '+"{:.2f}".format(dt)+'\n')
    
    frames=[]
    msg=[]
   
    
    
    if cfg.LG == 1:
        msg.append('...image centre...')
    else:
        msg.append('...center image...')

    if kend!=0 :
       
        for i in range(1,kend) :
            
            if cfg.LG == 1:
                msg.append('...décalage image '+str(range_dec[i])+"...")
            else:
                msg.append('...image shift '+str(range_dec[i])+"...")

    #msg=('...image centre...', '...image Doppler 1...', '...image Doppler 2...', '...image continuum...')
    
    
    """
    -------------------------------------------------------------------------
    -------------------------------------------------------------------------
    Boucle sur les disques de la liste de decalage range_dec
    -------------------------------------------------------------------------
    -------------------------------------------------------------------------
    """
    
    for k in range(0,kend):
        
        logme(' ')
        logme(msg[k])
        
        iw=Disk[k].shape[1]
        ih=Disk[k].shape[0]
        img=Disk[k]
        
        flag_noline = True
        
        if flag_noline == False :
            """
            --------------------------------------------------------------------
            #---- Calcul des mauvaises lignes et de la correction geometrique
            --------------------------------------------------------------------
            """       
        
            #print("bad lines")
            y1,y2=detect_bord (img, axis=1,offset=5, flag_disk=True)    # bords verticaux
            
            #detection de mauvaises lignes
            
            # somme de lignes projetées sur axe Y
            #img_blur=cv2.GaussianBlur(img,(3,3),0)
            ysum_img=np.mean(img,1)
            
            
            # ne considere que les lignes du disque avec marge de 15 lignes 
            marge=15
            ysum=ysum_img[y1+marge:y2-marge]
            
            # choix methode par fit polynome et division
            # ou detection des lignes et mediane 
            methode_poly= False
            debug=False
     
            # filtrage sur fenetre de 41 pixels, polynome ordre 3 (etait 101 avant)
            yc=savgol_filter(ysum,41, 3)
            """
            #polynome
            xval=np.arange(0,y2-y1-2*marge,1)
            p = np.polyfit(xval,ysum,6)
            fit=[]
            for x in xval:
                fitv = p[0]*x**6+p[1]*x**5+p[2]*x**4+p[3]*x**3+p[4]*x**2+p[5]*x+p[6]
                fit.append(fitv)
            #print('Coef poly ',p)
            #fp = ysum-np.array(fit)
            yc=np.array(fit)
            """
            if debug :
                # affichage debug
                plt.plot(yc)
                plt.plot(ysum)
                plt.show()
          
        
            # divise le profil somme par le profil filtré pour avoir les hautes frequences
            hcol1=np.divide(ysum,yc)
            hcol=np.copy(hcol1)
    
    
            # met à zero les pixels dont l'intensité est inferieur à 1.02 (2%)
            hcol[abs(hcol1-1)<= 0.02]=0
    
    
            if debug :
                # affichage debug
                plt.plot(hcol1)
                plt.plot(hcol)
                plt.show()
    
            # tableau de zero en debut et en fin pour completer le tableau du disque
            a=[0]*(y1+marge)
            b=[0]*(ih-y2+marge)
            hcol=np.concatenate((a,hcol,b))
            
            if debug :
                # affichage debug
                plt.plot(hcol1)
                plt.plot(hcol)
                plt.show()
            
    
            # creation du tableau d'indice des lignes a corriger
            l_col=np.where(hcol!=0)
            listcol=l_col[0]
            #print(listcol)
            
            # doit etre abandonné
            # demarre a x constant et gain trop fort
            # et poly ne convient pas pour disque tronqué
        
            if methode_poly:
                # on ne clamp pas le profil
                hcol=1+((ysum-yc)/yc)
                hcol[hcol-1>-0.03]=1
                # tableau de zero en debut et en fin pour completer le tableau du disque
                a=[1]*(y1+marge)
                b=[1]*(ih-y2+marge)
                hcol=np.concatenate((a,hcol,b))
                
                if debug :
                    #affichage debug
                    plt.plot(hcol)
                    plt.title('hcol_poly')
                    plt.show()
                
                # Geénère tableau image de flat 
                flat_hf=[]
                for i in range(0,iw):
                    flat_hf.append(hcol)
                    
                np_flat_hf=np.asarray(flat_hf)
                flat_hf = np_flat_hf.T
                
                # Evite les divisions par zeros...
                flat_hf[flat_hf==0]=1
                
                # somme des lignes mauvaises le long de x
                #hsum=np.mean(img[listcol],0)
                m=np.mean(img,0)
                hsum=np.max(m)/m
                hsum[hsum>2]=1
                #hsum=savgol_filter(hsum,31, 3)
                
                if debug :
                    plt.plot(hsum)
                    plt.show()
                
                # flat pondere en x 
                
                for i in listcol :
                    flat_hf[i,:]=flat_hf[i,:]/hsum
                   
                if debug :
                    # affichage debug
                    plt.imshow(flat_hf)
                    plt.show()
                    
                # Divise image par le flat
                b=np.divide(img,flat_hf)
                img=np.array(b, dtype='uint16')
                
            else :
                origimg=np.copy(img)
                # correction de lignes par filtrage median 13 lignes, empririque
                if debug : print(listcol)
                lines=[]
                lines_collection=[]
                try :
                    lines.append(listcol[0])
                    for i in range(1,len(listcol)-1):
                        if listcol[i]-listcol[i-1] ==1 :
                            lines.append(listcol[i])
                        else :
                            lines_collection.append(lines)
                            lines=[]
                            lines.append(listcol[i])
                    if debug : print(lines_collection)
                except:
                    pass
                
                for c in listcol:
                    c=c
                    m=origimg[c-11:c+10,] #now refer to original image
                    # calcul de la mediane sur 11-10 lignes de part et d'autres
                    # suppression de 2 colonnes en debut et fin qui peuvent être à zéro
                    s=np.median(m[:,2:-2],0)
    
                    #on prepare un patch de 2 valeurs 
                    a=[m[0][3],m[0][-3]]
                    #on ajoute les patchs
                    s=np.concatenate((a,s,a))
                    
                    #on remplace la ligne defectueuse
                    #s=65000 #pour visualiser la ligne idetifiée comme defectueuse
                    img[c:c+1,]=s #fix bug ecart une ligne
                
              
                """
                for c in lines_collection :
                    d=2
                    m1 = origimg[c[0]-d,]
                    m2 = origimg[c[len(c)-1]+d,]
                    m=np.array([m1,m2])
                    # suppression de 2 colonnes en debut et fin qui peuvent être à zéro
                    s=np.median(m[:,2:-2],0)
                    #on prepare un patch de 2 valeurs 
                    a=[s[0],s[-1]]
                    #on ajoute les patchs
                    s=np.concatenate((a,s,a))
                    
                    #on remplace la ligne defectueuse
                    #s=0
                    for i in c :
                        img[i:i+1,]=s #fix bug ecart une ligne
                """ 
                
            debug=False
            
            if debug :
                # affichage debug
                plt.plot(ysum_img)
                plt.plot(np.mean(img,1))
                plt.title('avant-apres')
                plt.show()
                # test sauve image apres correction pour debug
                DiskHDU=fits.PrimaryHDU(img,header=hdr)
                DiskHDU.writeto(basefich+img_suff[k]+'_line.fits',overwrite='True')
    
        
        """
        --------------------------------------------------------------
        #---- Correction de non uniformity - basse freq
        --------------------------------------------------------------
        """
        # en cours d'experimentation... non utilisé
        Flags['flat'] = False # False = disk
      
        frame=np.copy(img)      
        
        # on cherche la projection de la taille max du soleil en Y
        #print("non uniformity")
        y1,y2=detect_bord(frame, axis=1,offset=0, flag_disk=True) 
        #x1,x2=detect_bord(frame, axis=0,offset=0, flag_disk=True)
        if cfg.LG == 1:
            logme('Limites verticales y1,y2 : '+str(y1)+' '+str(y2))
        else:
            logme('Vertical limits y1,y2 : '+str(y1)+' '+str(y2))
        
    
        flag_nobords=detect_noXlimbs(frame)
            
        if Flags['flat'] :
            
            sb, seuil_haut = pic_histo (frame)
            myseuil=seuil_haut*0.3
            chull = (frame >= myseuil)
            frame_neg=np.where(chull, np.nan, frame)
            y_horsdisk=np.nanmedian(frame_neg,1)
            Smoothed_neg=savgol_filter(y_horsdisk,315, 3) 
            hf_neg=np.divide(y_horsdisk,Smoothed_neg) 
            flat = np.tile(hf_neg[:, np.newaxis], (1, iw))
            f_mask=flat*np.logical_not(chull)
            flat[f_mask==0]=1                
            
            debug = False
            if debug:
                plt.imshow(flat)
                plt.show()
    
            # Divise image par le flat
            BelleImage=np.divide(frame,flat)
            BelleImage[BelleImage>65535]=65535 # bug saturation
            frame=np.array(BelleImage, dtype='uint16')
            
        else :
            
            debug=False
            
            w1= ((y2-y1)//4)
            if w1%2 == 0 :
                w1 = w1 + 1 
            w2= int(0.3*w1)
            if w2%2 == 0 :
                w2 = w2 + 1 
            
            #print(w1, w2)
            
            for winterp in [w2,w1] :
                
                # si mauvaise detection des bords en x alors on doit prendre toute l'image
                if flag_nobords:
                    ydisk=np.median(img,1)
                    offset_y1=0
                    offset_y2=0
                else:
                    sb, seuil_haut = pic_histo (frame)
                    #print(seuil_haut*0.5)
               
                    myseuil=seuil_haut*0.5 # seuillage pour segmentation disque solaire 0.5 - 0.4 pour Ken H
                    #print("myseuil : ", myseuil)
                
                
                    if cfg.LowDyn : 
        
                        #guillaume
                        myseuil=seuil_haut*0.7 # faible dynamique, seuil pour segmenter le disque solaire was 0.6
                        #print(seuil_haut, myseuil)
                    
                    # filtre le profil moyen en Y en ne prenant que le disque
                    # pixel est dans disque si intensité supérieure à la moitié du seuil haut pic_histo
                    ydisk=np.empty(ih+1)
                    offset_y1=0
                    offset_y2=0
                    
                    #creation d'un mask
                    #chull=np.where((frame < (seuil_haut*0.7)) & (frame > seuil_haut) ,False,True)
                    chull = (frame >= seuil_haut*0.5)
                    
                    frameg=np.where(chull, frame, np.nan)

                    
                    for j in range(0,ih):
                        temp=np.copy(frameg[j,:])
                        temp=temp[temp>myseuil]
                        if len(temp)!=0:
                            ydisk[j]=np.nanmedian(temp)
                          
                        else:
                            # manage poor disk intensities inside disk
                            # avoid line artefact
                            if j>=y1 and j<=y2:
                                ydisk[j]=myseuil
                                if abs(j-y1) < abs(j-y2):
                                    offset_y1=offset_y1+1
                                else:
                                    offset_y2=offset_y2-1
        
                            else:
                                ydisk[j]=myseuil

            
                # ne prend que le profil des intensités pour eviter les rebonds de bords
                
                # manage flat range application
                y1=y1+offset_y1
                y2=y2+offset_y2
                
                ToSpline= ydisk[y1:y2]
                
                # traitement du cas ou le disque est saturé - cas de la couronne solaire raie verte
                moy_profil=np.mean(ToSpline)
                #print ('moy profil : ', moy_profil)
                if moy_profil <64000 :
                
                    #winterp=801 # en test, 101 ne corrige pas de nombreuses lignes en pack 
                    # test pour sunscan
                    #winterp=101 #version 6.7
                    if len(ToSpline)<winterp :
                        
                        if cfg.LG == 1:
                            logme('Hauteur du disque anormalement faible : '+str(y1)+' '+str(y2))
                        else:
                            logme('Disk Height abnormally low : '+str(y1)+' '+str(y2))
                        
                        if len(ToSpline)%2==0 :
                            winterp=len(ToSpline)-10+1
                        else:
                            winterp=len(ToSpline)-10
                    
                    Smoothed2=savgol_filter(ToSpline,winterp, 3) # window size, polynomial order
                    # Divise le profil reel par son filtre ce qui nous donne le flat
                    hf=np.divide(ToSpline,Smoothed2)             
                                                          
                    
                    
                    # Elimine possible artefact de bord lie au filtrage
                    mg = (winterp-1)//4
                    # mg was 5 6.7
                    hf=hf[mg:-mg]
                    #reconstruit le tableau du profil complet an completant le debut et fin
                    a=[1]*(y1+mg)
                    b=[1]*(ih-y2+mg)
                    hf=np.concatenate((a,hf,b))
                    
                    ToSpline=np.concatenate((a,ToSpline,b))
                    Smoothed2=np.concatenate((a,Smoothed2,b))
                    
                    if debug:
                        plt.plot(ToSpline)
                        plt.plot(Smoothed2)
                        plt.show()
                        plt.plot(hf)
                        plt.show()
                    
                    
                    # optimisation
                    flat = np.tile(hf[:, np.newaxis], (1, iw))
                    """
                    # Génère tableau image de flat 
                    flat=[]
                    for i in range(0,iw):
                        flat.append(hf)
                        
                    np_flat=np.asarray(flat)
                    flat = np_flat.T
                    """
                    # gere le flat sur le disque uniquement avec le mask chull
                    f_mask=flat*chull
                    
                    # Evite les divisions par zeros...
                    flat[f_mask==0]=1                
                    
                    if debug:
                        plt.imshow(flat)
                        plt.show()
            
                    # Divise image par le flat
                    BelleImage=np.divide(frame,flat)
                    BelleImage[BelleImage>65535]=65535 # bug saturation
                    frame=np.array(BelleImage, dtype='uint16')
                    
                                  
                else:
                    # pas de correction de flat on reprend l'image
                    logme ("Pas de correction de flat, profil saturé")
                    frame=np.copy(img)
                    break
                # repart pour boucle winterp
                    
            if debug:
                # sauvegarde de l'image deflattée
                #DiskHDU=fits.PrimaryHDU(frame,header=hdu.header)
                DiskHDU=fits.PrimaryHDU(frame,header=hdr)
                DiskHDU.writeto(basefich+img_suff[k]+'_flat.fits',overwrite='True')  
  
        
        # on sauvegarde les bords haut et bas pour les calculs doppler et cont
        if k==0: 
            y1_img=y1
            y2_img=y2
        
        t4=time.time()
        dt=t4-t3
        if debug_time : tim.write('flat : '+"{:.2f}".format(dt)+'\n')
       
        """
        -----------------------------------------------------------------------
        #---- Calcul du tilt 
        Si on voit les bords du soleil
        sinon on n'applique pas de correction de tilt,
        on applique un facteur SY/SX=0.5
        et on renvoit a ISIS
        -----------------------------------------------------------------------
        """
        
        img2=np.copy(frame)
        EllipseFit=[]
        crop=0
        
        if k != 0:
            flag_force = True
    
        #if float(ang_tilt)==0:
              
       
  
        if not(flag_nobords):
            
            # methode fit ellipse pour calcul de tilt
            # zone d'exclusion des points contours zexcl en pourcentage de la hauteur image
            #print("detect edge tilt")
            X = detect_edge (img2, zexcl=0.1, crop=crop, disp_log=False)
            EllipseFit,XE=fit_ellipse(img2, X,disp_log=False)
            
            # correction de tilt uniquement si on voit les limbes droit/gauche
            # trouve les coordonnées y des bords du disque dont on a les x1 et x2 
            # pour avoir les coordonnées y du grand axe horizontal
            # on cherche la projection de la taille max du soleil en Y et en X

           
            # methode calcul angle de tilt avec XE ellipse fit
            elX=XE.T[0]
            elY=XE.T[1]
            el_x1=np.min(elX)
            el_x2=np.max(elX)
            el_ind_x1= np.argmin(elX)
            el_ind_x2= np.argmax(elX)
            el_y_x1=elY[el_ind_x1]
            el_y_x2=elY[el_ind_x2]
            #print('ellipse x1,x2 : ', el_x1, el_x2)
            #print('ellipse y_x1,y_x2 : ', el_y_x1, el_y_x2)
            
            if k==0 :
                el_x1_img=el_x1
                el_x2_img=el_x2
                
            #if float(ang_tilt)==0:
            if flag_force != True :
                # calcul l'angle de tilt ellipse
                dy=(el_y_x2-el_y_x1)
                dx=(el_x2-el_x1)
                TanAlpha=(-dy/dx)
                AlphaRad=math.atan(TanAlpha)
                AlphaDeg=math.degrees(AlphaRad)
            
        
            else :
                AlphaDeg=float(ang_tilt)
                AlphaRad=math.radians(AlphaDeg)
                TanAlpha=math.tan(AlphaRad)
                
            
            if cfg.LG == 1:
                logme('Angle de Tilt : '+"{:+.4f}".format(AlphaDeg))
            else:
                logme('Tilt angle : '+"{:+.4f}".format(AlphaDeg))
            
            #on force l'angle de tilt pour les prochaines images
            ang_tilt=AlphaDeg

    
            #decale lignes images par rapport au centre
            colref=round(((el_x2-el_x1)/2)+el_x1) # bug 27 mai 2026 mauvais calcul du centre
            dymax=int(abs(TanAlpha)*(colref))
            
            # calcul du background
            # Taille de l'image
            #ih, iw = img2.shape
            
            # Pré-calcul du percentile bas si nécessaire
            if y1 <= 15 or y2 > ih - 15:
                perc15 = np.percentile(img2, 15)
            
            # Partie haute
            if y1 <= 15:
                img_dark1 = np.full((2, iw), perc15)
            else:
                img_dark1 = img2[0:10, :] 
            
            # Partie basse
            if y2 > ih - 15:
                img_dark2 = np.full((2, iw), perc15)
            else:
                img_dark2 = img2[-10:, :] 
            
            # Calculs des moyennes de remplissage (1D)
            #img_fill1 = img_dark1.mean(axis=0)
            #img_fill2 = img_dark2.mean(axis=0)
            
            # Fond basé sur l'ensemble des zones sombres
            # Astuce : np.vstack évite un appel à np.concatenate avec tuple
            img_dark = np.vstack((img_dark1, img_dark2))
            background = np.percentile(img_dark, 55) 
            
              
            # Ajout de padding en haut et en bas
            img2 = np.pad(img2, ((dymax, dymax), (0, 0)), mode='constant', constant_values=background)
            ih, iw = img2.shape
            crop=int((TanAlpha)*iw)
            
            # Grille de coordonnées
            y = np.arange(ih)[:, None]   # forme (ih, 1)
            x = np.arange(iw)[None, :]   # forme (1, iw)
            dy = (x - colref) * TanAlpha # décalage vertical par colonne
            
            # Coordonnées à échantillonner
            y_src = y - dy                     # (ih, iw)
            x_src = np.broadcast_to(x, y_src.shape)
            
            # Interpolation bilinéaire (ordre=1)
            NewImg = map_coordinates(img2, [y_src, x_src], order=1, cval=background)
            
            # Clamp à 0
            NewImg[NewImg <= 0] = 0
        
            img2=np.copy(NewImg)
            

                   
        debug_tilt=False
        
        if debug_tilt==True:
            # sauvegarde en fits de l'image tilt
            img2=np.array(img2, dtype='uint16')
            DiskHDU=fits.PrimaryHDU(img2,header=hdr)
            DiskHDU.writeto(basefich+img_suff[k]+'_tilt.fits', overwrite='True')
        
        
        t5=time.time()
        dt=t5-t4
        if debug_time : tim.write('tilt : '+"{:.2f}".format(dt)+'\n')
        
        """
        ----------------------------------------------------------------
        #--- Calcul du parametre de scaling SY/SX
        ----------------------------------------------------------------
        """
        
       
        
        if flag_nobords:
            ratio_fixe=0.5
            
        if k==1:
            ratio_fixe=ratio_fixe_d1
            
        #if ratio_fixe==0:
        if flag_force != True and not flag_nobords:
            # methode fit ellipse pour calcul du ratio SY/SX
            X = detect_edge (img2, zexcl=0.1,crop=crop, disp_log=False)
            EllipseFit,XE=fit_ellipse(img2, X,disp_log=False)
            
            ratio=EllipseFit[2]/EllipseFit[1]
       
            
            if cfg.LG == 1:
                logme('Facteur d\'echelle SY/SX : '+"{:+.4f}".format(ratio))
            else:
                logme('Scaling SY/SX : '+"{:+.4f}".format(ratio))
            
            NewImg, newiw = circularise2_opt(img2,iw,ih,ratio)
          
        
        else:
            # Forcer le ratio SY/SX
            NewImg, newiw=circularise2_opt(img2,iw,ih,ratio_fixe)

            if cfg.LG == 1:
                logme('Facteur d\'echelle fixe SY/SX : '+"{:+.4f}".format(ratio_fixe))
            else:
                logme('Fixed Scaling SY/SX : '+"{:+.4f}".format(ratio_fixe))

        
        if k==0:
            #if ratio_fixe==0:
            if flag_force != True and not flag_nobords:
                ratio_fixe_d1=ratio
            else:
                ratio_fixe_d1=ratio_fixe
        
        
        frame=np.array(NewImg, dtype='uint16')
        #print('shape',frame.shape)
        
        debug_SYSX=False
        
        if debug_SYSX==True:
            # sauvegarde en fits de l'image tilt
            img2=np.array(NewImg, dtype='uint16')
            DiskHDU=fits.PrimaryHDU(NewImg,header=hdr)
            DiskHDU.writeto(basefich+img_suff[k]+'_sysx.fits', overwrite='True')        
       
       
        """
        ----------------------------------------------------------------------
        #---- Sanity check, second iteration et calcul des parametres du disk occulteur
        ----------------------------------------------------------------------
        """
        if flag_nobords:
            cercle=[0,0,0,0]
            r=0
            
        else:
            # fit ellipse pour denier check
            # zone d'exclusion des points contours zexcl en pourcentage de la hauteur image 
    
            X = detect_edge (frame, zexcl=0.1, crop=crop, disp_log=False)
            EllipseFit,XE=fit_ellipse(frame, X, disp_log=False)
            
            ratio=EllipseFit[2]/EllipseFit[1]
            
            if abs(ratio-1)>=1 and k==0 and ratio_fixe==0: #0.01
                logme('Ratio iteration2 :'+ str(ratio))
                NewImg, newiw=circularise2_opt(frame,newiw,ih,ratio)
                frame=np.array(NewImg, dtype='uint16')
                X= detect_edge (frame, zexcl=0.1,crop=crop, disp_log=False)
                EllipseFit,XE=fit_ellipse(frame, X, disp_log=False)
                
            if k==0:
                xc=round(EllipseFit[0][0])
                yc=round(EllipseFit[0][1])
                wi=round(EllipseFit[1]) # diametre
                he=round(EllipseFit[2])
                #print('cercle '+str(EllipseFit[0][0])+' '+str(EllipseFit[0][1]))
                cercle=[xc,yc,wi,he]              
                
                #if ratio_fixe==0 :
                if flag_force != True :
                    cercle=[xc,yc,wi,he]

                r=round(min(wi-5,he-5)-4)
                if cfg.LG == 1:
                    logme('Final SY/SX :'+ "{:+.3f}".format(he/wi))
                    #logme('Centre xc,yc et rayon : '+str(xc)+' '+str(yc)+' '+str(int(r)))
                else:
                    logme('Final SY/SX :'+ "{:+.3f}".format(he/wi))
                    #logme('xc,yc center and radius : '+str(xc)+' '+str(yc)+' '+str(int(r)))
    
        t6=time.time()
        dt=t6-t5
        if debug_time : tim.write('SX/SY : '+"{:.2f}".format(dt)+'\n')
        
        if k==0:
            cercle0=np.copy(cercle)
            geom.append(ratio_fixe_d1)
            geom.append(ang_tilt)
            
        x0=cercle0[0]
        y0=cercle0[1]
        wi=round(cercle0[2])
        he=round(cercle0[3])
    
        
        sensNS_corr_ang_tilt=1
        sensEW_corr_ang_tilt=1
             
        # retourne horizontalement image si flag de flip est vrai
        if flag_flipRA:
             frame=np.flip(frame, axis=1)
             # il faut recentrer le centre du disque en X...
             x0=frame.shape[1]-x0
             cercle[0]=x0
             #if k==0 : cercle0[0]=x0
             sensEW_corr_ang_tilt=-1
             if cfg.LG == 1:
                 logme("Inversion EW")
             else:
                 logme("EW inversion")
         
        if flag_flipNS:
             frame=np.flip(frame, axis=0)
             # il faut recentrer le centre du disque en Y...
             y0=frame.shape[0]-y0-1
             cercle[1]=y0
             #if k==0: cercle0[1]=y0
             sensNS_corr_ang_tilt=-1
             if cfg.LG == 1:
                 logme("Inversion NS")
             else:
                 logme("NS inversion")
          
        """
        -----------------------------------------------------------------------
        #---- Correction rotation angle P et tilt
        -----------------------------------------------------------------------
        """
        # corrige angle P de rotation si non nul
        # on corrige aussi de l'angle de tilt dans geom[1]
        
        # recalcul angle de tilt apres circularisation
        # calcule angle de tilt apres rescaling
        
        if geom[1] != 0 :
            Arad=math.radians(geom[1])
            TanA=math.tan(Arad)
            TanA_scaled= TanA/geom[0]
            A2Rad=math.atan(TanA_scaled)
            Tilt_scaled_Deg=math.degrees(A2Rad)
        else :
            Tilt_scaled_Deg = 0
        

        angle_rot = ang_P+sensNS_corr_ang_tilt*sensEW_corr_ang_tilt*Tilt_scaled_Deg
        
        
        if angle_rot !=0 :
            
            
            # on copie l'image de travail et on recupere hauteur et largeur
            fr_avant_rot=np.copy(frame)
            hh,ww=fr_avant_rot.shape[:2]
            
            
            if cercle[2]*2 > cam_height :
                logme('Disque partiel ')
                ang_P = 0
                angle_rot = ang_P+sensNS_corr_ang_tilt*sensEW_corr_ang_tilt*Tilt_scaled_Deg
               
            if cfg.LG == 1:
                #logme("Angle de correction : "+ str(angle_rot))
                logme("Angle P : "+str(ang_P))
            else:
                #logme("Correction angle : "+ str(angle_rot))
                logme("P angle : "+str(ang_P))


                        
            # et on met a jour dimensions de l'image avant rotation
            h,w=fr_avant_rot.shape[:2]
            # calcul de la matrice de rotation, angle en degre
            rotation_mat=cv2.getRotationMatrix2D((int(cercle[0]),int(cercle[1])),float(angle_rot),1.0)
                        
            # application de la matrice de rotation
            fr_rot=cv2.warpAffine(fr_avant_rot,rotation_mat,(w,h),flags=cv2.INTER_LINEAR)
            frame=np.array(fr_rot, dtype='uint16')
            
            sb,sh = pic_histo(fr_avant_rot)
            frame[frame<=10]=sb
        
        else :
            h, w = frame.shape
            
        # modif du 11 mars 2025 sur remplacement de r qui est trop faible (-5)
        """
        if cfg.LG == 1:
            logme('Centre xc,yc et rayon : '+str(cercle[0])+' '+str(cercle[1])+' '+str(int(r)))
        else:
            logme('xc,yc center and radius : '+str(cercle[0])+' '+str(cercle[1])+' '+str(int(r)))
        """
        if cfg.LG == 1:
            logme('Centre xc,yc et rayon : '+str(cercle[0])+' '+str(cercle[1])+' '+str(cercle[2]))
        else:
            logme('xc,yc center and radius : '+str(cercle[0])+' '+str(cercle[1])+' '+str(cercle[2]))
            
        
        t7=time.time()
        dt=t7-t6
        if debug_time : tim.write('rotation : '+"{:.2f}".format(dt)+'\n')
        
        # on croppe et on centre
        # Hauteur du capteur est dans cam_Heigth=scan.height
        
        debug_crop=False
        #auto_crop=True
        
        if flag_nobords :
            auto_crop=False
            cercle=[0,0,0,0]
        
        if auto_crop :
            
            #cercleC, crop_he, crop_wi, crop_img=auto_crop_img(cam_height, h, w, frame, cercle0, debug_crop)
            cercleC, crop_he, crop_wi, crop_img=auto_crop_img(cam_height, h, w, frame, cercle, debug_crop, param, flag_h20)
            logme("Dim H x L : "+str(crop_he)+' x '+str(crop_wi))
            if crop_he==0 :
                auto_crop=False
                if cfg.LG == 1:
                    print("Erreur AutoCrop")
                else:
                    print ("Error AutoCrop")
        
        # doit changer les coord du centre
        
        if auto_crop :
            
            # Sauvegarde en fits de l'image finale
            frame=np.array(crop_img, dtype='uint16')
            hdr['NAXIS1']=crop_wi
            hdr['NAXIS2']=crop_he
            
            hdr['INTI_XC'] = cercleC[0]
            hdr['INTI_YC'] = cercleC[1]
            hdr['INTI_R'] = cercleC[2]
            # duplication keywords pour BASS2000
            hdr['CENTER_X'] = cercleC[0]
            hdr['CENTER_Y'] = cercleC[1]
            hdr['SOLAR_R'] = cercleC[2]
            
            # recalcul des coordonnées haut et bas du disque
            y1_img,y2_img = detect_bord(frame, axis=1, offset=0,flag_disk=True)
            x1_img,x2_img = detect_bord(frame, axis=0, offset=0,flag_disk=True)
            hdr['INTI_Y1'] = y1_img
            hdr['INTI_Y2'] = y2_img
            hdr['INTI_X1'] = x1_img
            hdr['INTI_X2'] = x2_img
        
        else :
            # Sauvegarde en fits de l'image finale
            #cercleC=cercle0
            cercleC=cercle
            frame=np.array(frame, dtype='uint16')
            hdr['NAXIS1']=newiw

            hdr['INTI_XC'] = cercle[0]
            hdr['INTI_YC'] = cercle[1]
            hdr['INTI_R'] = cercle[2]
            # duplication keywords pour BASS2000
            hdr['CENTER_X'] = cercle[0]
            hdr['CENTER_Y'] = cercle[1]
            hdr['SOLAR_R'] = cercle[2]
            
            # recalcul des coordonnées haut et bas du disque
            y1_img,y2_img = detect_bord(frame, axis=1, offset=0,flag_disk=True)
            x1_img,x2_img = detect_bord(frame, axis=0, offset=0,flag_disk=True)
            hdr['INTI_Y1'] = y1_img
            hdr['INTI_Y2'] = y2_img
            hdr['INTI_X1'] = x1_img
            hdr['INTI_X2'] = x2_img

        
        hdr['FILENAME']= basefich+img_suff[k]+'_'+filename_suffixe+".fits"
        
        if auto_crop :
            # modif du 11 mars 2025 car r est trop petit de 5 pixels voir eclipse virtuelle
            """
            if cfg.LG == 1:
                logme('Centre xcc,ycc et rayon : '+str(cercleC[0])+' '+str(cercleC[1])+' '+str(int(r)))
                logme('Coordonnées y1,y2 et x1,x2 disque : '+str(y1_img)+','+str(y2_img)+' '+str(x1_img)+','+str(x2_img))
            else:
                logme('xcc,ycc center and radius : '+str(cercleC[0])+' '+str(cercleC[1])+' '+str(int(r)))
                logme('Coordinates y1,y2 and x1,x2 disk : '+str(y1_img)+','+str(y2_img)+' '+str(x1_img)+','+str(x2_img))
            """
            if cfg.LG == 1:
                logme('Centre xcc,ycc et rayon : '+str(cercleC[0])+' '+str(cercleC[1])+' '+str(cercleC[2]))
                logme('Coordonnees y1,y2 et x1,x2 disque : '+str(y1_img)+','+str(y2_img)+' '+str(x1_img)+','+str(x2_img))
            else:
                logme('xcc,ycc center and radius : '+str(cercleC[0])+' '+str(cercleC[1])+' '+str(cercleC[2]))
                logme('Coordinates y1,y2 and x1,x2 disk : '+str(y1_img)+','+str(y2_img)+' '+str(x1_img)+','+str(x2_img))
        try :
            logme('Solar data : '+str( solar_dict))
            hdr['SEP_LAT']=float(solar_dict['B0'])
            hdr['SEP_LON']=float(solar_dict['L0'])
            hdr['SOLAR_P']=float(ang_P)
            hdr['CAR_ROT']=float(solar_dict['Carr'])
        except:
            pass
        
        t8=time.time()
        dt=t8-t7
        if debug_time : tim.write('autocrop : '+"{:.2f}".format(dt)+'\n')
        
        #-----------------------------------------------------------------
        #---- sauve images _recon, profil, log
        
        # detecte si image unique avec shift pour ajouter _cont au nom
        if len(range_dec)==1 and shift !=0 :
            nomfich=basefich+img_suff[k]+'_cont.fits'
        else:
            
            if flag_dopcont and k == len(range_dec)-1 and not flag_weak:               
                # on est avec image de continuum
                nomfich=basefich+img_suff[k]+'_cont.fits'
                
            else :
                nomfich =basefich+img_suff[k]+'_recon.fits'
        
        DiskHDU=fits.PrimaryHDU(frame,header=hdr)
        DiskHDU.writeto(nomfich, overwrite='True')
        if flag_weak != True and flag_pol != True and data_entete[6] !='Manual':
            DiskHDU.writeto("BASS2000"+os.path.sep+hdr['FILENAME'], overwrite='True')
     
        with  open(basefich+'_log.txt', "w") as logfile:
            logfile.writelines(mylog)
        
        # sauve profile en dat
        if k==0 :
            
            np.savetxt("Complements"+os.path.sep+basefich+'.dat',pro_lamb)
        
        # ajoute l'image a la liste
        frames.append(frame)
        
        t9=time.time()
        dt=t9-t8
        if debug_time : tim.write('sauve : '+"{:.2f}".format(dt)+'\n')
        
        dt=t9-t0
        
        if debug_time : 
            tim.write('total recon : '+"{:.2f}".format(dt)+'\n')
            tim.close()
            
        print('total recon (s): '+"{:.2f}".format(dt))
        
    return frames, hdr, cercleC, range_dec, geom, poly
    