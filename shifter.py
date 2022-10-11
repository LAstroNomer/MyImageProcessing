#! /usr/bin/env python3

import numpy  as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from reproject import reproject_interp


# Данная программа сдвигает изображение
#  hdu1 к изображению hdu0
#  на основе функции reproject_interp
#  Для обоих кадров необходима WCS
#  После сдвига в итоговом изображении header полностью
#  берется из файла hdu0!
    


if __name__ == '__main__':
 
    home = os.getcwd()
    # Название галактики
    name = 'NGC7625'
    # Изображение к которому выравниваем
    os.chdir('res_v')
    hdu0 = fits.open('v1_sub.fts')
    header = hdu0[0].header
    os.chdir(home)
    
    # V
    os.chdir('res_v')
    data = []
    for filename in os.listdir():
        if filename.__contains__('sub'):
            hdu1 = fits.open(filename)
            array, footprint = reproject_interp(hdu1, hdu0[0].header)
            data.append(array)
    os.chdir(home)
    data = sum(data)/len(data)
    fits.writeto(name+'_v.fits', data, header=header,  overwrite=True)

    # B
    os.chdir('res_b')
    data = []
    for filename in os.listdir():
        if filename.__contains__('sub'):
            hdu1 = fits.open(filename)
            array, footprint = reproject_interp(hdu1, hdu0[0].header)
            data.append(array)
    os.chdir(home)
    data = sum(data)/len(data)
    fits.writeto(name+'_b.fits', data, header=header,  overwrite=True)
    '''
    # R
    os.chdir('res_r')
    data = []
    for filename in os.listdir():
        if filename.__contains__('sub'):
            hdu1 = fits.open(filename)
            array, footprint = reproject_interp(hdu1, hdu0[0].header)
            data.append(array)
    os.chdir(home)
    data = (sum(data))/len(data)
    fits.writeto(name+'_r.fits', data, header=header,  overwrite=True)
    '''
