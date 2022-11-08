#! /usr/bin/env python3

# Данная программа предназачена для проведения первичной обработки изображений
# Примечания по работе кода:
#   1. n_bias :: число кадров bias для суммирования. Я предполагаю, 
#     что все исходные bias файлы лежат в той же дирректории, 
#     где находится скрип и называются "bias_%s.fits", 
#     где %s номер фитса начиная с 1. Результирующий файл "bias.fits   
#     HEADER не сохраняется!
#
#   2. filtrs :: Пути к папкам с фильтрами.Для каждого свой.
#     Из неё берутся все файлы с расширением fits, 
#     если все кадры лежат в одной папке, 
#     то этот кусок придется переписать.
#
#   3. flats :: Файлы плоского поля с путями к ним. 
#      У меня лежат в той же папке, где скрипт
#
#   4. Фон неба я учитывал по аналогии с тем как это делает SExtractor
#      Мне не нравится то, что получается, но как довести до идеала я не знаю.
#      Возможно будет лучше закоментить этот этап и сделать оценку фона более 
#      привычным методом.
#
#   5. Результаты обработки каждого фильтра сохраняются в свои папки res_%s,
#      где %s название фильтра. Эти папки надо создать руками. 
#      Можно подправить код, чтобы создавались автоматически.
#
#

from astropy.io import fits
import os
import numpy  as np
import matplotlib.pyplot as plt
from scipy import odr
from photutils.segmentation import detect_threshold, detect_sources
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from photutils.utils import circular_footprint
import subprocess
from scipy.signal import medfilt2d


def cre_bias(num, fname='bias'):
    # Функция усредняет все bias, что есть в папке  \ 1. \
    bias0 = 0
    for i in range(1, num+1):
        bias0 += fits.getdata(fname+'_%s.fts' %i)
    bias = bias0/num
    fits.writeto('bias.fts', bias, overwrite=True)
    return bias

def sub_sky(data):
    # Функция строит карту всех источников и находит 2D изображение  фона неба
    from astropy.stats import SigmaClip
    from photutils.background import Background2D, MedianBackground
   
    # Сырой поиск источников
    threshold = detect_threshold(data, nsigma=2.0)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    convolved_data = convolve(data, kernel, normalize_kernel=True)

    # Создание маски для повышения качества оценки фона
    segm = detect_sources(convolved_data, threshold, npixels=5)
    footprint = circular_footprint(radius=35)
    mask = segm.make_source_mask(footprint=footprint)
    print('before sky fix', sigma_clipped_stats(data, sigma=3.0, mask=mask))
    
    #Оценка фона
    sigma_clip = SigmaClip(sigma=3.)
    Bkg = Background2D(data, (64,64), mask=mask,  filter_size=(5,5), sigma_clip=sigma_clip, fill_value=0.0) 
    print('after sky fix', sigma_clipped_stats(data-Bkg.background, sigma=3.0, mask=mask))
    return Bkg, segm

def find_galaxy(segm):
    # Функция находит на маске галактику

    num_of_obj = np.max(segm.data)
    print('Objs = %s' %num_of_obj)
    max_area = 0
    gal_ind = 0
    for i in range(1, num_of_obj+1):
        area = segm.get_area(i)
        if area > max_area:
            gal_ind = i
            max_area = area
    print('Gal ind = %s' %gal_ind)
    segm.remove_label(gal_ind)
    return segm



def main(nbias, filtrs, flats):

    # \ 1. \
    bias = cre_bias(nbias)

    for i in range(len(filtrs)):
        filt = filtrs[i]

        # Поиск плохих пикселов. Медианная маскировка. Усреднение 
        flat = fits.getdata(flats[i])
        bady, badx = np.where((flat == 0))
        
        for y, x in zip(bady, badx):
            flat[y, x] = np.median(flat[y-1:y+1, x-1:x+1])
        
        flat = flat/np.mean(flat)
        

        # Загрузка галактики. Вычитание биаса. Деление на флат
        # Медианная фильтрация 3х3 пиксела
        # Сохранение файла результата
        # Анализ фона неба и сохранение результата в фитсе
        # Сохранение изображения без фона
        for filename in os.listdir(filt): 
            if filename.__contains__('.fits'):
                print(filename, flats[i])
                hdulist = fits.open('%s/%s' %(filt,filename))
                header  = hdulist[0].header
                data    = hdulist[0].data
                #t       = hdulist[0].header['EXPTIME']

                data_bias = data - bias

                data_res  = data_bias/flat#/t

                dara_res = medfilt2d(data_res)           
                
                name = filename.split('.fits')[0]
                fits.writeto('res_%s/%s_res.fts' %(filt, name), data_res, header=header, overwrite=True)
                       
                bkg, segm = sub_sky(data_res)
                fits.writeto('res_%s/%s_bkg.fts' %(filt, name), bkg.background, header=header, overwrite=True)
            
                data_sub = data_res - bkg.background
                fits.writeto('res_%s/%s_sub.fts' %(filt, name), data_sub, header=header, overwrite=True)

            
if __name__ == '__main__':

    n_bias = 10
    filtrs = [ 'b', 'v']
    flats = [ 'flat_%s.fts' %i for i in filtrs]

    main(n_bias, filtrs, flats)
