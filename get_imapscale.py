import os
import argparse
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib import pyplot as plt
import numpy as np
import warnings


def Enorm(arr):
    return np.sqrt(np.sum(arr**2))


def get_r(alpha, delta):
    return np.array([np.cos(delta)*np.sin(alpha), 
        np.cos(delta)*np.cos(alpha), 
        np.sin(delta)])


def get_sphere_dist(alpha1, delta1, alpha2, delta2):
    r1 = get_r(alpha1, delta1)
    r2 = get_r(alpha2, delta2)
    return np.arccos(np.dot(r1, r2)/Enorm(r1)/Enorm(r2))


def get_imscale(w, h, W):
    coord1, coord2 = W.wcs_pix2world(np.array([[0, 0],[h, w]]), 1)
    coord1 = np.radians(coord1)
    coord2 = np.radians(coord2)
    angle = np.degrees(get_sphere_dist(coord1[0], coord1[1], 
                                coord2[0], coord2[1]))*3600
    return angle/np.hypot(w, h)


def get_NE(w, h, W):
    coord = W.wcs_pix2world(np.array([[h//2, w//2]]), 1)[0]
    RA  = coord[0]
    DEC = coord[1]
    
    yx1, yx2 = W.wcs_world2pix(np.array([[RA, DEC], [RA, DEC+1]]), 1)
    N = np.array([yx2[0] - yx1[0], yx2[1] - yx1[1]])
    N = N/Enorm(N)
    
    yx1, yx2 = W.wcs_world2pix(np.array([[RA, DEC], [RA+1, DEC]]), 1)
    E = np.array([yx2[0] - yx1[0], yx2[1] - yx1[1]])
    E = E/Enorm(E)
    return N, E


def get_imapscale(fname, out, cbar, cd=0.1, cu=0.1, vmin=None,vmax=None, manual=False):

    try:
        hdu = fits.open(fname)
    except:
        print('Wrong fits')
        return

    data = hdu[0].data
    h, w = data.shape
    header = hdu[0].header

    # plot image
    cm ='Greys'
    fig = plt.figure()
    fig.add_subplot(111)
    mean = np.mean(data[np.where( False == np.isnan(data))])
    std  = np.std(data[np.where( False == np.isnan(data))])
    if manual:
        plt.imshow(data, origin='lower', vmin=vmin, vmax=vmax, cmap=cm)
    else:
        plt.imshow(data, cmap=cm, vmin = mean - std*cd, 
            vmax=mean+std*cu, origin = 'lower')
    plt.xticks([])
    plt.yticks([])

    # plot scale 5 arcmin
    try:
        with warnings.catch_warnings():
            # Ignore a warning on using DATE-OBS in place of MJD-OBS
            warnings.filterwarnings('ignore', 
                    message="'datfix' made the change",
                    category=FITSFixedWarning)
            W = WCS(header)
        imscale = get_imscale(w, h, W) 
        leng = 5*60*imscale
        y0 = 0.1*h
        x0 = w - 0.2*w
        plt.plot([x0, x0+leng], [y0, y0], color='red')
        plt.text(x0+leng//2-2, y0+10, "5'", color='red', fontsize=10)
        # plot equat axis
        y0 = 0.92*h
        x0 = 0.11*w
        N, E = get_NE(w, h, W)
        N = N*leng
        E = E*leng
        plt.arrow(x0, y0, N[0], N[1], width=0.1, head_width=10, 
                color='black')
        plt.text(x0+N[0], y0+N[1]+10, 'N', color='black')
        plt.arrow(x0, y0, E[0], E[1], width=0.1, head_width=10, 
                color='black')
        plt.text(x0+E[0]+10, y0+E[1]-30, 'E', color='black')
    except:
        print('No WCS')


    if cbar:
        plt.colorbar()

    plt.tight_layout()
    plt.savefig(out)
    plt.show()


def run(args):

    cbar   = args.cbar   == 'True'
    manual = args.manual == 'True'

    if args.fits:
        name = args.fits
    else:
        print('No input file!')
        return 

    if args.out:
        out = args.out
    else:
        out = args.fits.split('.')[0]+'.png'

   
    
    if manual:
        if args.vmin:
            vmin = args.vmin
        else:
            print('NO DOWN BORDER')
            return

        if args.vmax:
            vmax = args.vmax
        else:
            print('NO UP BORDER')
            return
        get_imapscale(fname=name, out=out, cbar=cbar, cd=args.cd,
                cu=args.cu, vmin=vmin, vmax=vmax, manual=manual)
        return 0 
    get_imapscale(fname=name, out=out, cbar=cbar, cd=args.cd,
            cu=args.cu, manual=manual)

    return 0

#----------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fits", type=str, help='Input fits file')
    parser.add_argument("--out", type=str, help='Output file')
    parser.add_argument("--cbar", default='False', type=str, help='Plot color bar')
    parser.add_argument("--cd", default=0.1, type=float, help='Coeficient vmin = mean - cd*std')
    parser.add_argument("--cu", default=0.5, type=float, help='Coeficient vmax = mean + cu*std')
    parser.add_argument("--manual", default='False', type=str, help='Handle set vmin/vmax')
    parser.add_argument("--vmin", default=None, type=float, help='manual vmin for plot')
    parser.add_argument("--vmax", default=None, type=float, help='manual vmax for plot')



    args = parser.parse_args()
    run(args)


    
