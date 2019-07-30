import numpy as np
fits_file = 'stacked_lockman_6st_csub.fits'  # the input fits filename

def foo():  # test fuction
    print('Anks')

def read_fits(fits_file):
    ''' This function opens the fits file and takes fits filename as parameter'''
    from astropy.io import fits
    HDU_file = fits.open(fits_file)
    return HDU_file

def index_data(lambda1, colour):
    '''This function converts wavelength to
        it's respective index. Colour key: 'b' for blue; 'r' for red '''
    if colour == 'b':
        i = (lambda1 - 3676.0) / 0.25
    elif colour == 'r':
        i = (lambda1 - 5772.0) / 0.25
    return int(i)

def index_data_original(lambda_1):
    '''This function converts wavelength to
        it's respective index.'''
    i1 = (lambda_1 - 1000.0) / 0.25
    return int(i1)

def wavelengthConvert(index,colour):
    '''Convert into Wavelength
        Colour key: 'b' for blue; 'r' for red '''
    if colour == 'b':
        w = (3676 + index * 0.25)
    elif colour == 'r':
        w = (5772 + index * 0.25)
    return w

def plot_fig(fig,subplot,xaxis,yaxis, title = '', label1 = '', colour1 = 'C0', alpha = 1):
    '''This function plots the figure. Specify subplot(rows, columns, target_index)
    and xaxis and yaxis'''
    ax = fig.add_subplot(subplot)
    ax.set_title(title)
    ax.plot(xaxis, yaxis, color = colour1 , label= label1, alpha = alpha)
    if len(label1) != 0:
            ax.legend()
    return ax

def makeGaus(mean, x, ht, std):
    y0 = np.exp(-0.5 * (((x - mean) / std) ** 2))
    y = y0 * ht
    return y

def prepareTemplate(key='SF'):
    if key == 'SF':
        galaxyLines = np.genfromtxt('LineRatio_nodust.txt', skip_header=1, usecols=(1))
        galaxyStr = np.genfromtxt('LineRatio_nodust.txt', skip_header=1, usecols=(2))
        # galaxyStr_Norm = galaxyStr / np.max(galaxyStr)
        galaxyStr_Norm = galaxyStr * 5
        gal_zipper = zip(galaxyLines, galaxyStr_Norm)
        gal_ziptmp = list(gal_zipper)
        gal_zip = sorted(gal_ziptmp,key=lambda x: x[1], reverse= True)
        # gal_zip.reverse()
    elif key == 'AGN':
        AGNLines = [1023., 1035, 1216, 1240, 1336, 1402, 1549, 1640, 1663, 1909, 2326, 2424, 2426, 2470, 2798, 3113,
                    3203, 3346, 3426, 2727, 2869, 3889, 3970, 3967, 4072, 4102, 4340, 4363, 4686, 4861, 4959, 5007,
                    6563]
        AGNstr = [85., 192, 3100, 154, 37, 163, 364, 318, 72, 177, 92, 41, 49, 41, 78, 11, 14, 20,
                  69, 364, 82, 19, 21, 26, 16, 22, 24, 8, 20, 100, 307, 866, 310]
        galaxyStr_Norm = AGNstr / np.max(AGNstr)
        gal_zipper = zip(AGNLines, galaxyStr_Norm)
        gal_ziptmp = list(gal_zipper)
        gal_zip = sorted(gal_ziptmp, key=lambda x: x[1], reverse=True)
        #gal_zip.reverse()
    return gal_zip

def makeShiftedTemplate(all_Lines_Template, xaxis, z, colour = ''):
    ySum = 0
    if colour == 'r':
        for everyLine in all_Lines_Template:
            if 5772 < everyLine[0] * (z + 1) < 9594.25:
                shifted = index_data(everyLine[0] * (z + 1), 'r')
                template = makeGaus(shifted, xaxis, everyLine[1], std=10)
                ySum += template
                # print(everyLine[0] * (z + 1), '\t',  everyLine[1])
        return ySum
    if colour == 'b':
        for everyLine in all_Lines_Template:
            if 3676.0 < everyLine[0] * (z + 1) < 6088.25:
                shifted = index_data(everyLine[0] * (z + 1), 'b')
                template = makeGaus(shifted, xaxis, everyLine[1], std=10)
                ySum += template
                # print(everyLine[0] * (z + 1), '\t',  everyLine[1])
        return ySum

def makeOriginalTemplate(all_Lines_Template, xaxis):  # function not in use
    ySum = 0
    for everyLine in all_Lines_Template:
        template = makeGaus(index_data_original(everyLine[0]), xaxis, everyLine[1],std=10)
        ySum += template
    return ySum