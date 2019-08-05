import numpy as np
import matplotlib.pyplot as plt
import funcfile_emission as fn  # read function module
import warnings
warnings.filterwarnings("ignore")
#OPR3_10660_blue
fits_file_b = './data/OPR3_60557423_blue.fits'  # the input fits filename # good1
fits_file_r = './data/OPR3_60557423_red.fits'  # the input fits filename
# fits_file_b = './data/OPR3_78138835_blue.fits'  # the input fits filename # bad
# fits_file_r = './data/OPR3_78138835_red.fits'  # the input fits filename
# fits_file_b = './data/OPR3_59990194_blue.fits'  # the input fits filename #good2
# fits_file_r = './data/OPR3_59990194_red.fits'  # the input fits filename

# Good List: 60557423 (0.205076), 61504585 (0.3062), 61524537, 59990194 (0.0178),

# fits_file_b = './data/OPR3_78267792_blue.fits'  # the input fits filename
# fits_file_r = './data/OPR3_78267792_red.fits'  # the input fits filename

z=0.205076 #good1
# z = 1.0164 # bad
# z= 0.0178 #good2

# z = 1.0166
file_b = fn.read_fits(fits_file_b)
file_r = fn.read_fits(fits_file_r)

# print(file_b[1].header)
blue = file_b[1]
blue_std = file_b[2]
red = file_r[1]
red_std = file_r[2]

b_list = blue.data[:]
# b_scaled = blue.data[:]* 10e17
b_list_std = blue_std.data[:]
b_div = b_list / b_list_std

r_list = red.data[:]
# r_scaled = red.data[:]* 10e17
r_list_std = red_std.data[:]
r_div = r_list / r_list_std

wavelength_b = np.arange(9649)*0.25 + 3676
wavelength_r = np.arange(15289)*0.25 + 5772
# wavelength = np.arange(80000)*0.25 + 3676

# fig1 = plt.figure(figsize=(8,6))
# ax1= fn.plot_fig(fig1,211,wavelength_b[:],b_scaled,'Blue Spectra','Spectrum')
# ax2 = fn.plot_fig(fig1,212,wavelength_b[:],b_div[:], '','Spectrum / Error')


# log_w = np.geomspace(5772,9594.25,15289)
# fig2 = plt.figure(figsize=(8,6))
# fn.plot_fig(fig2,211,wavelength_r[:],r_scaled,'Red Spectra','Spectrum',colour1='r')
# fn.plot_fig(fig2,212,wavelength_r[:],r_div[:],'','Spectrum / Error',colour1='r')

xaxis_b = np.arange(0, 9649,1)
xaxis_r = np.arange(0, 15289,1)
xaxis_ori = np.arange(0, 36000, 1)
# print(xaxis[:10])
gal_zip = fn.prepareTemplate(key='SF')

tem_r = fn.makeShiftedTemplate(gal_zip, xaxis_r, z,'r')
tem_b = fn.makeShiftedTemplate(gal_zip, xaxis_b, z, 'b')
tem2 = fn.makeOriginalTemplate(gal_zip, xaxis_ori)


# fig3 = plt.figure(figsize=(8,6))
# plt.plot(wavelength_r,tem_r, alpha = 0.7)
# fn.plot_fig(fig1,211,wavelength_b,tem_b,colour1='r')
# fn.plot_fig(fig2,211,wavelength_r[:],tem_r)
# fn.plot_fig(fig2,212,wavelength_r[:],r_div[:],'','Spectrum / Error',colour1='r')

# fig3 = plt.figure(figsize=(8,6))
# ax3 = fig3.add_subplot(111)
# ax3.xcorr(blue.data[:],blue_std.data[:])
# ax3.xcorr(tem,r_div[:])
# xaxis1 = np.arange(-7644.5, 7644.5, 1)
xaxis1 = np.arange(-15289, 15288, 1)
# xaxis2 = np.geomspace(-1,1)
y_r = np.correlate(tem_r, r_div, 'full') # better CCF
# y = np.correlate(tem, r_div, 'full')
# ax3.plot(xaxis1,y_r)
# ax3.set_xscale('log')

# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)
xaxis3 = np.arange(-9649, 9648, 1)
y_b = np.correlate(tem_b, b_div, 'full') # better CCF



# ax4.plot(xaxis3,y_b)
# a = np.log10(5772)
# b = np.log10(9594.25)
# bins = np.geomspace(a,b,15289)
# print((bins[50]-bins[49]))

# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)
# ax4.plot(bins, r_div)

# wavelength = np.arange(36000)*0.25 + 1000.0
# fig5 = plt.figure(figsize=(8,6))
# ax5 = fig5.add_subplot(111)
# ax5.plot(wavelength, tem2)

'''Logspace and Rebinning'''
e = np.e
ln = np.log
# print(log(e))
log_wvlngth = np.logspace(ln(3676), ln(9594.25), 5745, base = e) # correct
# print(ln(log_wvlngth[100])-ln(log_wvlngth[99]) )
# print(ln(log_wvlngth[5600])-ln(log_wvlngth[5599]) )
print(log_wvlngth[:15])
logunit_wvlngth = ln(log_wvlngth) # works , correct

bins = np.arange(5745)
#  start with blue spectra only. map log scale to log_wvlngth and put value in bin[i]. Refer Bull notes
#  For putting the value in the bin, take weighted mean. At each edge you have to check for fraction of value in that
#   bin and with that fraction*value

calc_log = np.where(np.logical_and((wavelength_b >= log_wvlngth[3]), (wavelength_b < log_wvlngth[4])))
# print(wavelength_b[know[-1:]])
calc_log_index = (np.asarray(calc_log)).flatten()
print(calc_log)
print(wavelength_b[calc_log_index])
frac_r = (log_wvlngth[4] - wavelength_b[calc_log_index[-1]]) / 0.25
print(wavelength_b[calc_log_index[0]])
frac_l = (wavelength_b[calc_log_index[0]] - log_wvlngth[3]) / 0.25
print('fracl',frac_l)
print('fracr', frac_r) # This frac has to be mul by the next lin index val and put in present log bin

num_sum  = 0 ##
den_sum = 0
blue_sum = 0

for i in calc_log_index:
    num_sum = b_list[i]*(1/np.square(b_list_std[i]))
    den_sum = (1/np.square(b_list_std[i]))

blue_sum = num_sum/den_sum
print(blue_sum)
num_sum = num_sum + frac_l * b_list[calc_log_index[0]-1] * (1 / np.square(b_list_std[calc_log_index[0]-1]))
den_sum = den_sum + frac_l * (1 / np.square(b_list_std[calc_log_index[0]-1]))
blue_sum = num_sum/den_sum
print(blue_sum)
num_sum = num_sum + frac_r * b_list[calc_log_index[-1]+1] * (1 / np.square(b_list_std[calc_log_index[-1]+1]))
den_sum = den_sum + frac_r * (1 / np.square(b_list_std[calc_log_index[-1]+1]))
blue_sum = num_sum/den_sum
print(blue_sum)
print(b_list[7:11]) 
# ok 

#put these vals in bins


fig6 = plt.figure(figsize=(14.4,10.8))
fn.plot_fig(fig6,221,wavelength_b[:], b_div, label1='Spectrum')
ax6_1 = fn.plot_fig(fig6,221,wavelength_b, tem_b, colour1='k', title='Blue Spectra', label1='Template', alpha=0.6)
ax6_1.set_xlabel('Wavelength')
fn.plot_fig(fig6,222,wavelength_r[:],r_div, title='Red Spectra',label1='Spectrum',colour1='r')
ax6_2 = fn.plot_fig(fig6,222,wavelength_r[:],tem_r, colour1='k', label1='Template', alpha=0.6, title='Red Spectra')
ax6_2.set_xlabel('Wavelength')
ax6_3 = fn.plot_fig(fig6,223,xaxis3,y_b)
ax6_3.set_ylabel('Blue - CCF')
ax6_3.set_xlabel('Lag (Pixels)')
ax6_4= fn.plot_fig(fig6,224,xaxis1,y_r)
ax6_4.set_ylabel('Red - CCF')
ax6_4.set_xlabel('Lag (Pixels)')
# plt.show()




