import numpy as np
import matplotlib.pyplot as plt
import funcfile_emission_x1 as fn  # read function module
import os
import warnings
warnings.filterwarnings("ignore")
'''File made for accessing x1 dataset and uses all sub-files'''
files = os.listdir('./data/')
file_number_id = 73079830
z = 0.842465
print('Redshift = ', z)
list_filenames = []
targetfiles = []
b_list = 0
b_list_std = 0
r_list = 0
r_list_std = 0
for filename in files:
    # break
    if filename.endswith('.txt') and filename.startswith('OPR3'):
        fileid = (filename.split('.txt')[0])
        file_number = fileid.split('OPR3_')[1]
        list_filenames.append((int(file_number), filename))

for file_id in list_filenames:
    if file_id[0] == file_number_id:
        text_file_id = file_id[1]

with open(('./data/' + str(text_file_id)), "r") as thisfile:
    for row in thisfile.readlines():
        targetfiles.append(row.split('\n')[0])
numberof_obsv = int(len(targetfiles)/2)
print(numberof_obsv, ': #Observation Counts')

for loop in range(numberof_obsv):
    fits_file_b = './data/' + str(targetfiles[0 + loop])
    fits_file_r = './data/' + str(targetfiles[numberof_obsv + loop])
    file_b = fn.read_fits(fits_file_b)
    file_r = fn.read_fits(fits_file_r)
    blue = file_b[1]
    blue_std = file_b[2]
    red = file_r[1]
    red_std = file_r[2]
    b_list += blue.data[:]
    b_list_std += blue_std.data[:]
    r_list += red.data[:]
    r_list_std += red_std.data[:]

b_list[fn.index_data(6000,'b'):] = 0
b_list_std[fn.index_data(6000,'b'):] = (1.0e+09)*numberof_obsv


print(len(b_list), len(r_list))
wavelength_b = np.arange(9649)*0.25 + 3676
wavelength_r = np.arange(15289)*0.25 + 5772


xaxis_b = np.arange(0, 9649,1)
xaxis_r = np.arange(0, 15289,1)

gal_zip = fn.prepareTemplate(key='SF')

tem_r = fn.makeShiftedTemplate(gal_zip, xaxis_r, z,'r')
tem_b = fn.makeShiftedTemplate(gal_zip, xaxis_b, z, 'b')
# tem2 = fn.makeOriginalTemplate(gal_zip, xaxis_ori)
'''Simulating and adding noise to template'''
rand_b = np.random.normal(size=9649)
rand_r = np.random.normal(size=15289)
empty_b = np.where(b_list_std> 1e8)
empty_r = np.where(r_list_std> 1e8)
tem_rand_b = tem_b + b_list_std*rand_b
tem_rand_r = tem_r + r_list_std*rand_r
tem_rand_b[empty_b] = 0
tem_rand_r[empty_r] = 0
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

# xaxis2 = np.geomspace(-1,1)


# y = np.correlate(tem, r_div, 'full')
# ax3.plot(xaxis1,y_r)
# ax3.set_xscale('log')

# fig4 = plt.figure(figsize=(8,6))
# ax4 = fig4.add_subplot(111)


# ax4.plot(xaxis3,y_b)
# a = np.log10(5772)
# b = np.log10(9594.25)
# bins = np.geomspace(a,b,15289)
# print((bins[50]-bins[49]))

# fig = plt.figure(figsize=(8,6))
# # ax4 = fig4.add_subplot(111)
# # ax4.plot(bins, r_div)
#
# # wavelength = np.arange(36000)*0.25 + 1000.0
# # fig5 = plt.figure(figsize=(8,6))
# # ax5 = fig5.add_subplot(111)
# # ax5.plot(wavelength, tem2)4
log_wvlngth_b, rebin_val_b, rebin_ivar_b = fn.rebin(3676, 6088.25, 1515, wavelength_b, b_list, b_list_std) # 1515

log_wvlngth_r, rebin_val_r, rebin_ivar_r = fn.rebin(5772, 9594.25, 1526, wavelength_r, r_list, r_list_std) #1526
# print(rebin_val_r[:50])
'''template log bin'''
log_wvlngth_temb, rebin_val_temb, rebin_ivar_temb = fn.rebin(3676, 6088.25, 1515, wavelength_b, tem_rand_b, b_list_std)
log_wvlngth_temr, rebin_val_temr, rebin_ivar_temr = fn.rebin(5772, 9594.25, 1526, wavelength_r, tem_rand_r, r_list_std)

pix_b = (np.log(6088.25) - np.log(3676))/1515
pix_r = (np.log(9594.25) - np.log(5772))/1526  # good above 0.9
# print(pix_r)
# log_wvlngth_temb, rebin_val_temb, rebin_ivar_temb = fn.rebin(3676, 6088.25, 1515, wavelength_b, tem_b, )
'''correlate'''
xaxis_corr_b = np.arange(-1514, 1515, 1) #(-1514, 1515, 1)
ccf_b = np.correlate(rebin_val_temb, rebin_val_b, 'full') # better CCF
xaxis_corr_r = np.arange(-1524, 1527, 1) #(-1524, 1527, 1)
ccf_r = np.correlate(rebin_val_temr, rebin_val_r, 'full') # better CCF
print(np.argmax(ccf_r))
print(((len(ccf_r)-1)/2 )-1)
pix_set1 = 2.144698713e-4  # good around 0.5
pix_set2 = 2.648662384e-4  # good at 0.7
lag_r = (np.argmax(ccf_r) - (((len(ccf_r)-1)/2 )-1))
print(lag_r)

z_guess_r = (10)**(pix_r * lag_r)-1
#z_guess_r = (10)**(pix_r * lag_r)*(1 + (30 / 3e5))-1
lag_list_b = []
lag_list_r = []
'''for plotting'''
for i in xaxis_corr_r:
    lag = (10) ** (pix_r * i) - 1
    lag_list_r.append(lag)
for i in xaxis_corr_b:
    lag = (10) ** (pix_r * i) - 1
    lag_list_b.append(lag)
print(round(z_guess_r, 5))
print(((np.ndarray.max(ccf_b)*np.ndarray.mean(ccf_b))/np.std(ccf_b))*1e32,'-- Blue')
print(((np.ndarray.max(ccf_r)*np.ndarray.mean(ccf_r))/np.std(ccf_r))*1e32, '-- Red')
# r = 1.5 is the min val we can trust

# fig9 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig9,111,log_wvlngth_b, rebin_val_b, label1='Log')
#
# fig10 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig10,111,wavelength_b, b_list, label1='Linear')
#
# fig91 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig91,111,log_wvlngth_r, rebin_val_r, colour1='r', label1='Log')
#
# fig101 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig101,111,wavelength_r, r_list,colour1='r', label1='Linear')
# plt.show()

fig6 = plt.figure(figsize=(14.4,10.8))
# ax6_1 = fn.plot_fig(fig6,111,wavelength_b[:], b_list, label1='Spectrum')
ax6_1 = fn.plot_fig(fig6,111,log_wvlngth_b, rebin_val_b, label1='Spectrum')
if isinstance(tem_b,int):
    tem_b = np.zeros(9649)
fn.plot_fig(fig6,111,log_wvlngth_temb, rebin_val_temb, colour1='k', title='Blue Spectra', label1='Template', alpha=0.6)
ax6_1.set_xlabel('Wavelength')
ax6_1.set_ylim(top=0.85e-16, bottom=-0.5e-16)
# fn.plot_fig(fig6,222,wavelength_r[:],r_div, title='Red Spectra',label1='Spectrum',colour1='r')
fig7 = plt.figure(figsize=(16,12))
# fn.plot_fig(fig7,111,wavelength_r[:], r_list, label1='Spectrum', colour1='r')
fn.plot_fig(fig7,111,log_wvlngth_r, rebin_val_r, title='Red Arm Spectra',label1='Spectrum',colour1='r')
if isinstance(tem_r,int):
    tem_r = np.zeros(15289)
ax6_2 = fn.plot_fig(fig7,111,log_wvlngth_temr,rebin_val_temr,
                    colour1='k', label1='Computer Generated Template', alpha=0.6)

ax6_2.annotate('Redshift (z) = 0.842', xy=(0.02, 0.75), xycoords='axes fraction', fontsize = 15)
ax6_2.annotate(r'$[O III]$', xy=(0.85, 0.875), xycoords='axes fraction', fontsize = 15)
ax6_2.annotate(r'$[O II]$', xy=(0.29, 0.41), xycoords='axes fraction', fontsize = 15)
ax6_2.annotate(r'$[H_\gamma]$', xy=(0.56, 0.46), xycoords='axes fraction', fontsize = 15)
ax6_2.annotate(r'$[H_\beta]$', xy=(0.785, 0.51), xycoords='axes fraction', fontsize = 15)
ax6_2.set_title('Red Arm Spectra', fontsize = 22)
ax6_2.set_xlabel('Wavelength (in $\AA$)',  fontsize = 22)
ax6_2.set_ylabel('Flux (in $erg/cm^2/s$)',  fontsize = 22)
ax6_2.set_ylim(top=1.2e-16, bottom=-0.1e-16)
ax6_2.legend(loc= 'upper left', fontsize = 'xx-large')
# ax6_3 = fn.plot_fig(fig6,223,xaxis3,y_b)
# ax6_3.set_ylabel('Blue - CCF')
# ax6_3.set_xlabel('Lag (Pixels)')
# ax6_4= fn.plot_fig(fig6,224,xaxis1,y_r)
# ax6_4.set_ylabel('Red - CCF')
# ax6_4.set_xlabel('Lag (Pixels)')

'''CCF Plots'''
# fig8 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig8,111,lag_list_b,ccf_b)
# fig9 = plt.figure(figsize=(14.4,10.8))
# fn.plot_fig(fig9,111,lag_list_r,ccf_r)
#

plt.show()




