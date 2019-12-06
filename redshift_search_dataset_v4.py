import numpy as np
import os
import funcfile_emission_x1 as fn  # read function module
import ccf_algo_v3 as ccf  # ccf algorithm file
from astropy.table import Table, Column
trunc = lambda f: f - f % 0.01
from astropy.io import fits
import warnings
# warnings.filterwarnings("ignore")
'''This file is made for accessing x1_dataset. Designed to run full dataset on the Linux Herts Machine'''
'''Allow ~more than a day to run. Made on Sept 5th, 2019'''

'''To run only for SF galaxies'''
fits_results= 'results_default.fits'
list_file = fits.open(fits_results)
target_id = list_file[1].data[:1][0][0]
agnType = list_file[1].data[:1][0][9]
z_True = list_file[1].data[:1][0][1]
z_SciEst = list_file[1].data[:1][0][2] #Estimated redshift from results file
pointing = list_file[1].data[:1][0][6]
SFtype = list_file[1].data[:1][0][8]
# print(target_id[:])
'''...................................................'''
list_filenames = []
z_True_list = []
z_SciEst_list = []
pointing_list = []
SFtype_list = []
# thefile = open('redshift_output.txt', 'w+')
file_count = 0
files = os.listdir('./data/')
search_fileList = []
for filename in files:
    # break
    if filename.endswith('.txt') and filename.startswith('OPR3'):
        fileid = (filename.split('.txt')[0])
        file_number = fileid.split('OPR3_')[1]
        list_filenames.append((int(file_number), filename))
for randf in list_filenames:
    # break
    search_fileid = randf[1].split('.txt')[0]
    search_fileList.append(search_fileid)
print('Start')
listof_fileid = []
listof_redshifts = []
count = 0
for fileid in search_fileList:
    # break
    located = np.argwhere(fileid==target_id)
    if agnType[located] == 0 and z_True[located] <= 2:  # agn type
        # print(fileid.split('OPR3_')[1])
        # print(slicepart[located])
        count += 1

        '''Big '2-tab' Loop starts here'''
# file_number = 83769217
        file_number = fileid.split('OPR3_')[1] #for dataset run uncomment this
        # print(file_number)
        # print(count)

        fits_file_b = './data/OPR3_' + str(file_number) + '_blue.fits'  # the input fits filename
        fits_file_r = './data/OPR3_' + str(file_number) + '_red.fits'  # the input fits filename


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
        # b_div = b_list / (b_list_std)

        r_list = red.data[:]
        # r_scaled = red.data[:]* 10e17
        r_list_std = red_std.data[:]
        # r_div = r_list / r_list_std
        # print(len(b_list), len(r_list))
        wavelength_b = np.arange(9649)*0.25 + 3676
        wavelength_r = np.arange(15289)*0.25 + 5772
        # wavelength = np.arange(80000)*0.25 + 3676


        xaxis_b = np.arange(0, 9649,1)
        xaxis_r = np.arange(0, 15289,1)
        # print(xaxis[:10])
        gal_zip = fn.prepareTemplate(key='SF')
        pix_r = (np.log(9594.25) - np.log(5772))/1526
        pix_b = (np.log(6088.25) - np.log(3676))/1515
        z_guess_flag = 0
        z_guess_list1 = []
        z_guess_list2 = []
        z_guess_list3 = []
        final_redshifts = []
        for i in range(21):  # checking till z = 2
            z = round(i * 1e-1,2)
            # print(i,z)
            z_lag_r, r_blue, r_red = ccf.algorithm(z,gal_zip,xaxis_b,xaxis_r,b_list,b_list_std,
                                                   r_list, r_list_std, wavelength_b, wavelength_r)
            if(z_lag_r>= 0.0 and z_lag_r<=0.2 and r_red>0.2):
                # print(z, 'Redshift ', r_blue, 'Blue ', r_red, 'Red', z_lag_r, 'z_lag')
                z_guess_list1.append((z, z_lag_r))
            z_guess_flag = 1


        if(z_guess_flag == 1):
            # print('Round Two')
            for trial_z, trialz_lag_r in z_guess_list1:
                z = trial_z - trialz_lag_r + 0.03
                # if trial_z < 0.5:
                #     z = trial_z - trialz_lag_r*0.6
                if z <= 0:
                    continue
                # print(z, 'trial z')

                z_lag_r, r_blue, r_red = ccf.algorithm(z,gal_zip,xaxis_b,xaxis_r,b_list,b_list_std,
                                                   r_list, r_list_std, wavelength_b,wavelength_r)
                # if (z_lag_r >= -0.06 and z_lag_r <= 0.09 and r_red>0.2):
                if (z_lag_r > 0 and z_lag_r <= 0.09 and r_red > 0.2):
                    # print(z, 'Redshift ', r_blue, 'Blue ', r_red, 'Red', z_lag_r, 'z_lag', (z - z_lag_r), 'final z R2')
                    z_guess_list2.append((z, z_lag_r))
            z_guess_flag = 2

        if(z_guess_flag == 2):
            # print('Round Three')
            for trial_z, trialz_lag_r in z_guess_list2:
                z = trial_z - trialz_lag_r
                # if trial_z < 0.5:
                #     z = trial_z - trialz_lag_r*0.6
                if z <= 0:
                    continue
                # print(z, 'trial z')

                z_lag_r, r_blue, r_red = ccf.algorithm(z,gal_zip,xaxis_b,xaxis_r,b_list,b_list_std,
                                                   r_list, r_list_std, wavelength_b,wavelength_r)
                # if (z_lag_r >= -0.06 and z_lag_r <= 0.09 and r_red>0.2):
                # print(z, 'Redshift ', r_blue, 'Blue ', r_red, 'Red', z_lag_r, 'z_lag', (z - z_lag_r), 'final z R2')
                if (z_lag_r > -0.01 and z_lag_r <= 0.009 and r_red > 0.2):
                    # print(z, 'Redshift ', r_blue, 'Blue ', r_red, 'Red', z_lag_r, 'z_lag', (z - z_lag_r), 'final z R2')
                    z_guess_list3.append((z, z_lag_r))
            z_guess_flag = 3
        # print(z_guess_list3)
        test_append = []
        if(z_guess_flag == 3):
            # print('Round Four')
            for trial_z, trialz_lag_r in z_guess_list3:
                z = trial_z - trialz_lag_r
                # print(z)
                if z <= 0:
                    continue
                z_lag_r = 999  # arbitrary number
                count = 0
                while (z_lag_r > 0.001 or z_lag_r < -0.001 and count <6):
                    count += 1
                    z_lag_r, r_blue, r_red = ccf.algorithm(z, gal_zip, xaxis_b, xaxis_r, b_list, b_list_std,
                                                           r_list, r_list_std, wavelength_b, wavelength_r)
                    z = z - z_lag_r
                    # print(count)
                    # print(z, 'Redshift ', z_lag_r, 'z_lag', (z - z_lag_r), 'final z R2')
                if count < 6:
                    # print(z_lag_r, 'z_lag', (z - z_lag_r), 'Final Estimated Redshift')
                    est_redshift = z - z_lag_r
                    if round(trunc(est_redshift), 2) not in test_append:
                        test_append.append(round(trunc(est_redshift),2))
                        final_redshifts.append((round(est_redshift, 5)))
        file_count += 1
        while len(final_redshifts) < 8:
            final_redshifts.append(-99)
        listof_redshifts.append(final_redshifts)
        listof_fileid.append(fileid)
        z_True_list.append(z_True[located])
        z_SciEst_list.append(z_SciEst[located])
        pointing_list.append(pointing[located])
        SFtype_list.append(SFtype[located])
        print(file_number, final_redshifts)
        print('File Count :', file_count)
        # thefile.write('{0} {1}\n'.format(fileid, final_redshifts))
        # print(listof_fileid)
        # print(listof_redshifts)

'''Big '2-tab' Loop ends here'''
# thefile.close()
table = Table(rows=listof_redshifts, names=['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8'])
new_col = Column(SFtype_list, name='SFTYPE')
table.add_column(new_col, index=0)
new_col = Column(pointing_list, name='POINTING')
table.add_column(new_col, index=0)
new_col = Column(z_SciEst_list, name='Z_SCI_EST')
table.add_column(new_col, index=0)
new_col = Column(z_True_list, name='ZTRUE')
table.add_column(new_col, index=0)
new_col = Column(listof_fileid, name='ID')
table.add_column(new_col, index=0)
# print(table)
table.write('table_results_ankur.fits', format='fits')
exit()