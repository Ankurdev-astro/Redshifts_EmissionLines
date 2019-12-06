from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
'''Used to make the statistical plots'''
results_file = fits.open('FITS-results/table_results_ankur_full.fits')
object_count = results_file[1].header[4]
# print(results_file[1].header[8])
# print(results_file[1].header[10])
# print(results_file[1].header[13])
# print(results_file[1].data[0])
# print(results_file[1].data[0][0])  # id
# print(results_file[1].data[0][1])  # z_true
# print(results_file[1].data[0][2])  # z_est
# print(results_file[1].data[0][5:13])  # choices
z_true = []
z_est_ankur = []
z_slinefit = []
id = []
count = 0
best_count = 0
all = []
good = []
good_slinefit = []
good_notslinefit = []
bin_size = 0.01
num_bins = 200
for i in range(object_count):
    id.append(results_file[1].data[i][0])
    z_true_var = results_file[1].data[i][1]

    # print(z_true_var[0][0], 'Hi')
    ankur_est = [val for val in results_file[1].data[i][5:13]]
    # print(type(ankur_est))
    # z_best_ankur = [z_true_var - j for j in results_file[1].data[i][5:13]]
    # [print(abs(z_true_var - j)) for j in results_file[1].data[i][5:13]]
    list = [abs(z_true_var[0][0] - j) for j in results_file[1].data[i][5:13]]
    # list = []
    # for j in results_file[1].data[i][5:13]:
    #     list.append(abs(z_true_var[0][0]))

    # print(list[np.argwhere()])
    # print(z_true_var[0][0], ankur_est[np.argmin(list)],min(list))
    if min(list) < 20:
        count += 1
        if min(list) <= 0.003:
            best_count += 1
        z_true.append(z_true_var[0][0])
        z_slinefit.append(results_file[1].data[i][2][0][0])
        z_est_ankur.append(ankur_est[np.argmin(list)])
        all.append(z_true_var[0][0])
        # exit()
        if min(list) <= 0.003:
            good.append(ankur_est[np.argmin(list)])

        if abs(results_file[1].data[i][2][0][0] - z_true_var[0][0]) <= 0.003:
            good_slinefit.append(results_file[1].data[i][2][0][0])

        if ((abs(results_file[1].data[i][2][0][0] - z_true_var[0][0]) > 0.003) and
                min(list) <= 0.003):
            good_notslinefit.append(ankur_est[np.argmin(list)])

    # print(max(list))
    # print(z_best_ankur)


print(len(all))
print(len(good))
print(len(good_slinefit))
print(len(good_notslinefit))
x = np.arange(0,2,0.001)
# print(x)
plt.figure(figsize=(16,12))
plt.plot(x,x, 'k')
plt.scatter(z_true, z_est_ankur, marker = '+', label = 'Sources')
plt.xlabel(r'True Redshift ($z_{true}$)', fontsize=22)
plt.ylabel(r'Estimated Redshift ($z_{estimated}$)', fontsize = 22)
plt.legend(loc= 'upper left', fontsize = 'xx-large')
# plt.scatter(z_slinefit,z_true, marker = '*', color = 'r')
plt.figure(figsize=(16,12))
plt.hist(all, num_bins, cumulative=True, histtype='step', range = (0,2),
         label = 'All Sources', color = 'k', linewidth=1.2)
plt.hist(good, num_bins, cumulative=True, histtype='step', range = (0,2),
         label=r'Estimated redshifts using CCF ($\Delta z < 0.003$)', color = 'r', linewidth=1.2, linestyle='--')
plt.hist(good_slinefit, num_bins, cumulative=True, histtype='step', range = (0,2),
         label= 'Redshifts estimated by slinefit',linestyle=':', color = 'b', linewidth=1.2)
plt.hist(good_notslinefit, num_bins, cumulative=True, histtype='step', range = (0,2),
         label='Redshifts estimated by CCF and not slinefit',linewidth=2, color = 'r', linestyle='-')
plt.yscale('log')
plt.ylabel('Number of Redshifts N(z)', fontsize=22)
plt.xlabel('Redshift (z)', fontsize = 22)
plt.ylim(1, 10000)
plt.xlim(0, 2)
# plt.figure()
# # plt.hist(res, num_bins)
# plt.bar(all, res.cumcount, width=bin_size)
plt.legend(loc= (0.51,0), fontsize = 'x-large')
plt.show()
