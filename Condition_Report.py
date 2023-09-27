import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time as aTime


def Condition_Report(Object, Time, fil, Path2Save):

    plt.switch_backend('pdf')  # для фикса какой-то тупой ошибки в нарнии pyplot
    fig, axs = plt.subplots(7, 1, sharex=True, figsize=(6, 7), dpi=125)
    Title = Object
    Title += '\n'
    Title += 'Robophot, filter=' + fil

    fig.suptitle(Title, fontsize=8)
    pos = [0.125, 0.79, 0.8, 0.1]
    # axs[0].plot(Time['JD'], Time['AIRMASS'], 'b.')
    # axs[0].set_ylabel('airmass', fontsize=6)

    axs[1].plot(Time['JD'], Time['Sky'], 'r.')
    axs[1].set_ylabel('sky (ADU)', fontsize=6)

    axs[2].plot(Time['JD'], Time['EXTINCT'], 'b.')
    axs[2].set_ylabel('extinction (mag)', fontsize=6)

    # axs[3].plot(Time['JD'], Time['SEXFWHM']*1.85, 'r.')
    # axs[3].set_ylabel('FWHM (arcsec)', fontsize=6)

    # axs[4].plot(Time['JD'], Time['SEXELL'], 'b.')
    # axs[4].set_ylabel('ellipticity', fontsize=6)

    axs[5].plot(Time['JD'], Time['Max'], 'r.')
    axs[5].set_ylabel('Max count (ADU)', fontsize=6)

    axs[6].plot(Time['JD'], Time['X']-np.mean(Time['X']), 'r.', label='X Shift')
    axs[6].plot(Time['JD'], Time['Y']-np.mean(Time['Y']), 'b.', label='Y Shift')
    axs[6].legend(loc=0, fontsize=6)
    axs[6].set_ylabel('shift (pix)', fontsize=6)

    locs, labels = plt.xticks()
    t = aTime(locs, format='jd')
    x_ticks_labels = []
    for x in t:
        x_ticks_labels.append(str(x.iso).split(' ')[1].split('.')[0])
    axs[6].set_xticklabels(x_ticks_labels, rotation='vertical', fontsize=5)
    axs[6].set_xlabel('Date-Time (UTC), ' + Time['DATE-OBS'][0].split('T')[0], fontsize=6)

    for ii in range(0, len(axs)):
        axs[ii].set_position(pos)
        pos[1] = pos[1]-0.115
        # axs[ii].axvspan(T0.jd, T1.jd, facecolor='k', alpha=0.2)
        axs[ii].tick_params(axis='both', labelsize=6, direction='in')
        axs[ii].grid()
        
    plt.savefig(Path2Save+r'\condition.pdf')
    # plt.show()
