from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._0errors import *
from ._0imports import *


class Functions:

    def __init__(self):
        pass

    @staticmethod
    def open_dict(file_path):
        return pickle.load(open(file_path, 'rb'))

    @staticmethod
    def save_dict(dict_data, file_path):
        pickle.dump(dict_data, open(file_path, 'wb'), protocol=2)

    @staticmethod
    def sci(hdulist):
        return np.where(np.array([ff.name for ff in hdulist]) == 'SCI')[0]

    @staticmethod
    def err(hdulist):
        return np.where(np.array([ff.name for ff in hdulist]) == 'ERR')[0]

    def sci_err(self, hdulist):
        return np.swapaxes([self.sci(hdulist), self.err(hdulist)], 0, 1)

    @staticmethod
    def fits_like(hdulist):
        copy_fits = [ff.copy() for ff in hdulist]
        return pf.HDUList(copy_fits)

    @staticmethod
    def save_figure(directory, name, figure=None, transparent=True):
        if not os.path.isdir(directory):
            os.makedirs(directory)
        if not name:
            raise IraclisFileError('No name for figure given')
        else:
            if not figure:
                plt.savefig(directory + '/{0}.eps'.format(name), bbox_inches='tight', transparent=transparent)
                plt.close()
            else:
                figure.savefig(directory + '/{0}.eps'.format(name), bbox_inches='tight', transparent=transparent)

    @staticmethod
    def adjust_ticks():

        xlim1, xlim2 = plt.xlim()
        xticks = plt.xticks()[0]
        digit = max([len(ff) for ff in str(np.abs(xticks)-np.int_(np.abs(xticks)))[1:-1].split()]) - 2
        xticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in xticks]
        dxticks = xticks[1] - xticks[0]

        if xticks[1] - dxticks / 2 < xlim1:
            new_xlim1 = xticks[1] - dxticks / 2
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]
        else:
            new_xlim1 = xticks[0] - dxticks / 2

        if xticks[-2] + dxticks / 2 > xlim2:
            new_xlim2 = xticks[-2] + dxticks / 2
            xticks = xticks[:-1]
            xticklabels = xticklabels[:-1]
        else:
            new_xlim2 = xticks[-1] + dxticks / 2

        plt.xticks(xticks, xticklabels, fontsize=15)
        plt.xlim(new_xlim1, new_xlim2)

        ylim1, ylim2 = plt.ylim()
        yticks = plt.yticks()[0]
        digit = max([len(ff) for ff in str(np.abs(yticks)-np.int_(np.abs(yticks)))[1:-1].split()]) - 2
        yticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in yticks]
        dyticks = yticks[1] - yticks[0]

        if yticks[1] - dyticks / 2 < ylim1:
            new_ylim1 = yticks[1] - dyticks / 2
            yticks = yticks[1:]
            yticklabels = yticklabels[1:]
        else:
            new_ylim1 = yticks[0] - dyticks / 2

        if yticks[-2] + dyticks / 2 > ylim2:
            new_ylim2 = yticks[-2] + dyticks / 2
            yticks = yticks[:-1]
            yticklabels = yticklabels[:-1]
        else:
            new_ylim2 = yticks[-1] + dyticks / 2

        plt.yticks(yticks, yticklabels, fontsize=15)
        plt.ylim(new_ylim1, new_ylim2)

    @staticmethod
    def adjust_ticks_ax(ax):

        xlim1, xlim2 = ax.get_xlim()
        xticks = ax.get_xticks()
        digit = max([len(ff) for ff in str(np.abs(xticks)-np.int_(np.abs(xticks)))[1:-1].split()]) - 2
        xticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in xticks]
        dxticks = xticks[1] - xticks[0]

        if xticks[1] - dxticks / 2 < xlim1:
            new_xlim1 = xticks[1] - dxticks / 2
            xticks = xticks[1:]
            xticklabels = xticklabels[1:]
        else:
            new_xlim1 = xticks[0] - dxticks / 2

        if xticks[-2] + dxticks / 2 > xlim2:
            new_xlim2 = xticks[-2] + dxticks / 2
            xticks = xticks[:-1]
            xticklabels = xticklabels[:-1]
        else:
            new_xlim2 = xticks[-1] + dxticks / 2

        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=15)
        ax.set_xlim(new_xlim1, new_xlim2)

        ylim1, ylim2 = ax.get_ylim()
        yticks = ax.get_yticks()
        digit = max([len(ff) for ff in str(np.abs(yticks)-np.int_(np.abs(yticks)))[1:-1].split()]) - 2
        yticklabels = ['{1:.{0}f}'.format(digit, ii) for ii in yticks]
        dyticks = yticks[1] - yticks[0]

        if yticks[1] - dyticks / 2 < ylim1:
            new_ylim1 = yticks[1] - dyticks / 2
            yticks = yticks[1:]
            yticklabels = yticklabels[1:]
        else:
            new_ylim1 = yticks[0] - dyticks / 2

        if yticks[-2] + dyticks / 2 > ylim2:
            new_ylim2 = yticks[-2] + dyticks / 2
            yticks = yticks[:-1]
            yticklabels = yticklabels[:-1]
        else:
            new_ylim2 = yticks[-1] + dyticks / 2

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=15)
        ax.set_ylim(new_ylim1, new_ylim2)


functions = Functions()
