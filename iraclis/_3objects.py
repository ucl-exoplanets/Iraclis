from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._2variables import *


class DataSet:

    def __init__(self, input_data, direct_image=None):

        if direct_image is not None:
            if isinstance(direct_image, str):
                if os.path.isfile(direct_image):
                    direct_image = pf.open(direct_image)
                else:
                    raise IraclisFileError('No such file {0}'.format(input_data))
            else:
                raise IraclisFileError('Please give a file name for the direct image or leave it empty.')

        if not isinstance(input_data, str):
            raise IraclisFileError('Please give a file or directory name for the input data.')

        elif os.path.isfile(input_data):

            self.file_names = []
            self.spectroscopic_images = [pf.open(input_data)]
            if not direct_image:
                self.direct_image = None
            else:
                self.direct_image = direct_image
            self.splitted = False
            self._data_set_directory_path = None

            if self.spectroscopic_images[0][0].header[variables.observation_type.keyword] != 'SPECTROSCOPIC':
                raise IraclisFileError('A single direct image is not s valid input dataset.')

        elif os.path.isdir(input_data) and len(glob.glob(os.path.join(input_data, '*', ''))) == 0:

            nsamp = []
            final_list = []
            direct_image = None

            files = sorted(glob.glob(os.path.join(input_data, '*.fits')))
            for i in files:
                with pf.open(i) as j:
                    if j[0].header[variables.observation_type.keyword] == 'SPECTROSCOPIC':
                        final_list.append([j[0].header[variables.exposure_start.keyword],
                                           os.path.split(i)[1], plc.copy_fits(j)])
                        nsamp.append(j[0].header[variables.total_samples.keyword])
                    elif not direct_image:
                        direct_image = pf.open(i, mode='update')

            nsamps = [int(np.median(np.array(nsamp)))]

            final_list.sort()
            list_of_times, list_of_files, list_of_fits = np.swapaxes(final_list, 0, 1)

            outliers = True
            while outliers:
                outliers = False
                for i in range(len(list_of_fits)):
                    if list_of_fits[i][0].header[variables.total_samples.keyword] not in nsamps:
                        list_of_fits = np.delete(list_of_fits, i)
                        list_of_files = np.delete(list_of_files, i)
                        outliers = True
                        break

            self.file_names = list_of_files
            self.spectroscopic_images = list_of_fits
            self.direct_image = direct_image
            self.splitted = False
            self._data_set_directory_path = input_data

        elif os.path.isdir(input_data) and len(glob.glob(os.path.join(input_data, '*', ''))) > 0:

            self.file_names = []
            self.spectroscopic_images = []
            self.direct_image = []
            self.splitted = True
            self._data_set_directory_path = input_data

            for input_data in sorted(glob.glob(os.path.join(input_data, '*', ''))):

                nsamp = []
                final_list = []
                direct_image = None

                files = sorted(glob.glob(os.path.join(input_data, '*.fits')))
                for i in files:
                    with pf.open(i) as j:
                        if j[0].header[variables.observation_type.keyword] == 'SPECTROSCOPIC':
                            final_list.append([j[0].header[variables.exposure_start.keyword],
                                               os.path.split(i)[1], plc.copy_fits(j)])
                            nsamp.append(j[0].header[variables.total_samples.keyword])
                        elif not direct_image:
                            direct_image = pf.open(i, mode='update')

                nsamps = [int(np.median(np.array(nsamp)))]

                final_list.sort()
                list_of_times, list_of_files, list_of_fits = np.swapaxes(final_list, 0, 1)

                outliers = True
                while outliers:
                    outliers = False
                    for i in range(len(list_of_fits)):
                        if list_of_fits[i][0].header[variables.total_samples.keyword] not in nsamps:
                            list_of_fits = np.delete(list_of_fits, i)
                            list_of_files = np.delete(list_of_files, i)
                            outliers = True
                            break

                self.file_names = list_of_files
                self.spectroscopic_images.append(list_of_fits)
                self.direct_image = direct_image

        else:
            raise IraclisFileError('No such file or directory: {0}'.format(input_data))

        if not self.direct_image:
            raise IraclisFileError('A direct image is necessary.')

    def save(self, export_directory, arrange=True, export_pipeline_variables_file='variables.txt'):

        if os.path.isdir(export_directory):
            backup = '{0}_{1}'.format(export_directory, time.strftime('%y-%m-%d_%H-%M-%S'))
            shutil.copytree(export_directory, backup)
            shutil.rmtree(export_directory)

        os.mkdir(export_directory)

        if arrange:
            for i in range(len(self.file_names)):
                date = str(self.spectroscopic_images[i][0].header['DATE-OBS'])
                obs_time = str(self.spectroscopic_images[i][0].header['TIME-OBS'])
                obs_time = '-'.join(obs_time.split(':'))
                if self.file_names[i].split('_')[0] != date or self.file_names[i].split('_')[1] != obs_time:
                    self.file_names[i] = '{0}_{1}_{2}'.format(date, obs_time, os.path.split(self.file_names[i])[1])

        for i in range(len(self.file_names)):
            plc.copy_fits(self.spectroscopic_images[i]).writeto(
                os.path.join(export_directory, self.file_names[i]), output_verify='fix')

        plc.copy_fits(self.direct_image).writeto(
            os.path.join(export_directory, 'direct_image.fits'), output_verify='fix')

        if export_pipeline_variables_file:
            variables.save(os.path.join(export_directory, export_pipeline_variables_file))

    def copy_split(self, split_number):

        x = DataSet()
        x.spectroscopic_images = self.spectroscopic_images[split_number]

        return x


# pipeline counter


class PipelineCounter:

    def __init__(self, task, total_iterations, show_every=1):

        self.task = '{0}{1}'.format(task, '.' * (15 - len(task)))
        self.current_iteration = 0
        self.total_iterations = int(total_iterations)
        self.start_time = time.time()
        self.show = 0
        self.show_every = int(show_every)

        if self.total_iterations == 1:
            self.show_every = 10

    def update(self):

        self.current_iteration += 1
        self.show += 1.0 / self.show_every

        out_of = '{0}{1} / {2}'.format(' ' * (len(str(self.total_iterations)) - len(str(self.current_iteration))),
                                       str(self.current_iteration), str(self.total_iterations))

        delta_time = time.time() - self.start_time

        time_left = str(datetime.timedelta(
            seconds=int((self.total_iterations - self.current_iteration) * delta_time / self.current_iteration)))

        total_time = str(datetime.timedelta(seconds=int(delta_time)))

        if int(self.show):
            sys.stdout.write('\r\033[K')
            sys.stdout.write('{0}: {1}   time left: {2}   total time: {3}'.format(
                self.task, out_of, time_left, total_time))
            sys.stdout.flush()
            self.show = 0

        if self.current_iteration == self.total_iterations and self.total_iterations > 1:
            print('')
