
import os
import ssl
import time
import shutil
import urllib

from urllib.request import urlretrieve
from pylightcurve.tools_files import open_dict, save_dict

version = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), '__version__.txt')).read()
build_in_databases_pickle_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), '__databases__.pickle')


ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def download(link, destination):
    try:
        urlretrieve(link, destination)
        return True
    except:
        try:
            with urllib.request.urlopen(link, context=ctx) as u, \
                    open(destination, 'wb') as f:
                f.write(u.read())

            return True
        except:
            print('Could not download {0}'.format(link))
            return False


class IraclisData:

    def __init__(self, _reset=False, _test=False):

        self.package_name = 'iraclis'

        self.databases_directory_path = os.path.join(os.path.abspath(os.path.expanduser('~')),
                                                     '.{0}'.format(self.package_name))

        self.databases_pickle_path = os.path.join(self.databases_directory_path, 'databases.pickle')

        # initiate databases

        if not os.path.isdir(self.databases_directory_path):
            os.mkdir(self.databases_directory_path)

        if not os.path.isfile(self.databases_pickle_path):
            shutil.copy(build_in_databases_pickle_path, self.databases_pickle_path)

        # check for updates in the databases (identified on github)

        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################

        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO

        # the link should be changed
        # and the line below deleted

        if not os.path.isfile(os.path.join(self.databases_directory_path, time.strftime('%y%m%d') + '.txt')):

            # download('https://github.comm/ucl-exoplanets/pylightcurve/blob/master/__databases__.pickle?raw=true',
            #          self.databases_pickle_path)

            shutil.copy(build_in_databases_pickle_path, self.databases_pickle_path)

        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO

        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################

        # load databases

        self.databases = open_dict(self.databases_pickle_path)

        self.hstwfc3_loaded = self._setup_database('hstwfc3')

        w = open(os.path.join(self.databases_directory_path, time.strftime('%y%m%d') + '.txt'), 'w')
        w.close()

    def hstwfc3(self):
        return self.hstwfc3_loaded

    def _setup_database(self, database_name):

        # define paths

        directory_path = os.path.join(self.databases_directory_path, database_name)
        pickle_path = os.path.join(self.databases_directory_path, database_name + '.pickle')
        pickle_path_new = os.path.join(self.databases_directory_path, database_name + '_new.pickle')
        last_update_file_path = os.path.join(self.databases_directory_path,
                                             '{0}_last_update.txt'.format(database_name))

        # define paths

        # check if folder exists, if not clean database

        clean_last_update = False

        if not os.path.isdir(directory_path):

            clean_last_update = True

            os.mkdir(directory_path)

            try:
                os.remove(pickle_path)
            except:
                pass

        if not os.path.isfile(pickle_path):
            if not download(self.databases[version][database_name], pickle_path):
                print('\n{0} features cannot be used.'.format(database_name))
                return False

        # check if folder exists

        # check for updates, remove files that need to be updated

        if not os.path.isfile(os.path.join(self.databases_directory_path, time.strftime('%y%m%d') + '.txt')):

            current_database = open_dict(pickle_path)

            if download(self.databases[version][database_name], pickle_path_new):

                new_database = open_dict(pickle_path_new)

                for dbx_file in new_database['files']:

                    if dbx_file in current_database['files']:
                        if new_database['files'][dbx_file]['link'] != current_database['files'][dbx_file]['link']:
                            try:
                                os.remove(os.path.join(self.databases_directory_path,
                                                       new_database['files'][dbx_file]['local_path']))
                            except:
                                pass

                shutil.move(pickle_path_new, pickle_path)

        # check for updates, remove files that need to be updated

        # download missing files

        print('Checking {0} database...'.format(database_name))

        final_check = True

        current_database = open_dict(pickle_path)
        dbx_files = current_database['files']
        frequency = current_database['frequency']

        for dbx_file in dbx_files:
            if not os.path.isfile(os.path.join(self.databases_directory_path,
                                               dbx_files[dbx_file]['local_path'])):
                clean_last_update = True
                print('\tDownloading: ', dbx_file)
                download(dbx_files[dbx_file]['link'], os.path.join(self.databases_directory_path,
                                                                   dbx_files[dbx_file]['local_path']))

                if not os.path.isfile(os.path.join(self.databases_directory_path,
                                                   dbx_files[dbx_file]['local_path'])):
                    final_check = False

        # download missing files

        # update files from external links

        if clean_last_update:
            try:
                os.remove(last_update_file_path)
            except:
                pass

        if frequency:

            try:
                last_update_date = int(open(last_update_file_path).read())
            except:
                last_update_date = 0

            must_update_date = int(time.strftime('%y%m%d')) - frequency

            if last_update_date < must_update_date:

                for dbx_file in dbx_files:
                    if 'external_link' in dbx_files[dbx_file]:
                        print('\tUpdating: ', dbx_file)
                        download(dbx_files[dbx_file]['external_link'],
                                 os.path.join(self.databases_directory_path, dbx_files[dbx_file]['local_path']))

                        if not os.path.isfile(os.path.join(self.databases_directory_path,
                                                           dbx_files[dbx_file]['local_path'])):
                            final_check = False

                w = open(last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

        # update files from external links

        if not final_check:
            print('\n{0} features cannot be used.'.format(database_name))
            return False
        else:
            return directory_path


iraclis_data = IraclisData()
