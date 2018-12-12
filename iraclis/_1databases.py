from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._0errors import *
from ._0imports import *


class Database:

    def __init__(self, database_name, vital=False, force_update=False):

        self.database_name = database_name
        self.database_info_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), '_0database.pickle')
        self.package_database_location = os.path.join(os.path.expanduser('~'), '.iraclis')
        if not os.path.isdir(self.package_database_location):
            os.mkdir(self.package_database_location)

        self.database_location = os.path.join(self.package_database_location, '{0}_database'.format(self.database_name))
        self.database_last_update_file_path = os.path.join(self.package_database_location,
                                                           '{0}_database_last_update.txt'.format(self.database_name))

        if os.path.isdir(self.database_location):
            if force_update or len(glob.glob(os.path.join(self.database_location, '*'))) == 0:
                shutil.rmtree(self.database_location)
                os.mkdir(self.database_location)
                update = True
            else:
                if not os.path.isfile(self.database_last_update_file_path):
                    update = True
                elif int(open(self.database_last_update_file_path).readlines()[0]) < 181212:
                    update = True
                else:
                    update = False
        else:
            os.mkdir(self.database_location)
            update = True

        if update:
            # noinspection PyBroadException
            try:
                print('\nDownloading {0} database...'.format(self.database_name))

                dbx_files = pickle.load(open(self.database_info_file_path, 'rb'))
                dbx_files = dbx_files['{0}_database'.format(self.database_name)]

                for i in glob.glob(os.path.join(self.database_location, '*')):
                    if os.path.split(i)[1] not in dbx_files:
                        os.remove(i)

                for i in dbx_files:
                    if not os.path.isfile(os.path.join(self.package_database_location, dbx_files[i]['local_path'])):
                        print(i)
                        urlretrieve(dbx_files[i]['link'],
                                    os.path.join(self.package_database_location, dbx_files[i]['local_path']))

                w = open(self.database_last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

            except:
                print('\nDownloading {0} database failed. A download will be attempted next time.'.format(
                    database_name))
                pass

        if (not os.path.isdir(self.database_location) or
                len(glob.glob(os.path.join(self.database_location, '*'))) == 0):
            if vital:
                raise IraclisLibraryError('{0} database not available.'.format(database_name))
            else:
                print('\n{0} features cannot be used.'.format(self.database_name))
                self.database_location = False


class Databases:

    def __init__(self):

        self.wfc3 = Database('wfc3', vital=True).database_location


databases = Databases()
