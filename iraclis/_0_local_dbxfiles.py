from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import dropbox
import pickle
import os

package_name = 'iraclis'

os.chdir(os.path.abspath(os.path.dirname(__file__)))

dbx = dropbox.Dropbox('ThsvEzueGyAAAAAAAAAWqcWQ6EDURhvUMqJRm_1Y_fQZW5sy5Wv7ak2efTfBoGqS')

dbx_directories = [str(ff.name) for ff in dbx.files_list_folder('/Storage/{0}'.format(package_name)).entries]

dbx_files = {ff: {str(gg.name): {'path': '/Storage/{2}/{0}/{1}'.format(ff, gg.name, package_name), 'link': ''} for
                  gg in dbx.files_list_folder('/Storage/{1}/{0}'.format(ff, package_name)).entries} for
             ff in dbx_directories}

for i in dbx_files:
    for j in dbx_files[i]:
        dbx_files[i][j]['link'] = str(dbx.sharing_create_shared_link(dbx_files[i][j]['path']).url).replace('?dl=0',
                                                                                                           '?dl=1')
        print(i, j, dbx_files[i][j]['link'])

for i in dbx_files:
    for j in dbx_files[i]:
        dbx_files[i][j]['local_path'] = dbx_files[i][j]['path'].replace('/Storage/{0}/'.format(package_name), '')

pickle.dump(dbx_files, open('_0database.pickle', 'wb'))
