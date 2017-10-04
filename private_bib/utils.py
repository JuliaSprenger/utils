# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:59:01 2015

@author: julia sprenger
"""
import os
import re
import glob


#class utilities():

local = ''
# if os.path.isdir('/home/julia/projects'):
#     local = 'laptop'
# el
if os.path.isdir('/mnt/Transcend/Datasets'):
    local = 'disc'
elif os.path.isdir('/home/julia/'):
    local = 'laptop'
elif os.path.isdir('/users/sprenger/python/modules'):
    local = 'hambach'
elif os.path.isdir('/home/j.sprenger'):
    local = 'blaustein'

print('Working on data saved at {}'.format(local))

spp_struct = re.compile('.{4}-.{2}-.{2}_.{2}-.{2}-.{2}')

paths = {'laptop':{'hamburg':{'data':'/home/julia/data/SPP/data',
                              'metadata':'/home/julia/data/SPP/metadata',
                              'cache':'/home/julia/cache',
                              'sortdir':'/home/julia/projects/SPP1665/analysis/spikesortings'}, # This is not the best solution as there is also the hambach sorting locally available...
                   'marseille':{'data':'/home/julia/data/marseille/ReachGrasp',
                                'metadata':'/home/julia/data/marseille/ReachGrasp'},
                   'result':'/home/julia/projects/SPP1665/analysis'},
         'disc':{  'hamburg':{'data':'/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data',
                              'metadata':'/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/metadata',
                              'cache': '/home/julia/cache',
                              'sortdir':'/home/julia/projects/SPP1665/analysis/spikesortings'},  # This is not the best solution as there is also the hambach sorting locally available...
                   'marseille':{'data':'/mnt/Transcend/Datasets/marseille/ReachGrasp',
                                'metadata':'/mnt/Transcend/Datasets/marseille/ReachGrasp'},
                   'result':'/home/julia/projects/SPP1665/analysis'},
         'hambach':{'hamburg':{'data':'/datasets/devel-circuits/data/data',
                               'metadata':'/datasets/devel-circuits/data/metadata',
                               'cache': '~/cache'},
                    'marseille':{'data':'/datasets/marseille/congloue/data/DataGrasp',
                                 'metadata':'/datasets/marseille/congloue/data/DataGrasp'},
                    'result':'~/projects/SPP1665/analysis'}, #????
         'blaustein':{'hamburg':{'data':'/home/j.sprenger/data/SPP/data',
                                 'metadata':'/home/j.sprenger/data/SPP/metadata',
                                 'cache': '/home/j.sprenger/cache',
                                 'sortdir': '/home/j.sprenger/results/spikesortings'},
                      'marseille':{'data':'/home/j.sprenger/data/SPP/data',
                                   'metadata':'/home/j.sprenger/data/SPP/metadata'},
                      'result':'/home/j.sprenger/results',
                      }
         }


# datapath = {'hamburg':{'laptop':'/home/julia/data/SPP1665/Data/devel-circuits/data',
#                         'disc':'/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data',
#                         'hambach':'/datasets/devel-circuits/data/data',
#                         'blaustein':'/home/j.sprenger/data/SPP/data'},
#
#             'marseille':{'laptop':'/home/julia/data/marseille/ReachGrasp',
#                         'disc':'/mnt/Transcend/Datasets/marseille/ReachGrasp',
#                         'hambach':'/datasets/marseille/congloue/data/DataGrasp',
#                         'blaustein':'/home/j.sprenger/data/marseille/data'}}



monkey_translator = {'l':'Lilou','n':'Nikos','i':'Nikos2','t':'Tanya','s':'Sana'}

def get_data_location(filename=None):
    group = 'hamburg' if spp_struct.match(filename) else 'marseille'

    # walking into monkey directory
    refined_path = '/Data%s'%(monkey_translator[filename[0]]) if (group == 'marseille' and filename != None) else ''

    return paths[local][group]['data'] + refined_path


def get_metadata_location(filename):
    group = 'hamburg' if spp_struct.match(filename) else 'marseille'

    # walking into monkey meta data directory
    refined_path = '/MetaData%s'%(monkey_translator[filename[0]]) if (group == 'marseille' and filename != None) else ''

    return paths[local][group]['metadata'] + refined_path

def get_cache_location(filename):
    group = 'hamburg' if spp_struct.match(filename) else 'marseille'

    # walking into monkey meta data directory
    refined_path = '/MetaData%s'%(monkey_translator[filename[0]]) if (group == 'marseille' and filename != None) else ''

    if 'cache' not in paths[local][group]:
        raise ValueError('No cache directory saved for %s %s'%(local,group))
    return paths[local][group]['cache'] + refined_path

def get_result_location():
    return paths[local]['result']

def get_cache_location(filename):
    group = 'hamburg' if spp_struct.match(filename) else 'marseille'
    return paths[local][group]['cache']

def get_sort_location(filename):
    group = 'hamburg' if spp_struct.match(filename) else 'marseille'
    return paths[local][group]['sortdir']

def get_STA_predata_location(filename):
    if local:
        return "/home/julia/projects/Synchrofacts/synchrofact_manuscript/predata/%s"%(filename)
    else:
        return "/python/projects/synchrofacts/predata/%s"%(filename)
    

def get_unsrt_unit255_sessions():
    return ['n130515-003', 'n130726-002', 'n130510-002', 'n130726-003', 'n130509-003',
            't120310-002', 't290310-003', 't120310-001', 's131209-004', 's131209-003']





def get_ns6_noise_sessions(monkey):
    if monkey == 'Sana':
        path = "/home/julia/data/marseille/NoiseDataprelim/2016_01_14_Sana"
    elif monkey == 'Saline':
        path = "/home/julia/data/marseille/NoiseDataprelim/2016_01_15_Saline"

    ns6_files = glob.glob(path + '/*.ns6')
    sessions = [name.rstrip('.ns6') for name in ns6_files]
    return sessions



