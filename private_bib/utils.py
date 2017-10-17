# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:59:01 2015

@author: julia sprenger
"""
import os
import re
import glob


#class utilities():

# local = ''
# # if os.path.isdir('/home/julia/projects'):
# #     local = 'laptop'
# # el
# if os.path.isdir('/mnt/Transcend/Datasets'):
#     local = 'disc'
# elif os.path.isdir('/home/julia/'):
#     local = 'laptop'
# elif os.path.isdir('/users/sprenger/python/modules'):
#     local = 'hambach'
# elif os.path.isdir('/home/j.sprenger'):
#     local = 'blaustein'
#
# print('Working on data saved at {}'.format(local))

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

locs = list(paths)
# print('Locations {}'.format(locs))

# datapath = {'hamburg':{'laptop':'/home/julia/data/SPP1665/Data/devel-circuits/data',
#                         'disc':'/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data',
#                         'hambach':'/datasets/devel-circuits/data/data',
#                         'blaustein':'/home/j.sprenger/data/SPP/data'},
#
#             'marseille':{'laptop':'/home/julia/data/marseille/ReachGrasp',
#                         'disc':'/mnt/Transcend/Datasets/marseille/ReachGrasp',
#                         'hambach':'/datasets/marseille/congloue/data/DataGrasp',
#                         'blaustein':'/home/j.sprenger/data/marseille/data'}}


monkey_translator = {'l': 'Lilou', 'n': 'Nikos', 'i': 'Nikos2', 't': 'Tanya',
                     's': 'Sana', 'a': 'Tanya2', 'e': 'Enya'}

def _find_loc(session, group, type=None):
    loc = None
    for id in range(len(locs)):
        loc = locs[id]
        if type is None:
            path = paths[loc][group]
        else:
            path = paths[loc][group][type]
        if (glob.glob(os.path.join(path, '{}*'.format(session)))
                + glob.glob(os.path.join(path, '**', '{}*'.format(session)))):
            break
        loc = None

    if loc is None:
        raise ValueError('Can not determine location of {} {} for '
                         'session {}.'.format(group, type, session))
    return loc

def get_data_location(session=None):
    group = 'hamburg' if spp_struct.match(session) else 'marseille'

    # walking into monkey directory
    refined_path = '/Data%s'%(monkey_translator[session[0]]) if (group == 'marseille' and session != None) else ''

    loc = _find_loc(session=session, group=group, type='data')
    return paths[loc][group]['data'] + refined_path


def get_metadata_location(session):
    group = 'hamburg' if spp_struct.match(session) else 'marseille'

    loc = _find_loc(session=session, group=group, type='metadata')

    # walking into monkey meta data directory
    refined_path = '/MetaData%s'%(monkey_translator[session[0]]) if (group == 'marseille' and session != None) else ''

    return paths[loc][group]['metadata'] + refined_path

def get_cache_location(session):
    group = 'hamburg' if spp_struct.match(session) else 'marseille'

    # walking into monkey meta data directory
    refined_path = '/MetaData%s'%(monkey_translator[session[0]]) if (group == 'marseille' and session != None) else ''

    # using same location for caching as for data
    loc = _find_loc(session=session, group=group, type='data')

    if 'cache' not in paths[loc][group]:
        raise ValueError('No cache directory saved for %s %s'%(loc, group))
    return paths[loc][group]['cache'] + refined_path

def get_result_location():
    loc = _find_loc(group='result')
    return paths[loc]['result']

def get_sort_location(session):
    group = 'hamburg' if spp_struct.match(session) else 'marseille'
    loc = _find_loc(group=group, type='sortdir')
    return paths[loc][group]['sortdir']

# def get_STA_predata_location(session):
#     if local:
#         return "/home/julia/projects/Synchrofacts/synchrofact_manuscript/predata/%s"%(session)
#     else:
#         return "/python/projects/synchrofacts/predata/%s"%(session)
    

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



