# -*- coding: utf-8 -*-
"""
Created on Mon May 18 15:42:34 2015

@author: sprenger
"""

import quantities as pq

import developmentio as DIO

sessiondir = '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/data/2015-02-03_14-22-50'
mdatadir = '/mnt/Transcend/Datasets/SPP1665/Data/devel-circuits/metadata'
cachedir = '/home/julia/cache'

IO = DIO.DevelopmentIO(sessiondir=sessiondir,mdatadir=mdatadir,cachedir=cachedir,use_cache='datesize',print_diagnostic=True)

block = IO.read_block(t_starts=0*pq.s,t_stops=30*pq.s,electrode_list=[1],units=None)

print block

pass