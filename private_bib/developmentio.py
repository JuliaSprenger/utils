# coding=utf-8
'''
Causemech IO functions

This module can be used to load data from the optogenetic Causative Mechanism
Experiment via the Neuralynx loading IO of neo.

Authors: Julia Sprenger
Based on the reach-grasp IO
'''

import copy
import os

import numpy as np
import odml.tools
import quantities as pq
import warnings

import neo

from neo.core import ChannelIndex



class DevelopmentIO(neo.io.NeuralynxIO):
    """
    Derived class to handle a managing a Neuralynx recording session from the
    optogenetic causative mechanisms experiment.

    """

    # put global definitions / translation tables here



    def __init__(self, sessiondir, mdatadir=None, cachedir=None,use_cache='hash',
                 print_diagnostic=False):
        """
        Constructor

        Arguments:
            sessiondir : the directory the files of the recording session are
                            collected. Default 'None'.
            print_diagnostic: indicates, whether information about the loading of
                            data is printed in terminal or not. Default 'False'.
            check_files: check associated files on consistency. This is highly
                            recommended to ensure correct operation. Default 'True'
        """

        # Remember choice whether to print diagnostic messages or not
        self._print_diagnostic = print_diagnostic

        # Initialize session
        neo.io.NeuralynxIO.__init__(self, sessiondir=sessiondir,cachedir=cachedir,use_cache=use_cache,
                                print_diagnostic=self._print_diagnostic)
        if mdatadir == None:
            mdatadir = sessiondir

        # Save location of metadata directory
        self.metadata_dir = mdatadir

        # Get odML_doc if available
        self.odML_avail = False
        self.odML_filename = None
        self.odML_doc = None

        session = os.path.basename(sessiondir)


        # Looking for odml file
        if os.path.isdir(os.path.join(self.metadata_dir ,'odml_builds')):
            if os.path.isfile(os.path.join(self.metadata_dir ,'odml_builds','%s.odml'%session)):
                self.odML_avail = True
                self.odML_filename = session
            elif os.path.isfile(os.path.join(self.metadata_dir ,'odml_builds','%s_prefilled.odml'%session)):
                self.odML_avail = True
                self.odML_filename = session + '_prefilled'

        if self.odML_filename != None:
            # reading odml file
            self.odML_doc = odml.tools.xmlparser.load(os.path.join(self.metadata_dir,
                                                               'odml_builds',
                                                               '%s.odml'%self.odML_filename))
            #check if this odml actually belongs to the recording session
            odml_rec_file = self.odML_doc.get_property_by_path(
                            'Experiment/Recording:RecordingSession').value.data
            if  odml_rec_file != os.path.basename(os.path.normpath(session)):
               raise ValueError('Associated Odml file belongs to different'
                               ' recording session (%s)'%odml_rec_file)




        if self.odML_filename != None:
            self._diagnostic_print('Using odml file "%s"'%(self.odML_filename))



    def read_block(self, lazy=False, cascade=True, t_starts=[None], t_stops=[None],
                    electrode_list=[], unit_list=None, analogsignals=True,
                   events=False,
                    waveforms = False):
        """
        Reads data in a requested time window and returns block with single segment
        containing these data.

        Arguments:
            lazy : Postpone actual reading of the data files. Default 'False'.
            cascade : Do not postpone reading subsequent neo types (segments).
                            Default 'True'.
            t_starts : list of quantities or quantity describing the start of the
                            requested time window to load. If None or [None]
                            the complete session is loaded. Default '[None]'.
            t_stops : list of quantities or quantity describing the end of the
                            requested time window to load. Has to contain the
                            same number of values as t_starts. If None or [None]
                            the complete session is loaded. Default '[None]'.
            channel_list : list of integers containing the IDs of the requested
                            to load. If [] all available channels will be loaded.
                            Default: [].
            unit_list : list of integers containing the IDs of the requested
                            units to load. If [] all available units will be
                            loaded.
                            Default: None.
            events : Loading events. If True all available events in the given
                            time window will be read. Default: False.
            load_waveforms : Load waveform for spikes in the requested time
                            window. Default: False.

        Returns: Block object containing the requested data in neo structures.

        Usage:
            from neo import io
            import quantities as pq
            import matplotlib.pyplot as plt

            session_folder = '../Data/2014-07-24_10-31-02'
            NIO = io.NeuralynxIO(session_folder,print_diagnostic = True)
            block = NIO.read_block(lazy = False, cascade = True,
                                   t_starts = 0.1*pq.s, t_stops = 0.2*pq.s,
                                   channel_list = [1,5,10], unit_list = [1,2,3],
                                   events = True, load_waveforms = True)
        """



        # Load neo block
        block = neo.io.NeuralynxIO.read_block(self, lazy=lazy, cascade=cascade,
                                       t_starts=t_starts, t_stops=t_stops,
                                       electrode_list=electrode_list,
                                       unit_list=unit_list,
                                       analogsignals=analogsignals,
                                       events=events, waveforms=waveforms)



        # TODO: odML <-> data files consistency checks? Low priority

        # Add annotations of odML meta data info
        if self.odML_avail:

            """
            TODO:
            * Add electroporation area to recording channel group
            """


            area_dict = self.get_electrodes_by_area()
            electroporated_areas = self.get_electroporation()
            channel_indexes = [r for r in block.channel_indexes
                               if r.name == 'all channels']

            for area, channels in area_dict.iteritems():
                electroporated, expression = False, None
                if area in electroporated_areas.keys():
                    electroporated = True
                    expression = electroporated_areas[area]
                    chidx = ChannelIndex(name='%s channels'%area,
                                   channel_indexes=channels,
                                   channel_names=['channel %i'%i for i in channels],
                                   electroporated=electroporated,
                                   expression=expression)
                    chidx.block = block
                    block.channel_indexes.append(chidx)

            # raise NotImplementedError('neo block annotation using odmls is not implemented yet.')



            # ########### Annotate information of 'Recording' Section ############
            # # Annotate Amplifier Information
            # amp_properties = ['LowpassCutoff','HighpassCutoff','SamplingRate']
            # ff = lambda x: x.name in amp_properties and 'Amplifier' in x.parent.get_path()
            # pobj = {p.name:p.value.data for p in self.odML_doc.iterproperties(filter_func=ff)}
            # block.annotate(amplifier= pobj)
            #
            # # Consistency Check with Analogsignal Sampling Rate
            # if any([pobj['SamplingRate'] != asa.annotations['SamplingFrequency']
            #             for asa in block.segments[0].analogsignalarrays]):
            #     raise ValueError('Inconsistent sampling rates detected in odml'
            #                 ' and original data files (%s / %s)'%(
            #                 pobj['SamplingRate'],
            #                 [asa.annotations['SamplingFrequency'] for asa in
            #                  block.segments[0].analogsignalarray]))
            #
            # # Annotate different Recording Areas
            # # Extracting Recording Area sections
            # ff = lambda x: 'RecordingArea' in x.name and 'Probe' in x.sections
            # recording_secs = [p for p in self.odML_doc.itersections(filter_func=ff)]
            # rec_properties = ['Hemisphere','Probe ID','Channels',
            #                   'SpikingChannels','BrokenChannels','Quality']
            # ff2 = lambda x: x.name in rec_properties
            # area_dict = {}
            # for recording_sec in recording_secs:
            #     # extracting specific properties of each recording area section
            #     area_dict[recording_sec.name] = {a.name:a.value.data for a in
            #                     recording_sec.iterproperties(filter_func=ff2)}
            #     # adding two 'area' properties manually as they have the same name
            #     area_dict[recording_sec.name]['RecordingArea'] = \
            #         recording_sec.properties['Area'].value.data
            #     area_dict[recording_sec.name]['ReferenceArea'] = \
            #         recording_sec.get_property_by_path('Reference:Area').value.data
            # block.annotate(recordingareas=area_dict)

        return block


    def read_segment(self,lazy=False, cascade=True, t_start=None, t_stop=None,
                        electrode_list=[], unit_list=None, analogsignals=True,
                        events=False, waveforms=False):
        """Reads one Segment.

        The Segment will contain one AnalogSignalArray for each channel
        and will go from t_start to t_stop.

        Arguments:


            lazy : Postpone actual reading of the data files. Default 'False'.
            cascade : Do not postpone reading subsequent neo types (SpikeTrains,
                            AnalogSignalArrays, Events).
                            Default 'True'.
            t_start : time (quantity) that the Segment begins. Default None.
            t_stop : time (quantity) that the Segment ends. Default None.
            electrode_list : list of integers containing the IDs of the requested
                            to load. If [] all available channels will be loaded.
                            Default: [].
            unit_list : list of integers containing the IDs of the requested
                            units to load. If [] all available units will be
                            loaded.
                            Default: None.
            analogsignals : boolean, indication whether analogsignals should be
                            read. Default: True.
            events : Loading events. If True all available events in the given
                            time window will be read. Default: False.
            waveforms : Load waveform for spikes in the requested time
                            window. Default: False.


        Returns:
            Segment object containing neo objects, which contain the data.
        """

        # Load neo segment
        seg = neo.io.NeuralynxIO.read_segment(self,lazy=lazy,
                                              cascade=cascade,
                                              t_start=t_start,
                                              t_stop=t_stop,
                                              electrode_list=electrode_list,
                                              unit_list=unit_list,
                                              analogsignals=analogsignals,
                                              events=events,
                                              waveforms=waveforms)

        # # Generate t_start and t_stop annotations of segments
        # seg.annotations['t_start'] = min([a.t_start for a in seg.analogsignalarrays + seg.spiketrains])
        # seg.annotations['t_stop'] = max([a.t_stop for a in seg.analogsignalarrays + seg.spiketrains])

        signaltype = []
        # Generate analogsignal classifications
        for sig in seg.analogsignals:
            for channel_idx in sig.annotations['channel_index']:
                if channel_idx<=32:
                    signaltype.append('neural')
                elif channel_idx in [32,35]:
                    signaltype.append('stimulation')
                else:
                    raise TypeError('Signal has unkown channel type (id %s)'%channel_idx)
            sig.annotations['signal_type'] = signaltype

        for sig in seg.spiketrains:
            if 'electrode_id' in sig.annotations and sig.annotations['electrode_id']<=32:
                sig.annotations['signaltype'] = 'neural'
            elif 'electrode_id' in sig.annotations and sig.annotations['electrode_id'] in [32,35]:
                sig.annotations['signaltype'] = 'stimulation'
            elif 'electrode_id' in sig.annotations:
                raise TypeError('Signal has unkown electrode_id annotation (%s)'%sig.annotations['electrode_id'])
            else:
                raise TypeError('Signal has no electrode_id annotation')

        if self.odML_avail:

            ################ STIMULATIONS ######################
            # Get stimulation periods from odml
            stimulations, stimarea = self.get_stimulations()

            if stimulations != None:
                # Add stimulations as epocharray
                s_start = seg.t_start if seg.t_start else 0*pq.s
                s_stop = seg.t_stop if seg.t_stop else self.parameters_global['t_stop'] - self.parameters_global['t_start']

                # Add only stimulations which are completely in the current segment
                stim = {k:v for k,v in stimulations.iteritems() if v['StimulationPeriodEnd']-self.parameters_global['t_start'] >s_start and v['StimulationPeriodStart']-self.parameters_global['t_start']<s_stop}

                start_times = [t['StimulationPeriodStart']-self.parameters_global['t_start'] for t in stim.itervalues()]
                durations = [t['StimulationPeriodEnd'] - t['StimulationPeriodStart'] for t in stim.itervalues()]

                if len(durations)>0 and len(start_times)>0:
                    start_times = pq.Quantity([t.rescale(start_times[0].units) for t in start_times],start_times[0].units)
                    durations = pq.Quantity([d.rescale(durations[0].units) for d in durations],durations[0].units)

                    labels = ['Stimulationperiod ' + str(i) for i in stim.iterkeys()]
                    stimtype = [t['StimulationType'] for t in stim.itervalues()]
                    stimfreq = [t['StimulusFrequency'] if 'StimulusFrequency' in t else None for t in stim.itervalues()]
                    stimpulseduration = [t['PulseDuration'] if 'PulseDuration' in t else None for t in stim.itervalues()]
                    stim_count = [t['N_Stimuli'] for t in stim.itervalues()]
                    stimoutput = [t['LaserOutput'] for t in stim.itervalues()]
                    stimpower = [t['LaserPower'] for t in stim.itervalues()]
                    stimquality = [t['StimulationQuality'] for t in stim.itervalues()]
                    # stimtimeunits = [t['StimulationTimes'][0].units if type(t['StimulationTimes'])==list else t['StimulationTimes'].units for t in stim.itervalues()]
                    stimtimeunit = stim.values()[0]['StimulationTimes'][0].units

                    stimtimes = [np.asarray(t['StimulationTimes'])*stimtimeunit - self.parameters_global['t_start'] for t in stim.itervalues()]


                    ep = neo.Epoch(times=start_times,
                                        durations=durations,
                                        labels=labels,
                                        name='Stimulation epoch',
                                        file_origin=self.odML_filename + '.odml',
                                        type="stimulation",
                                        stimtype=stimtype,
                                        definition='Times of optogenetic stimulation',
                                        stimarea = stimarea,
                                        stimfreq=stimfreq,
                                        stimpulseduration=stimpulseduration,
                                        stim_count=stim_count,
                                        stimpower=stimpower,
                                        stimoutput=stimoutput,
                                        stimquality=stimquality,
                                        stimtimes=stimtimes)

                    seg.epochs.append(ep)
                    seg.create_relationship()


            ################ SPINDLES ######################
            # Get stimulation periods from odml
            spindledet = self.get_spindles()

            # Add spindles as epocharray
            s_start = seg.t_start
            s_stop = seg.t_stop
            # Add only spindles which are contained in the current segment # if (spindle['SpindleEnd']-self.parameters_global['t_start']>s_start) and (spindle['SpindleStart']-self.parameters_global['t_start']<s_stop)
            # stim = {p:{c:{s:spindle for s,spindle in channel.iteritems() if (spindle['SpindleEnd']-self.parameters_global['t_start']>s_start) and (spindle['SpindleStart']-self.parameters_global['t_start']<s_stop)} for c,channel in period.iteritems()} for p,period in spindledet.iteritems()}

            # convert spindle to
            for p,period in spindledet.iteritems():
                for c,channel in period.iteritems():
                    valid_spindles = {}
                    for s, spindle in channel.iteritems():
                        if type(s)==int and (spindle['SpindleEnd']-self.parameters_global['t_start']>s_start) and (spindle['SpindleStart']-self.parameters_global['t_start']<s_stop):
                            valid_spindles.update({s:spindle})

                    if len(valid_spindles.keys())>0:


                        start_times = [t['SpindleStart']-self.parameters_global['t_start'] for t in valid_spindles.itervalues()]
                        start_times = pq.Quantity([t.rescale(start_times[0].units) for t in start_times],start_times[0].units)
                        durations = [t['SpindleEnd'] - t['SpindleStart'] for t in valid_spindles.itervalues()]
                        durations = pq.Quantity([d.rescale(durations[0].units) for d in durations],durations[0].units)
                        
                        labels = ['Spindleperiod %i in channel %i in time period %s'%(i,c,p) for i in valid_spindles.iterkeys()]
                        spindamplitude = [t['SpindleAmplitude'] for t in valid_spindles.itervalues()]
                        spindmaxamplitude = [t['SpindleMaxAmplitude'] for t in valid_spindles.itervalues()]

                        ep = neo.Epoch(times=start_times,
                                            durations=durations,
                                            labels=labels,
                                            name='Spindle epoch',
                                            file_origin=self.odML_filename + '.odml',
                                            type="spindle oscillation",
                                            channel_id = c,
                                            spindperiod = p,
                                            definition='Times of LFP spindle oscillations',
                                            spindamplitude=spindamplitude,
                                            spindmaxamplitude=spindmaxamplitude)

                        seg.epochs.append(ep)
                        seg.create_relationship()
        return seg


############### Supplementory functions for data extraction from odml #################
    def get_stimulations(self):
        ff = lambda x:  x.name.startswith('Stimulation_') and \
                        (x.parent.name == 'Stimulations')
        sobjs = [s for s in self.odML_doc.itersections(filter_func=ff)]

        if len(sobjs)==0:
            return None, None

        sdict = {}
        for sobj in sobjs:
            temp = {}
            t = self._extract_quantities(sobj,['StimulationPeriodStart',
                                               'StimulationPeriodEnd',
                                               'StimulationDuration',
                                               'StimulusFrequency',
                                               'PulseDuration',
                                               'LaserPower',
                                               'StimulationTimes'])
            temp.update(t)
            for prop_name in ['N_Stimuli','StimulationType','LaserOutput','StimulationQuality']:
                prop = sobj.properties[prop_name]
                temp[prop_name] = prop.value.data

            stim_id = int(sobj.name.split('_')[-1])

            sdict[stim_id] = temp.copy()

        return sdict, sobjs[0].parent.properties['StimulationArea'].value.data


    def get_spindles(self):
        ff = lambda x:  x.name.startswith('SpindleDetection_') and \
                        (x.parent.name == 'Preprocessing')
        sobjs = [s for s in self.odML_doc.itersections(filter_func=ff)]

        sdict = {}
        for sobj in sobjs:
            periodname = sobj.properties['DetectionPeriodName'].value.data
            sdict[periodname] = {}
            ff = lambda x:  x.name.startswith('Channel_')
            cobjs = [c for c in sobj.itersections(filter_func=ff)]
            for cobj in cobjs:
                chid = cobj.properties['Channel_ID'].value.data
                sdict[periodname][chid] = {}
                t = self._extract_quantities(cobj,['SpindleRate',
                                                   'MeanSpindleAmplitude','MeanSpindleAmplitudeSD','MeanSpindleAmplitudeSEM',
                                                   'MeanSpindleMaxAmplitude','MeanSpindleMaxAmplitudeSD','MeanSpindleMaxAmplitudeSEM',
                                                   'MeanSpindleDuration','MeanSpindleDurationSD','MeanSpindleDurationSEM'])
                sdict[periodname][chid].update(t)

                ff = lambda x: x.name.startswith('Spindle_')
                spindobjs = [s for s in cobj.itersections(filter_func=ff)]
                for spindobj in spindobjs:
                    sid = int(spindobj.name.strip('Spindle_'))
                    sdict[periodname][chid][sid] = {}
                    t = self._extract_quantities(spindobj,['SpindleStart','SpindleEnd','SpindleDuration','SpindleAmplitude','SpindleMaxAmplitude'])
                    sdict[periodname][chid][sid].update(t)

        return sdict

    def get_electrodes_by_area(self):
        ff = lambda x: x.name.startswith('RecordingArea') and \
                'Area' in x.properties
        sobjs = [s for s in self.odML_doc.itersections(filter_func=ff)]

        result = {}
        for sobj in sobjs:
            area = sobj.properties['Area'].value.data
            channels = [v.data for v in sobj['Probe'].properties['ChannelIDs'].values]
            result[area] = channels

        return result


    def get_electroporation(self):
        ff = lambda x: x.name == 'InUteroElectroporation' and \
                        x.parent.name == 'Modification'
        sobjs = [s for s in self.odML_doc.itersections(filter_func=ff)]

        result = {}
        for sobj in sobjs:
            result[sobj.properties['TargetArea'].value.data] = sobj.properties['Expression'].value.data

        return result


    def _diagnostic_print(self, text):
        '''
        Print a diagnostic message.

        Args:
            text (string):
                Diagnostic text to print.

        Returns:
            -
        '''

        if self._print_diagnostic:
            print('DevelopmentIO: ' + text)


########### Supplementory functions ######################################
    def _extract_quantities(self,sobj,qlist):
        temp = {}
        for quant_prop_name in qlist:
            if quant_prop_name in sobj.properties:
                quant_prop = sobj.properties[quant_prop_name]
                quant_vals = quant_prop.values
                if len(quant_vals) == 1:
                    temp[quant_prop_name] = pq.Quantity(quant_vals[0].data,
                                                        units=quant_vals[0].unit,
                                                        dtype=quant_vals[0].dtype)
                else:
                    temp[quant_prop_name] = [pq.Quantity(val.data,
                                                     units=val.unit,
                                                     dtype=val.dtype)
                                             for val in quant_vals]
                #
        return temp

########### Utilities for neo objects ####################################
def slice_analogsignal_by_epochs(anasig,epoch,pre=0*pq.s,post=0*pq.s,invert=False):
    data = []
    t_starts = epoch.times - pre
    t_stops = epoch.times + epoch.durations + post

    # inverting time periods by appending t_start and t_stop of anasig and swapping t_starts and t_stops
    if invert:
        t_sta=t_starts.copy()
        t_sto=t_stops.copy()
        #testing units
        anasig.t_start.rescale(t_sta.units)
        #appending t_start
        t_starts = pq.Quantity(np.append(anasig.t_start.rescale(t_sto.units),t_sto).magnitude,t_sto.units)
        #testing unit
        anasig.t_stop.rescale(t_sto.units)
        #appending t_stop
        t_stops = pq.Quantity(np.append(t_sta,anasig.t_stop.rescale(t_sta.units)).magnitude,t_sta.units)

    for i in range(len(t_starts)):
        data.append(anasig.time_slice(t_starts[i],t_stops[i]))
    return data



########### Synchrofact Handling ##########################################

def keep_synchrofacts(block, n=2, dt=0, dt2=1, unit_type=['sua']):
    """
    Keeps only spike artefacts in spiketrains of given unit type in given block.
    Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
    unit_type = [] selects all available spike trains
    """
    del_spikes(block, n=n, dt=dt, dt2=dt2, invert=False, unit_type=unit_type)


def remove_synchrofacts(block, n=2, dt=0, dt2=1, unit_type=['sua']):
    """
    Removes all spike artefacts in spiketrains of block.
    Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
    unit_type = [] selects all available spike trains
    """
    del_spikes(block, n=n, dt=dt, dt2=dt2, invert=True, unit_type=unit_type)


# _____________ Helper Functions for synchrofact handling ________________
def del_spikes(block, n, dt, dt2, invert=True, unit_type=['sua']):
    '''
    Given block with spike trains, delete all spikes engaged
    in synchronous events of size *n* or higher. If specified, delete
    spikes close to such syncrhonous events as well.

    *Args*
    ------
    block [list]:
        a block containing neo spike trains

    n [int]:
        minimum number of coincident spikes to report synchrony

    dt [int. Default: 0]:
        size of timestep lag for synchrony. Spikes closer than *dt* are
        considered synchronous. Groups of *n* or more synchronous spikes
        are deleted from the spike trains.

    dt2 [int. Default: 1]:
        maximum distance of timesteps for two spikes to be "close". Spikes
        "close" to synchronous spikes are eliminated as well.

    invert [bool. Default: True]:
        selects synchronous events to be deleted (Default:True). False keeps only
        synchrofacts and definded neigboring spikes in spiketrains

    unit_type [list of strings. Default ['sua']]:
        selects only spiketrain of certain units / channels for synchrofact extraction.
        unit_type = [] considers all provided spiketrains
        Accepted unit types: 'sua','mua','idx' (whereas x is the id number requested)
        Warning: for invert=True non-selected spiketrains are emptied, for invert=False
                no change occurs in non-selected spiketrains
    '''

    # data check
    if len(block.segments[0].spiketrains) == 0:
        warnings.warn(
            'Attempted synchronous event extraction with empty block')
        return block

    dt = dt * block.segments[0].spiketrains[0].units
    dt2 = dt2 * block.segments[0].spiketrains[0].units

    # extracting spiketrains which should be used for synchrofact extraction based on given unit type
    # possible improvement by using masks for different conditions and adding
    # them up
    neo_spiketrains, index = [], []
    for i in range(len(block.segments[0].spiketrains)):
        take_it = False
        for utype in unit_type:
            if utype[:2] == 'id' and block.segments[0].spiketrains[i].annotations['unit_id'] == int(utype[2:]):
                take_it = True
            elif (utype == 'sua' or utype == 'mua') and utype in block.segments[0].spiketrains[i].annotations and block.segments[0].spiketrains[i].annotations[utype]:
                take_it = True
        if take_it:
            neo_spiketrains.append(
                copy.deepcopy(block.segments[0].spiketrains[i]))
            index.append(i)

    # considering all spiketrains for unit_type == []
    if unit_type == []:
        neo_spiketrains = copy.deepcopy(block.segments[0].spiketrains)
        index = range(len(block.segments[0].spiketrains))

    # if no spiketrains were selected
    if len(neo_spiketrains) == 0:
        warnings.warn(
            'No matching spike trains for given unit selection criteria %s found' % (unit_type))
        times, dt = np.array([]), 0 * pq.s
        time_dt2 = np.array([])
    else:
        # find times of synchrony of size >=n
        times, dt = detect_syn_spikes(neo_spiketrains, n=n, dt=dt)
        # hstack only works for input != []
        times_dt2 = np.hstack([np.arange(
            j - dt2.base, j + dt2.base + 1) for j in times.magnitude]) if times != [] else []

    j = 0  # index of pre selected sts
    # iterating over original spike trains
    for idx in range(len(block.segments[0].spiketrains)):
        # for j, idx in enumerate(index):  # and copy in it the original ones
        # devoided of the synchrony times
        if idx in index:

            annotations = block.segments[0].spiketrains[idx].annotations

            mask = np.in1d(neo_spiketrains[j], times_dt2)
            if invert:
                mask = np.invert(mask)

            block.segments[0].spiketrains[
                idx] = neo_spiketrains[j][mask.nonzero()]

            if neo_spiketrains[j].waveforms != None:
                block.segments[0].spiketrains[idx].waveforms = neo_spiketrains[
                    j].waveforms[mask.nonzero()]

            block.segments[0].spiketrains[idx].annotations = annotations
            j += 1
        # case st was not selected -> removal of spikes in case of 'keep
        # synchrofacts'
        else:
            if not invert:
                block.segments[0].spiketrains[idx] = neo.SpikeTrain(np.array([]) * pq.ms,
                                                        t_start=block.segments[
                                                            0].spiketrains[idx].t_start,
                    t_stop=block.segments[
                                                            0].spiketrains[idx].t_stop,
                    **block.segments[0].spiketrains[idx].annotations)


# self, x, t_stop=10, t_start=0):
def spiketrains2gdf(neo_spiketrains, ids=[]):
    """
    Converts a list of spike trains to gdf format.

    Gdf is a 2-column data structure containing neuron ids on the first
    column and spike times (sorted in increasing order) on the second
    column. Information about the time unit, not preserved in the float-
    like gdf, is returned as a second output

    *Args*
    ------
    sts [list]:
        a list of neo spike trains.

    ids [list. Default to []]:
        List of neuron IDs. Id[i] is the id associated to spike train
        sts[i]. If empty list provided (default), dds are assigned as
        integers from 0 to n_spiketrains-1.

    *Returns*
    ---------
    gdf [ndarray of floats with shape (n_spikes, 2)]:
        ndarray of unit ids (first column) and
    """

#        t_stop = max([train.t_stop for train in neo_spiketrains])
#        t_start = min([train.t_start for train in neo_spiketrains])

    # Find smallest time unit
    time_unit = neo_spiketrains[0].units
    for st in neo_spiketrains[1:]:
        if st.units < time_unit:
            time_unit = st.units

#        self.unit = time_unit
#        self.times = neo.SpikeTrain(self[:, 1] * pq.second, t_start=t_start, t_stop=t_stop)
#        self.dtypes = (int, neo.SpikeTrain)

    # By default assign integers 0,1,... as ids of sts[0],sts[1],...
    if len(ids) == 0:
        ids = range(len(neo_spiketrains))

    gdf = np.zeros((1, 2))
    # Rescale all spike trains to that time unit, and add to the gdf
    for st_idx, st in zip(ids, neo_spiketrains):
        to_be_added = np.array([[st_idx] * len(st),
            st.view(pq.Quantity).rescale(time_unit).magnitude]).T
        gdf = np.vstack([gdf, to_be_added])

    # Eliminate first row in gdf and sort the others by increasing spike
    # times
    gdf = gdf[1:]
    gdf = gdf[np.argsort(gdf[:, 1])]

    # Return gdf and time unit corresponding to second column
    return gdf, time_unit


def detect_syn_spikes(neo_spiketrains, n=2, dt=0 * pq.ms, ids=[]):
    """
    Given a list *sts* of spike trains, returns the times where *n* or more
    spike trains have exactly synchronous spikes, and the neuron ids
    firing at those times

    *Args*
    ------
    sts [list]:
        a list of neo spike trains

    n [int]:
        minimum number of coincident spikes to report synchrony

    dt [Quantity. Default to 0 ms]:
        size of time lag for synchrony. Starting from the very first spike,
        spike trains are binned at a binsize *dt* (bins half-open to the
        right), and spikes in the same bin are considered synchronous. If 0
        (default)

    ids [list. Default to []]:
        List of neuron IDs. Id[i] is the id associated to spike train
        sts[i]. If empty list provided (default), ids are assigned as
        integers from 0 to n_spiketrains-1.


    *Returns*
    ---------
    tracts [list of ndarrays]:
        a list of transactions, each as an array of neuron ids

    times [Quantity]
        a list of coincidence times. coinc[i] is the time of tracts[i]
    """
    # Convert list of spike trains to time-sorted gdf
    gdf, time_unit = spiketrains2gdf(neo_spiketrains, ids=ids)
    dt_dimless = dt.rescale(time_unit).magnitude  # make dt dimension-free
    # if dt_dimless is 0, set to half the minimum non-zero ISI
    if dt_dimless == 0:
        dt_dimless = np.diff(np.unique(gdf[:, 1])).min() / 2.

    # TODO: Clean up comments
    # tracts, times = [], []  # Initialize transactions and times to be
    # returned to empy lists
    idx_synch = []  # numpy.empty((0,2))
    # Set the initial time for the synchrony search to first available spike
    # time
    time = gdf[0, 1]
    # idx_start, idx_stop = 0, 0 # Starting from the first row in the gdf
    idx_start, idx_stop = 0, 0  # starting from the very first spike in the gdf
    while idx_stop <= gdf.shape[0] - 2:  # until end of gdf is reached,
        # Until the next spike falls in [time, time+dt] (symmetrical!)
        while time <= gdf[idx_stop + 1, 1] <= time + dt_dimless:
            idx_stop += 1  # Include that spike in the transaction
            if idx_stop >= gdf.shape[0] - 1:
                break  # And stop if end of gdf reached
        # If at least n spikes fall between idx_start and idx_stop
        if idx_stop >= idx_start + n - 1:
            idx_synch.extend(range(idx_start, idx_stop + 1))
        idx_start += 1  # Set new idx_start to the next spike
        idx_stop = idx_start  # and idx_stop to idx_start
        time = gdf[idx_stop, 1]  # and set the new corresponding spike time.
    # check if last entries also fulfill synchrony condition
    # If at least n spikes fall between idx_start and idx_stop
    if idx_stop >= idx_start + n - 1:
        idx_synch.extend(range(idx_start, idx_stop + 1))

    idx_synch = np.array(np.unique(idx_synch), dtype=int)

    # Return transactions of >=n synchronous spikes, and the times of these
    # transactions (first spike time in each transaction)
    # #return tracts, times*time_unit, dt_dimless*time_unit
    return gdf[idx_synch][:, 1] * time_unit, dt_dimless * time_unit






