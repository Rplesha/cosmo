import plotly.graph_objs as go
import datetime

from itertools import repeat
from astropy.time import Time, TimeDelta

from monitorframe import BaseMonitor
from .osm_data_models import OSMDataModel

COS_MONITORING = '/grp/hst/cos2/monitoring'


class FUVADOSMShiftMonitor(BaseMonitor):
    data_model = OSMDataModel
    output = COS_MONITORING
    labels = ['ROOTNAME', 'LIFE_ADJ', 'FPPOS', 'PROPOSID', 'OBSET_ID']
    subplots = True
    subplot_layout = (2, 1)

    def track(self):
        # TODO: Define interesting things to track
        pass

    def filter_data(self):
        return self.data[self.data.DETECTOR == 'FUV']

    @staticmethod
    def compute_start_times(df):
        start_time = Time(df.EXPSTART, format='mjd')
        lamp_dt = TimeDelta(df.TIME, format='sec')
        lamp_time = start_time + lamp_dt

        return start_time, lamp_time

    def plot(self):  # TODO: Add A - B subplot
        groups = self.filtered_data.groupby(['OPT_ELEM', 'CENWAVE'])
        lp4_move = datetime.datetime.strptime('2017-10-02', '%Y-%m-%d')

        for i, group_info in enumerate(groups):
            name, group = group_info

            start_time, lamp_time = self.compute_start_times(group)

            self.figure.append_trace(
                go.Scattergl(
                    x=lamp_time.to_datetime(),
                    y=group.SHIFT_DISP,
                    name='-'.join([str(item) for item in name]),
                    mode='markers',
                    marker=dict(
                        cmax=18,
                        cmin=0,
                        color=list(repeat(i, len(group))),
                        colorscale='Viridis',
                        symbol=group.FPPOS * 4,
                        size=[
                            10 if time > lp4_move and lp == 3 else 6
                            for lp, time in zip(group.LIFE_ADJ, start_time.to_datetime())
                        ]
                    ),
                    text=group.hover_text
                ), 1, 1
            )

        segment_diff = []
        time = []
        hover_text = []
        root_groups = self.filtered_data.groupby('ROOTNAME')
        for rootname, group in root_groups:
            if 'FUVA' in group.SEGMENT.values and 'FUVB' in group.SEGMENT.values:
                _, lamp_time = self.compute_start_times(group[group.SEGMENT == 'FUVA'])
                segment_diff.extend(
                    group[group.SEGMENT == 'FUVA'].SHIFT_DISP.values - group[group.SEGMENT == 'FUVB'].SHIFT_DISP.values
                )
                hover_text.extend(group[group.SEGMENT == 'FUVA'].hover_text.values)
                time.extend(lamp_time.to_datetime())

        self.figure.append_trace(
            go.Scattergl(
                x=time,
                y=segment_diff,
                name='FUVA - FUVB',
                mode='lines+markers',
                text=hover_text,
                xaxis='x'
            ), 2, 1
        )

        layout = go.Layout(
                xaxis=dict(title='Datetime'),
                yaxis=dict(title='AD Shift [pix]', anchor='x', domain=[0.3, 1]),
                yaxis1=dict(title='Shift Difference A - B [pix]', anchor='x2', domain=[0, 0.18])
            )

        self.figure['layout'].update(layout)

    def store_results(self):
        # TODO: decide on how to store results and what needs to be stored
        pass
