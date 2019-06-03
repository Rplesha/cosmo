import plotly.graph_objs as go
import datetime
import pandas as pd
import os

from itertools import repeat

from monitorframe import BaseMonitor
from cosmo.monitor_helpers import compute_lamp_on_times
from .osm_data_models import OSMDataModel

COS_MONITORING = '/grp/hst/cos2/monitoring'

LP_MOVES = lp_moves = {
        i + 2: datetime.datetime.strptime(date, '%Y-%m-%d')
        for i, date in enumerate(['2012-07-23', '2015-02-09', '2017-10-02'])
    }


def plot_fuv_osm_shift_cenwaves(df, shift):
    groups = df.groupby(['OPT_ELEM', 'CENWAVE'])

    fp_symbols = {
        1: 'circle',
        2: 'cross',
        3: 'triangle-up',
        4: 'x'
    }

    traces = []

    seg_diff_results = compute_segment_diff(df, shift)

    traces.append(
        go.Scattergl(
            x=seg_diff_results.lamp_time,
            y=seg_diff_results.seg_diff,
            name='FUVA - FUVB',
            mode='markers',
            text=seg_diff_results.hover_text,
            visible=False
        )
    )

    for i, group_info in enumerate(groups):
        name, group = group_info

        start_time, lamp_time = compute_lamp_on_times(group)

        traces.append(
            go.Scattergl(
                x=lamp_time.to_datetime(),
                y=group[shift],
                name='-'.join([str(item) for item in name]),
                mode='markers',
                text=group.hover_text,
                xaxis='x',
                yaxis='y2',
                visible=False,
                marker=dict(
                    cmax=len(df.CENWAVE.unique()) - 1,
                    cmin=0,
                    color=list(repeat(i, len(group))),
                    colorscale='Viridis',
                    symbol=[fp_symbols[fp] for fp in group.FPPOS],
                    size=[
                        10 if time > LP_MOVES[4] and lp == 3 else 6
                        for lp, time in zip(group.LIFE_ADJ, start_time.to_datetime())
                    ]
                ),
            )
        )

    layout = go.Layout(
        xaxis=dict(title='Datetime'),
        yaxis2=dict(title='Shift [pix]', anchor='x', domain=[0.3, 1], gridwidth=5),
        yaxis=dict(title='Shift (A - B) [pix]', anchor='x2', domain=[0, 0.18])
    )

    return traces, layout


def compute_segment_diff(df, shift):
    root_groups = df.groupby('ROOTNAME')

    results_list = []
    for rootname, group in root_groups:
        if 'FUVA' in group.SEGMENT.values and 'FUVB' in group.SEGMENT.values:
            _, lamp_time = compute_lamp_on_times(group[group.SEGMENT == 'FUVA'])

            fuva, fuvb = group[group.SEGMENT == 'FUVA'], group[group.SEGMENT == 'FUVB']

            diff_df = fuva.assign(seg_diff=fuva[shift].values - fuvb[shift].values).reset_index(drop=True)

            diff_df['hover_text'] = [
                item.hover_text + f'<br>CENWAVE    {item.CENWAVE}' for _, item in fuva.iterrows()
            ]

            diff_df['lamp_time'] = lamp_time.to_datetime()
            diff_df.drop(columns=['TIME', 'SHIFT_DISP', 'SHIFT_XDISP', 'SEGMENT', 'DETECTOR'], inplace=True)
            results_list.append(diff_df)

    results_df = pd.concat(results_list, ignore_index=True)

    return results_df


class FuvOsmShiftMonitor(BaseMonitor):
    data_model = OSMDataModel
    output = COS_MONITORING
    labels = ['ROOTNAME', 'LIFE_ADJ', 'FPPOS', 'PROPOSID']

    shift = None

    def track(self):
        return compute_segment_diff(self.filtered_data, self.shift)

    def filter_data(self):
        return self.data[self.data.DETECTOR == 'FUV']

    def plot(self):
        traces, layout = plot_fuv_osm_shift_cenwaves(self.filtered_data, self.shift)

        trace_length = len(traces)
        all_visible = list(repeat(True, trace_length))

        fp_groups = self.filtered_data.groupby('FPPOS')
        fp_trace_lengths = {}
        for fp, group in fp_groups:
            fp_traces, _ = plot_fuv_osm_shift_cenwaves(group, self.shift)
            traces.extend(fp_traces)
            all_visible.extend(list(repeat(False, len(fp_traces))))
            fp_trace_lengths[fp] = len(fp_traces)

        self.figure.add_traces(traces)
        self.figure['layout'].update(layout)

        all_hidden = list(repeat(False, len(all_visible)))

        # TODO: Find a better way to track which traces to turn "on" and "off". This is gross.
        fp1_visible = all_hidden.copy()
        fp2_visible = all_hidden.copy()
        fp3_visible = all_hidden.copy()
        fp4_visible = all_hidden.copy()

        fp1_visible[trace_length:trace_length + fp_trace_lengths[1]] = list(repeat(True, fp_trace_lengths[1]))

        fp2_visible[
            trace_length + fp_trace_lengths[1]:trace_length + fp_trace_lengths[1] + fp_trace_lengths[2]
        ] = list(repeat(True, fp_trace_lengths[2]))

        fp3_visible[
            trace_length + fp_trace_lengths[1] + fp_trace_lengths[2]:
            trace_length + fp_trace_lengths[1] + fp_trace_lengths[2] + fp_trace_lengths[3]
        ] = list(repeat(True, fp_trace_lengths[3]))

        fp4_visible[
            trace_length + fp_trace_lengths[1] + fp_trace_lengths[2] + fp_trace_lengths[3]:
        ] = list(repeat(True, fp_trace_lengths[4]))

        button_labels = ['All FPPOS', 'FPPOS 1', 'FPPOS 2', 'FPPOS 3', 'FPPOS 4']
        vsibilities = [all_visible, fp1_visible, fp2_visible, fp3_visible, fp4_visible]
        titles = [f'{self.name} All FPPOS'] + [f'{self.name} {label} Only' for label in button_labels[1:]]

        updatemenues = [
            dict(
                type='buttons',
                buttons=[
                    dict(
                        label=label,
                        method='update',
                        args=[
                            {'visible': visibility},
                            {'title': button_title}
                        ]
                    ) for label, visibility, button_title in zip(button_labels, vsibilities, titles)
                ]
            )
        ]

        lines = [
            {
                'type': 'line',
                'x0': lp_time,
                'y0': 0,
                'x1': lp_time,
                'y1': 1,
                'yref': 'paper',
                'line': {
                    'width': 2,
                }
            } for lp_time in LP_MOVES.values()
        ]

        self.figure['layout'].update({'shapes': lines, 'updatemenus': updatemenues})

    def store_results(self):
        outliers = self.results[self.outliers]
        outliers.to_csv(os.path.join(os.path.dirname(self.output), f'{self._filename}-outliers.csv'))


class FuvOsmShift1Monitor(FuvOsmShiftMonitor):
    shift = 'SHIFT_DISP'

    def find_outliers(self):
        return self.results.seg_diff.abs() > 10


class FuvOsmShift2Monitor(FuvOsmShiftMonitor):
    shift = 'SHIFT_XDISP'

    def find_outliers(self):
        return self.results.seg_diff.abs() > 5


class NuvOsmShiftMonitor(BaseMonitor):
    data_model = OSMDataModel
    output = COS_MONITORING
    labels = ['ROOTNAME', 'LIFE_ADJ', 'FPPOS', 'PROPOSID']

    shift = None

    def filter_data(self):
        return self.data[self.data.DETECTOR == 'NUV']

    def plot(self):
        groups = self.filtered_data.groupby(['OPT_ELEM', 'CENWAVE'])

        traces = []
        for i, group_info in enumerate(groups):
            name, group = group_info

            start_time, lamp_time = compute_lamp_on_times(group)

            traces.append(
                go.Scattergl(
                    x=lamp_time.to_datetime(),
                    y=group[self.shift],
                    name='-'.join([str(item) for item in name]),
                    mode='markers',
                    text=group.hover_text,
                    marker=dict(
                        cmax=len(self.filtered_data.CENWAVE.unique()) - 1,
                        cmin=0,
                        color=list(repeat(i, len(group))),
                        colorscale='Viridis',
                    ),
                )
            )

        layout = go.Layout(
            xaxis=dict(title='Datetime'),
            yaxis=dict(title='Shift [pix]', gridwidth=5),
        )

        self.figure.add_traces(traces)
        self.figure['layout'].update(layout)

    def store_results(self):
        # TODO: decide on what results to store and how
        pass

    def track(self):
        # TODO: decide on what to track
        pass


class NuvOsmShift1Monitor(NuvOsmShiftMonitor):
    shift = 'SHIFT_DISP'


class NuvOsmShift2Monitor(NuvOsmShiftMonitor):
    shift = 'SHIFT_XDISP'
