import os
import plotly.graph_objects as go

from plotly.subplots import make_subplots

# from ..filesystem import JitterFileData
# from ..monitor_helpers import absolute_time

from cosmo.filesystem import JitterFileData
from cosmo.monitor_helpers import absolute_time

# TODO: Jitter Monitor needs to include 2 (possibly 3) subplots of statistics on the jitter per exposure.
#  Plot 1: mean vs time with +/- std "line boundaries"
#  Plot 2: max vs time
#  Plot 3 (potentially): std vs time


def view_jitter(jitter_file: str, figure) -> go.Figure:
    """Plot Jitter data for a single jitter root or association."""
    # Get the data
    jitter_data = JitterFileData(
        jitter_file,
        primary_header_keys=('PROPOSID', 'CONFIG'),
        ext_header_keys=('EXPNAME',),
        table_keys=('SI_V2_AVG', 'SI_V3_AVG', 'SI_V2_RMS', 'SI_V3_RMS', 'Seconds'),
        get_expstart=True
    )

    for jit in jitter_data:
        for direction, position in zip(['V2', 'V3'], [(1, 1), (2, 1)]):
            # Set time
            time = absolute_time(expstart=jit['EXPSTART'], time=jit['Seconds'])

            # Lower bound
            figure.add_trace(
                go.Scatter(
                    x=time.to_datetime(),
                    y=jit[f'SI_{direction}_AVG'] - jit[f'SI_{direction}_RMS'],
                    mode='lines',
                    line={'width': 0},
                    showlegend=False,
                    legendgroup=jit['EXPNAME'],
                    name=f'{jit["EXPNAME"]} -RMS'
                ),
                row=position[0],
                col=position[-1]
            )

            # Upper bound
            figure.add_trace(
                go.Scatter(
                    x=time.to_datetime(),
                    y=jit[f'SI_{direction}_AVG'] + jit[f'SI_{direction}_RMS'],
                    mode='lines',
                    line={'width': 0},
                    fillcolor='rgba(68, 68, 68, 0.1)',
                    fill='tonexty',
                    showlegend=False,
                    legendgroup=jit['EXPNAME'],
                    name=f'{jit["EXPNAME"]} +RMS'
                ),
                row=position[0],
                col=position[-1]
            )

            # Jitter
            figure.add_trace(
                go.Scatter(
                    x=time.to_datetime(),
                    y=jit[f'SI_{direction}_AVG'],
                    mode='lines',
                    legendgroup=jit['EXPNAME'],
                    name=f'{jit["EXPNAME"]} - {direction}',
                    hovertemplate='Time=%{x}<br>Jitter=%{y} arcseconds'
                ),
                row=position[0],
                col=position[-1]
            )

    figure.update_layout(
        {
            'title': f'{fits.getval(jitter_file, "PROPOSID")}/{fits.getval(jitter_file, "ROOTNAME")[4:6]}',
            'xaxis': {'title': 'Datetime'},
            'yaxis': {'title': 'V2 Jitter (averaged over 3 seconds) [arceconds]'},
            'xaxis2': {'title': 'Datetime'},
            'yaxis2': {'title': 'V3 Jitter (averaged over 3 seconds) [arcseconds]'}
        }
    )

    #figure.show()
    return figure

#-------------------------------------------------------------------------------
def add_annotations(figure, data_dir, row_pos):

    all_files = np.sort(glob.glob(os.path.join(data_dir, '*rawacq*')) +
                        glob.glob(os.path.join(data_dir, '*x1d.fits*')))
    for i, f in enumerate(all_files):
        if i % 2 == 1:
            txt_pos = -0.006
        else:
            txt_pos = 0.006

        expstart = fits.getval(f, 'expstart', ext=1)
        expend = fits.getval(f, 'expend', ext=1)
        figure.add_shape(type="line",
                      x0=Time(expstart, format='mjd').datetime,
                      y0=0.05,
                      x1=Time(expstart, format='mjd').datetime,
                      y1=-0.05,
                      line=dict(color="black",
                                width=3),
                    row=row_pos,
                    col=1)

        figure.add_shape(type="line",
                      x0=Time(expend, format='mjd').datetime,
                      y0=0.05,
                      x1=Time(expend, format='mjd').datetime,
                      y1=-0.05,
                      line=dict(color="grey",
                                width=3),
                    row=row_pos,
                    col=1)

        if 'ACQ' in fits.getval(f, 'exptype'):
            # for acqs
            start_label = f"{fits.getval(f, 'exptype')} start"
        else:
            # for science data
            start_label = f"{fits.getval(f, 'cenwave')} start"

        if row == 1:
            figure.add_annotation(
                        x=Time(expstart, format='mjd').datetime,
                        y=txt_pos,
                        text=start_label,
                        font=dict(family="Courier New, monospace",
                                  size=12,
                                  color="#ffffff"),
                        align="center",
                        arrowhead=2,
                        arrowsize=1,
                        arrowwidth=2,
                        arrowcolor="#636363",
                        bordercolor="#c7c7c7",
                        borderwidth=2,
                        borderpad=4,
                        bgcolor="grey",
                        opacity=0.9)

            figure.update_annotations(dict(
                        xref="x",
                        yref="y",
                        showarrow=True,
                        arrowhead=7,
            ))
    return figure

if __name__ == "__main__":

    import glob
    from astropy.io import fits
    from astropy.time import Time
    import os
    import numpy as np

    # Start the figure
    figure = make_subplots(2, 1, shared_xaxes=True, shared_yaxes=True)

    data_dir = "/user/rplesha/front_end/hybrid_mode_onorbit/" #day_170/15955_16"
    #data_dir = "/grp/hst/cos2/cosmo/15773/"

    for f in np.sort(glob.glob(os.path.join(data_dir, '*/*/*jit*.fits'))):
        pid = fits.getval(f, 'PROPOSID')
        figure = view_jitter(f, figure)

    for row in [1, 2]:
        figure = add_annotations(figure, data_dir, row)

    figure.update_yaxes(range=[-0.01, 0.01])

    figure.show()
    #figure.write_html(os.path.join(data_dir, f'{pid}_jitter_vs_time.html'))
    figure.write_html(os.path.join(data_dir, f'all_jitter_vs_time.html'))
