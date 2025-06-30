import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# import dash_bio

from typing import List

import json
from math import log10
import pandas as pd
import numpy as np
from itertools import chain

from backend.types import MetaPepDbSearch, MetaPepDeNovo
from backend.type_operations import assert_db_search_compatibility, assert_de_novo_compatibility
from backend.post_processing import reference_score_distribution, reference_score_dist_peaks
from backend.spectral_ref_builder import *
from backend.utils.graph_utils import *
from backend.utils.pd_utils import *
from backend.io.import_spectra import *
from constants import *




def tic_over_rt_plot(spectral_dataset: pd.DataFrame,
                     mzml_peaks_dataset: str | None,
                     peaks_compression: str,
                     peaks_precision: str,
                     ms_level: int, 
                     secondary_param: str | None, 
                     sma_window: int, 
                     data_reduction_factor: int=1,
                     metapep_dataset: MetaPepDbSearch | MetaPepDeNovo | None = None,
                     confidence_threshold: int | float | None = None,
                     peak_int_threshold: int | float | None = None):
    # Compute moving average of intensity over rt
    spectral_dataset = spectral_dataset[spectral_dataset['MS level'] == ms_level]
    spectral_dataset = spectral_dataset[['scan number', 'total ion current', 'retention time', 'peaks count']]
    spectral_dataset['total ion current SMA'] = spectral_dataset['total ion current']\
        .rolling(sma_window, closed='both')\
        .mean()
    
    # cut out data from dataset in equally spaced parts for faster rendering
    spectral_dataset = spectral_dataset.iloc[::data_reduction_factor]

    # scattergl implementation
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scattergl(x=spectral_dataset['retention time'],
                               y=spectral_dataset['total ion current SMA'],
                               mode='lines',
                               name="total ion current",
                               line=dict(width=1, color=GraphConstants.primary_color)))
    fig.update_yaxes(title="Total ion current",
                     secondary_y=False)
    
    # add secondary parameter if defined
    if secondary_param == 'Confidence' and metapep_dataset is not None:
        metapep_data = metapep_dataset.data
        valid_rts = metapep_data[metapep_data[secondary_param] > confidence_threshold]['RT']

        name = f"{metapep_dataset.confidence_format} > {confidence_threshold}"
        
        if valid_rts.shape[0] != 0:
            fig.add_trace(go.Histogram(x=valid_rts,
                                       nbinsx=40,
                                       name=name,
                                       opacity=0.4),
                          secondary_y=True)
            fig.update_yaxes(title="Peptide matches",
                             secondary_y=True,
                             tickmode="sync")
    elif secondary_param == "Peak Count":
        # reprocess peak count by filtering intensities below threshold from peak array
        if peak_int_threshold is not None and \
            peak_int_threshold > 0 and\
            mzml_peaks_dataset is not None:
                
            mzml_peaks_dict: Dict = json.loads(mzml_peaks_dataset)

            def get_peaks_count_threshold(mzml_row: pd.Series) -> int:
                peak_arr = fetch_mzml_peaks_data(
                    mzml_row=mzml_row,
                    peaks_dict=mzml_peaks_dict,
                    peaks_compression=peaks_compression,
                    peaks_precision=peaks_precision,
                    include_mz=False
                )
                return len(peak_arr[peak_arr > peak_int_threshold]) #type:ignore
            
            # count peaks above specified threshold by processing peak intensities
            spectral_dataset['peaks count'] = spectral_dataset[["scan number", "peaks count"]].apply(
                get_peaks_count_threshold,
                axis=1,
                result_type="expand")
            trace_name = f'peak count > {peak_int_threshold}'
        else:
            trace_name = f'peak count'
            
        spectral_dataset['peaks count SMA'] = spectral_dataset['peaks count']\
            .rolling(sma_window, closed='both')\
            .mean()
            
        
        fig.add_trace(go.Scatter(x=spectral_dataset['retention time'], 
                                 y=spectral_dataset['peaks count SMA'], 
                                 mode='lines', 
                                 name=trace_name,
                                 line=dict(width=1, color="red")),
                      secondary_y=True)
        fig.update_yaxes(title="Peaks count",
                         secondary_y=True,
                         tickmode="sync")
    
    if secondary_param == "None":
        fig.update_layout(showlegend=False)
    
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth,
                     zerolinecolor="Black",
                     range=[0, None], 
                     nticks=4)
    fig.update_xaxes(gridcolor="rgba(0,0,0,0)", 
                     zerolinecolor="Black",
                     title="Retention time")
    return fig


def ms2_from_signal_arrays(mz_array, int_array):
    dataset = go.Bar(x=mz_array, y=int_array, width=2, marker={"color": "black"})
    
    fig = go.Figure()
    
    fig.update_xaxes(title="m/z", range=[0, max(mz_array)+100])
    fig.update_yaxes(title="Ion Intensity", range=[0, max(int_array)*1.1])

    # add signals as line in ms2 figure
    # for mz_val, int_val in :
    #     fig.add_shape(type="line",
    #                   x0=mz_val, x1=mz_val,
    #                   y0=0, y1=int_val,
    #                   line=dict(color="Black",
    #                             width=1,
    #                             )
    #                   )
    fig.update_layout(shapes=[{'type': 'line','y0':0,'y1': int_val,
                               'x0': mz_val, 'x1': mz_val,
                               'line': {'color': GraphConstants.primary_color,
                                        'width': 1}} 
                              for mz_val, int_val in zip(mz_array, int_array)])
    
    fig.update_shapes(dict(xref='x', yref='y'))

    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10))
    
    return fig


def mz_over_rt_plot(dataset: pd.DataFrame,
                    db_search_psm: MetaPepDbSearch | None = None,
                    de_novo: MetaPepDeNovo | None = None,
                    int_cutoff: int | float | None = None):
    # Compute moving average of intensity over rt
    tic_df = dataset[dataset['MS level'] == 1]
    tic_df = tic_df[['scan number', 'total ion current', 'retention time']]
    tic_df['total ion current EMA'] = tic_df['total ion current'].ewm(span=100).mean()

    dataset = dataset[['scan number', 'precursor m/z', 'retention time', 'total ion current']]
    dataset = dataset.astype(float)
    
    # add tic cutoff constraint to data if specified
    if int_cutoff is not None:
        dataset = dataset[dataset['total ion current'] > int_cutoff]
    
    fig = make_subplots(rows=1, cols=2, shared_yaxes=True,
                        column_widths=[0.15, 0.85], horizontal_spacing=0.02)

    fig.add_trace(go.Scatter(x=tic_df["total ion current EMA"],
                             y=tic_df["retention time"],
                             mode="lines",
                             showlegend=False,
                             line=dict(color=GraphConstants.primary_color,
                                       width=1)),
                  row=1, 
                  col=1)
    fig.update_xaxes(title="TIC", row=1, col=1, autorange="reversed",
                     showgrid=False, zerolinecolor='Black')
    fig.update_yaxes(title="retention time", zerolinecolor="Black", showgrid=False, row=1, col=1)
    #px.line(tic_df, x="retention time", y="total ion current EMA")
    #fig.update_layout(margin=dict(l=20, r=20, t=10, b=10))
    
    # specify for each cel:l if it matches to db_search_psm
    if db_search_psm is not None:
        dataset = match_db_search_psm(dataset, db_search_psm.data)
        # fig = px.scatter(dataset, x="precursor m/z", y="retention time", color="db search identified")
        
        palette = iter(GraphConstants.color_palette)
        for i, category in enumerate(dataset['db search identified'].unique()):
            cat_df = dataset[dataset['db search identified'] == category]
            # fig.add_trace(go.Histogram2dContour(
            #     x=cat_df["precursor m/z"], y=cat_df["retention time"], opacity=0.4, colorscale=scales[i]),
            #               row=1, col=2)
            fig.add_trace(go.Scattergl(x=cat_df["precursor m/z"], y=cat_df["retention time"],
                                       mode='markers',
                                       name=str(category),
                                       marker_size=2,
                                       marker_color=next(palette)),
                          row=1, col=2)
        
        fig.update_layout(legend={'title': 'db search identified', 'itemsizing': 'constant'})
        
    elif de_novo is not None:
        dataset = match_de_novo(dataset, de_novo.data)
        # fig = px.scatter(dataset, x="precursor m/z", y="retention time", color="de novo identified")
        
        palette = iter(GraphConstants.color_palette)
        for category in dataset['de novo identified'].unique():
            cat_df = dataset[dataset['de novo identified'] == category]
            fig.add_trace(go.Scattergl(x=cat_df["precursor m/z"], y=cat_df["retention time"],
                                       mode='markers',
                                       name=str(category),
                                       marker_size=2,
                                       marker_color=next(palette)),
                          row=1, col=2)

        fig.update_layout(legend= {'title': 'de novo identified','itemsizing': 'constant'})
        
    else:
        # fig = px.scatter(dataset, x="precursor m/z", y="retention time")
        # fig.add_trace(go.Histogram2dContour(
        #     x=dataset["precursor m/z"], y=dataset["retention time"], colorscale='Blues', showlegend=False),
        #               row=1, col=2)
        
        
        fig.add_trace(go.Scattergl(x=dataset["precursor m/z"], y=dataset["retention time"], mode='markers',
                                 showlegend=False),
                      row=1, col=2)
        fig.update_traces(marker=dict(size=2,
                                      color=GraphConstants.primary_color),
                          row=1,
                          col=2)
        
    fig.update_xaxes(title="precursor m/z",
                     gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth, 
                     zerolinecolor="Black",
                     row=1,
                     col=2)
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     zerolinecolor="Black",
                     row=1,
                     col=2)
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    return fig


def scan_tic_dist_plot(dataset: pd.DataFrame,
                       ms_level: int,
                       db_search_psm: MetaPepDbSearch | None = None,
                       de_novo: MetaPepDeNovo | None = None,
                       alc_cutoff: int = 0,
                       normalize_bars: bool = False):
    # get correct columns and match to db search and de novo
    dataset = dataset[['scan number', 'total ion current', 'MS level']]
    dataset = dataset[dataset['MS level'] == ms_level]

    barnorm = None if normalize_bars is False else 'fraction'
    
    if db_search_psm is not None and ms_level != 1:
        dataset = match_db_search_psm(dataset, db_search_psm.data)
    else:
        dataset['db search identified'] = np.nan
    
    if de_novo is not None and ms_level != 1:
        # fetch confidence column name for specific de novo source for axis naming
        de_novo_conf_format = de_novo.confidence_format
        de_novo_data = de_novo.data
        # filter peptides under confidence cutoff
        if alc_cutoff > 0:
            de_novo_data = de_novo_data[de_novo_data["Confidence"] > alc_cutoff]
        dataset = match_de_novo(dataset, de_novo_data)
    else:
        de_novo_conf_format = None
        dataset['de novo identified'] = np.nan
    
    # categorize scans by db search, de novo or no identification
    def partition_ident(row):
        if row["db search identified"] is True:
            return "db search"
        elif row["de novo identified"] is True:
            return "de novo only (> {} {})".format(alc_cutoff, de_novo_conf_format)
        else:
            return "unidentified"
    dataset['identification'] = dataset.apply(partition_ident, axis=1)
    
    # log transform tic
    dataset['total ion current'] = dataset['total ion current'].apply(log10)
    
    if (db_search_psm is not None or de_novo is not None) and ms_level != 1:
        fig = px.histogram(dataset, 
                           x="total ion current", 
                           color="identification", 
                           nbins=30,
                           barnorm=barnorm,
                           color_discrete_sequence=GraphConstants.color_palette,
                           category_orders={"identification": ["db search", 
                                                               "de novo only (> {} {})".format(alc_cutoff, de_novo_conf_format),
                                                               "unidentified"]})
    else:
        fig = px.histogram(dataset, 
                           x="total ion current", 
                           color_discrete_sequence=[GraphConstants.primary_color], 
                           nbins=30)
    
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    fig.update_xaxes(title="log10(tic)", 
                     showline=True, 
                     linecolor="Black")
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth)
    return fig


def ms1_int_over_ms2_int(dataset: pd.DataFrame,
                         peaks: str | None = None,
                         db_search_psm: MetaPepDbSearch | None = None,
                         de_novo: MetaPepDeNovo | None = None,
                         min_mz: int = 0,
                         peaks_compression: str = 'no compression',
                         peaks_precision: str = "64-bit float"):
    if min_mz <= 0:
        field_list = ['scan number', 'total ion current', 'precursor intensity', 'retention time']
    else:
        field_list = ['scan number', 'total ion current', 'precursor intensity', 'retention time',
                      'peaks count']
    
    dataset = dataset[dataset['MS level'] == 2]
    dataset = dataset[field_list]

    # obtain peak arrays
    if min_mz > 0 and peaks is not None:
        peaks_dict: Dict = json.loads(peaks)

        def fetch_peaks_data(df_row: pd.Series):
            spectrum_peaks = peaks_dict[str(df_row.name)]
            mz_arr = decode_mzml_peaks(spectrum_peaks["m/z array"],
                                       df_row['peaks count'],
                                       compression_type=peaks_compression,
                                       precision=peaks_precision)
            int_arr = decode_mzml_peaks(spectrum_peaks["intensity array"],
                                        df_row['peaks count'],
                                        compression_type=peaks_compression,
                                        precision=peaks_precision)
            return (mz_arr, int_arr)

        dataset[['m/z array', 'intensity array']] = dataset[["scan number", "peaks count"]].apply(
            fetch_peaks_data,
            axis=1, result_type="expand")
               
    # recompute totall ion current for peaks above mz cutoff in new column
    if min_mz > 0:
        y_col = 'total ion current > {}'.format(min_mz)
        # filter from intensity array values with m/z under cutoff, then compute sum
        dataset[y_col] = dataset.apply(
            lambda x: np.sum(x['intensity array'][np.where(x['m/z array'] > min_mz)]),
            axis=1)
    else:
        y_col = 'total ion current'
    
    fig = go.Figure()
    
    # match to peptide datases
    if db_search_psm is not None:
        dataset = match_db_search_psm(dataset, db_search_psm.data)
        
        palette = iter(GraphConstants.color_palette)
        for category in dataset['db search identified'].unique():
            cat_df = dataset[dataset['db search identified'] == category]
            fig.add_trace(go.Scattergl(x=cat_df["precursor intensity"], y=cat_df[y_col],
                                       name=str(category), mode='markers', opacity=0.5,
                                       marker_size=2,
                                       marker_color=next(palette)))
            fig.update_layout(legend= {'title': 'db search identified','itemsizing': 'constant'})
    elif de_novo is not None:
        dataset = match_de_novo(dataset, de_novo.data)
        
        palette = iter(GraphConstants.color_palette)
        for category in dataset['de novo identified'].unique():
            cat_df = dataset[dataset['de novo identified'] == category]
            fig.add_trace(go.Scattergl(x=cat_df["precursor intensity"], 
                                       y=cat_df[y_col],
                                       name=str(category), 
                                       mode='markers', 
                                       opacity=0.5,
                                       marker_size=2,
                                       marker_color=next(palette)))
            fig.update_layout(legend= {'title': 'de novo identified','itemsizing': 'constant'})
    else:
        fig.add_trace(go.Scattergl(x=dataset["precursor intensity"], 
                                   y=dataset[y_col], 
                                   mode='markers', 
                                   opacity=0.5,
                                   marker_size=2,
                                   marker_color=GraphConstants.primary_color))
    
    if min_mz > 0:
        ytitle = f"MS2 tic (>{min_mz})"
    else:
        ytitle = "MS2 tic"
    
    fig.update_layout(legend={'itemsizing': 'constant'})
    fig.update_xaxes(title="precursor intensity", 
                     type='log', 
                     gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth, 
                     zerolinecolor="Black", 
                     nticks=6)
    fig.update_yaxes(title=ytitle, 
                     type='log', 
                     gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth,
                     zerolinecolor="Black",
                     nticks=6)
    
    # if min_mz > 0 and peaks is not None:
    #     fig.update_layout(title="relation precursor tic to MS2 tic (>{} m/z)".format(min_mz))
    # else:
    #     fig.update_layout(title="relation precursor tic to MS2 tic")
    
    upper_bound = dataset[["total ion current", "precursor intensity"]].max(axis=None)
    fig.add_shape(type="line", x0=1000, y0=1000, x1=upper_bound, y1=upper_bound,
                  line=dict(width=2, dash="dot"))
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)')
    return fig
    

def confidence_dist_plot(db_search: MetaPepDbSearch | None,
                         de_novo: MetaPepDeNovo | None,
                         db_search_conf_cutoff: float,
                         de_novo_conf_cutoff: float,
                         plot_all: bool=False):
    db_search_ranking, de_novo_ranking = None, None
    if db_search is not None:
        db_search_ranking = db_search.data["Confidence"].sort_values(ascending=False)\
            .reset_index(drop=True)
        if db_search_conf_cutoff > 0:
            db_search_ranking = db_search_ranking[db_search_ranking > db_search_conf_cutoff]
    if de_novo is not None:
        de_novo_ranking = de_novo.data["Confidence"].sort_values(ascending=False)\
            .reset_index(drop=True)
        if de_novo_conf_cutoff > 0:
            de_novo_ranking = de_novo_ranking[de_novo_ranking > de_novo_conf_cutoff]


    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    # configure plot traces and axes based on presence of proteomics data
    if db_search is not None and de_novo is not None:
        # resize to both to length of shortest
        shape = min(db_search_ranking.shape[0], de_novo_ranking.shape[0]) #type:ignore
        db_search_ranking = db_search_ranking.head(shape) #type:ignore
        de_novo_ranking = de_novo_ranking.head(shape) #type:ignore
        
        confidence_series = db_search_ranking / de_novo_ranking    
        
        ytitle = "{}, {}".format(db_search.confidence_format, de_novo.confidence_format)
        y2title = "{} / {}".format(db_search.confidence_format, de_novo.confidence_format)
        
        if plot_all is True:
            fig.add_trace(
                go.Scattergl(x=confidence_series.index, y=confidence_series.values,
                            mode="lines", name="DB Search / De novo"),
                secondary_y=True
            )
            fig.update_yaxes(title_text=y2title, range=[0.0, confidence_series.max()], secondary_y=True)
            fig.update_yaxes(title_text=ytitle, range=[0.0, max(db_search_ranking.max(), 100)], secondary_y=False)
        else:
            fig.add_trace(
                go.Scattergl(x=confidence_series.index, y=confidence_series.values,
                            mode="lines", name="DB Search / De novo"))
            fig.update_yaxes(title_text=y2title, range=[0.0, confidence_series.max()])
    elif db_search is not None:
        fig.update_yaxes(title_text=db_search.confidence_format, range=[0.0, db_search_ranking.max()]) #type:ignore
    elif de_novo is not None:
        fig.update_yaxes(title_text=de_novo.confidence_format, range=[0.0, 100])
    else:
        return None
        
    if db_search_ranking is not None:
        if de_novo_ranking is None or plot_all is True:
            fig.add_trace(
                go.Scattergl(x=db_search_ranking.index, y=db_search_ranking.values, mode="lines", name="DB Search")
            )
    
    if de_novo_ranking is not None:
        if db_search_ranking is None or plot_all is True:
            fig.add_trace(
                go.Scattergl(x=de_novo_ranking.index, y=de_novo_ranking.values, mode="lines", name="De novo")
            )
        

    # fig.update_traces(line_color=primary_color)
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      yaxis2=dict(tickmode="sync")
                      )
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor,
                     gridwidth=GraphConstants.gridwidth, 
                     zerolinecolor="Black",
                     secondary_y=False, 
                     nticks=6)
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth, 
                     zerolinecolor="Black",
                     secondary_y=True, 
                     nticks=6)
    
    fig.update_xaxes(title="Peptide rank", gridcolor="rgba(0,0,0,0)", zerolinecolor="Black")
    
    return fig


def ref_score_dist_plot(stat_dict: dict,
                        sample_db_search: MetaPepDbSearch | None,
                        sample_de_novo: MetaPepDeNovo | None,
                        format: str,
                        normalize_scans: bool = False,
                        normalize_matches: bool = False,
                        sample_ms2_count: int | None = None):
    fig = go.Figure()
    
    if (format == "db search" and sample_db_search is None) or \
        sample_de_novo is None:
        plot_sample = False
    else:
        plot_sample = True
    
    if sample_ms2_count is None and plot_sample:
        normalize_scans = False
        print("No ms2 count in sample data, skip scan normalization...")

    if format == "db search":
        if normalize_scans is True:
            array_key = "db search confidence dist norm ms2"
        elif normalize_matches is True:
            array_key = "db search confidence dist norm"
        else:
            array_key = "db search confidence dist"
        score = stat_dict['metadata']['db search confidence format']
    elif format == "de novo":
        if normalize_scans is True:
            array_key = "de novo confidence dist norm ms2"
        elif normalize_matches is True:
            array_key = "de novo confidence dist norm"
        else:
            array_key = "de novo confidence dist"
        score = stat_dict['metadata']['de novo confidence format']
    else:
        raise ValueError(f"Invalid format supplied: {format}")
    
    # do not process db search and de novo if they are not compatible
    if sample_db_search is not None and \
        assert_db_search_compatibility(stat_dict, sample_db_search) is False:
        sample_db_search = None
    elif sample_de_novo is not None and\
        assert_de_novo_compatibility(stat_dict, sample_de_novo) is False:
        sample_de_novo = None
    
    # fetch arrays
    if normalize_scans is False and normalize_matches is False:
        mean_array = np.array(stat_dict['global'][array_key]["mean"])
        std_array =  np.array(stat_dict['global'][array_key]["std"])
        xvals = list(range(len(mean_array)))
    else:
        mean_array = np.array(stat_dict['global'][array_key]["mean"])
        std_array =  np.array(stat_dict['global'][array_key]["std"])
        xvals = np.array(stat_dict['global'][array_key]["x vals"])
    
    # return nothing if no data found in dict
    if len(mean_array) == 0:
        return

    # add confidence interval
    fig.add_trace(go.Scatter(x=xvals, y=mean_array, name="mean", marker=dict(color="#444")))
    fig.add_trace(go.Scatter(x=xvals, y=mean_array + std_array, name="upper bound", marker=dict(color="#444"), line=dict(width=0), showlegend=False))
    fig.add_trace(go.Scatter(x=xvals, y=np.clip(mean_array - std_array, 0, None), name="lower bound", marker=dict(color="#444"), line=dict(width=0),
                            fillcolor='rgba(68, 68, 68, 0.3)', fill='tonexty', showlegend=False))

    # overlay score distribution from sample
    if format == "db search" and sample_db_search is not None:
        if normalize_scans is True:
            y, x = pept_match_dist_normalize(
                sample_db_search.data,
                'Confidence',
                total_scans=sample_ms2_count
            )
        elif normalize_matches is True:
            y, x = pept_match_dist_normalize(
                sample_db_search.data,
                'Confidence'
            )
        else:
            sorted_scores = fetch_sort_column(sample_db_search.data, "Confidence")
            x = sorted_scores.index.to_numpy()
            y = sorted_scores.values

        fig.add_trace(go.Scatter(x=x,
                                 y=y,
                                 name="DB Search Import",
                                 marker=dict(color="red")))
    elif format == "de novo" and sample_de_novo is not None:
        if normalize_scans is True:
            y, x = pept_match_dist_normalize(
                sample_de_novo.data,
                'Confidence',
                total_scans=sample_ms2_count
            )
        elif normalize_matches is True:
            y, x = pept_match_dist_normalize(
                sample_de_novo.data,
                'Confidence'
            )
        else:
            sorted_scores = fetch_sort_column(sample_de_novo.data, "Confidence")
            x = sorted_scores.index.to_numpy()
            y = sorted_scores.values

        fig.add_trace(go.Scatter(x=x,
                                 y=y,
                                 name="De Novo Import",
                                 marker=dict(color="red")))

    if normalize_scans is True:
        xtitle = "score n'th peptide match (% of MS2)"
    elif normalize_matches is True:
        xtitle = "score n'th peptide match (% of matches)"
    else:
        xtitle = "n'th peptide match"
    
    fig.update_xaxes(title=xtitle, 
                     gridcolor="rgba(0,0,0,0)", 
                     zerolinecolor="Black")
    fig.update_yaxes(title=score, 
                     gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth, 
                     zerolinecolor="Black", 
                     nticks=6, 
                     rangemode='tozero')
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)'
                    )

    return fig


def ref_score_threshold_plot(stat_dict: dict,
                             formats: List[str],
                             normalize_psm: bool=False,
                             normalize_rt: bool=False,
                             sample_db_search: MetaPepDbSearch | None=None,
                             sample_de_novo: MetaPepDeNovo | None=None,
                             spectral_metadata: dict | None=None):
    db_search_score_unit = stat_dict['metadata']['db search confidence format']
    de_novo_score_unit = stat_dict['metadata']['de novo confidence format']
    
    # Define group label dict, each format has a prefix and suffix
    format_to_label = {"db search": [f"{db_search_score_unit}", ""],
                       "de novo": [f"{de_novo_score_unit}", ""],
                       "de novo only": [f"{de_novo_score_unit}", "d-only"]}

    # assert that stat dict is compatible with db search and de novo
    if sample_db_search is not None and\
        assert_db_search_compatibility(stat_dict, sample_db_search) is False:
        sample_db_search = None
    if sample_de_novo is not None and\
        assert_de_novo_compatibility(stat_dict, sample_de_novo) is False:
        sample_de_novo = None
    
    # check that formats are as expected
    if any([i not in ["db search", "de novo", "de novo only"] for i in formats]):
        raise ValueError("invalid format supplied")
    
    merged_df, _ = reference_score_distribution(stat_dict,
                                                formats,
                                                format_to_label,
                                                normalize_psm,
                                                normalize_rt)

    fig = make_subplots(specs=[[{"secondary_y": True}]])

    # set order of categories (db search, de novo, de novo only)
    cat_order = merged_df["ident method"].unique()

    for category in cat_order:
        cat_data = merged_df[merged_df["ident method"] == category]
        
        # add alc counts to a separate y axis
        if category == "db search confidence dist":
            sec_y = False
        else:
            sec_y = True
        
        # get data from category df
        x_vals = cat_data['x axis'].to_numpy()
        y_vals = cat_data['value'].to_numpy()
        
        fig.add_trace(go.Box(x=x_vals,
                             y=y_vals,
                             boxpoints='all',
                             fillcolor='rgba(255,255,255,0)',
                             line=dict(color='rgba(0,0,0,0)'),
                             marker=dict(color=GraphConstants.color_palette[0]),
                             pointpos=0,
                             #whiskerwidth=0.2,
                             marker_size=6,
                             jitter=1,
                             showlegend=False,
                             #line_width=1,
                             ),
                    secondary_y=sec_y)

    # only process sample if normalization method possible or no normalization
    if not (spectral_metadata is None and any(x is True for x in [normalize_psm,
                                                                  normalize_rt])):
        # based on normalization method, specify division factor within value counts
        div_factor = None
        sample_legend = False
        
        if normalize_psm is True and spectral_metadata is not None:
            div_factor = spectral_metadata.get('spectrum count')
        elif normalize_rt is True and spectral_metadata is not None:
            div_factor = spectral_metadata.get('total retention time')
            

        if sample_db_search is not None and "db search" in formats:
            db_search_data = sample_db_search.data
            lgp_counts, lgp_names = count_threshold_values(db_search_data,
                                                           'Confidence',
                                                           stat_dict['metadata']['db search thresholds'],
                                                           div_factor=div_factor)
            prefix, suffix = format_to_label["db search"]
            lgp_names = [f"{prefix} {i} {suffix}" for i in lgp_names]
            
            fig.add_trace(go.Scatter(
                y=lgp_counts,
                x=lgp_names,
                mode='markers',
                marker_symbol="x",
                marker_line_width=2, 
                marker_size=15,
                name="Sample",
                showlegend=sample_legend == False,
                marker_color=GraphConstants.sample_trace_color
            ))
            sample_legend = True
        
        # Only process sample data if relevant files present
        alc_counts, alc_only_counts = [], []
        alc_names, alc_only_names = [], []
        
        if sample_de_novo is not None and "de novo" in formats:
            de_novo_data = sample_de_novo.data
            alc_counts, alc_names = count_threshold_values(de_novo_data,
                                                           'Confidence',
                                                           stat_dict['metadata']['de novo thresholds'],
                                                           div_factor=div_factor)
            prefix, suffix = format_to_label["de novo"]
            alc_names = [f"{prefix} {i} {suffix}" for i in alc_names]
        
        if sample_db_search is not None and sample_de_novo is not None \
            and "de novo only" in formats:
            sample_de_novo_only = filter_denovo_only(sample_de_novo.data, sample_db_search.data)
        
            alc_only_counts, alc_only_names = count_threshold_values(sample_de_novo_only,
                                                                     'Confidence',
                                                                     stat_dict['metadata']['de novo thresholds'],
                                                                     div_factor=div_factor)
            prefix, suffix = format_to_label["de novo only"]
            alc_only_names = [f"{prefix} {i} {suffix}" for i in alc_only_names]
        
        y_vals = (alc_counts + alc_only_counts)
        x_vals = (alc_names + alc_only_names)

        if len(y_vals) > 0:
            fig.add_trace(go.Scatter(
                y=y_vals,
                x=x_vals,
                mode='markers',
                marker_symbol="x",
                marker_line_width=2, 
                marker_size=15,
                name="Sample",
                showlegend=sample_legend == False,
                marker_color=GraphConstants.sample_trace_color
            ), secondary_y=True)

    # specify y axis names, overwrite based on normalization settings
    alc_y_name = "de novo matches"
    lgp_y_name = "db search matches"
    if normalize_rt is True:
        alc_y_name = "de novo matches / min"
        lgp_y_name = "db search matches / min"
    if normalize_psm is True:
        lgp_y_name = "fraction db search matches"


    fig.update_layout(#title='Distribution of match scores at specific thresholds',
                      margin=dict(l=20, r=30, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      showlegend=True)
    fig.update_xaxes(gridcolor="rgba(0,0,0,0)", 
                     zerolinecolor="Black")
    fig.update_yaxes(title=lgp_y_name, 
                     gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth, 
                    zerolinecolor="Black", nticks=6, 
                    rangemode='tozero')
    fig.update_yaxes(title=alc_y_name, 
                     showgrid=False, 
                     secondary_y=True,
                     tickmode="sync")
    
    return fig


def ref_intensity_dist_plot(stat_dict: dict,
                            spectral_data: None | pd.DataFrame=None):
    # specify prefix string added to categories
    ms1_pref, ms2_pref = "MS1-", "MS2-"
    
    def compute_ms_int_dist(data_dict: dict, mslevel: int):
        if mslevel == 1:
            key_field = "ms1 intensity"
        elif mslevel == 2:
            key_field = "ms2 intensity"
        else:
            raise ValueError("Invalid MS level given, only '1' or '2' supported") 
        
        # store output in dict
        ms_int_dict = dict()

        # fetch data over samples
        for data in data_dict["samples"].values():
            ms_cats = data[key_field]["percentiles"]
            ms_int = data[key_field]["values"]
            
            # add data to output dict
            for cat, value in zip(ms_cats, ms_int):
                if cat not in ms_int_dict.keys():
                    ms_int_dict[cat] = [value]
                else:
                    ms_int_dict[cat].append(value)

        return ms_int_dict

    ms1_dict = compute_ms_int_dist(stat_dict, 1)
    ms2_dict = compute_ms_int_dist(stat_dict, 2)
    
    # if spectral data given, process percentiles of total ion current
    if spectral_data is not None:
        ref_percentiles = stat_dict['metadata']['intensity percentiles']
        ms2_df = spectral_data[spectral_data["MS level"] == 2]
        ms1_df = spectral_data[spectral_data["MS level"] == 1]

        # tuples: (perc. name, TIC value)
        sample_ms2_tuple = array_to_percentiles(
            ms2_df['total ion current'].to_numpy(),
            ref_percentiles)
        sample_ms1_tuple = array_to_percentiles(
            ms1_df['total ion current'].to_numpy(),
            ref_percentiles)
        
        sample_yvals = np.array(sample_ms1_tuple[0] + sample_ms2_tuple[0])
        sample_xvals = np.array(
            [ms1_pref + key for key in sample_ms1_tuple[1]] + \
            [ms2_pref + key for key in sample_ms2_tuple[1]]
        )
    else:
        sample_xvals, sample_yvals = None, None

    # create figure
    fig = go.Figure()


    for key in ms1_dict.keys():
            
        xrow = []
        yrow = []
        
        vals = ms1_dict[key] + ms2_dict[key]
        yrow += vals
        xrow += [ms1_pref + key]*len(ms1_dict[key]) + \
                [ms2_pref + key]*len(ms2_dict[key])
        # yrow += ["MS1"]*len(ms1_dict[key]) + ["MS2"]*len(ms2_dict[key])

        fig.add_trace(go.Box(y=yrow,
                             x=xrow,
                             boxpoints='all',
                             fillcolor='rgba(255,255,255,0)',
                             line=dict(color='rgba(0,0,0,0)'),
                             marker=dict(color=GraphConstants.color_palette[0]),
                             pointpos=0,
                             #whiskerwidth=0.2,
                             marker_size=6,
                             jitter=1,
                             showlegend=False,
                             #line_width=1,
                             ))

        # fig.add_trace(go.Box(
        #     x=xrow,
        #     y=yrow,
        #     boxpoints='all',
        #     #jitter=0.5,
        #     whiskerwidth=0.2,
        #     marker_size=3,
        #     line_width=1,
        #     name=key
        # ))
        
    # add sample data
    if sample_xvals is not None and sample_yvals is not None:
        # fetch values that are part of correct percentile
        #sample_xvals_perc = sample_xvals[sample_p == key]
        #sample_yvals_perc = sample_yvals[sample_p == key]
        
        fig.add_trace(go.Scatter(
            x=sample_xvals,
            y=sample_yvals,
            mode='markers',
            marker_symbol="x",
            marker_color=GraphConstants.sample_trace_color,
            marker_line_width=2, 
            marker_size=15,
            showlegend=True,
            name="Sample"
        ))

    x_order = [ms1_pref + key for key in ms1_dict.keys()] + \
              [ms2_pref + key for key in ms1_dict.keys()]

    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      #yaxis=dict(zeroline=False)
                      )
    
    fig.update_yaxes(title='Intensity', 
                     type="log", 
                     gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth, 
                     nticks=6)
    fig.update_xaxes(showline=True, 
                     linecolor="Black", 
                     linewidth=2, 
                     categoryorder='array',
                     categoryarray=x_order)
    fig.update_legends(title="Percentile")

    # fig.update_traces(orientation='h') # horizontal box plots
    
    return fig


def ref_transmission_scatter_plot(stat_dict: dict,
                                  spectral_data: pd.DataFrame | None=None,
                                  scale_ion_injection_time: bool=True) -> go.Figure:
    # create figure
    fig = go.Figure()

    yrow = []
    xrow = []

    # based on injection time scale, get correct data
    stat_dict_param = "transmission loss" if scale_ion_injection_time is False\
        else "transmission loss ion injection time scaled"

    for data in stat_dict["samples"].values():
        yrow += data[stat_dict_param]['values']
        xrow += data[stat_dict_param]['percentiles']

    fig.add_trace(go.Box(
        y=(1 / np.array(yrow))*100,
        x=xrow,
        boxpoints='all',
        jitter=1,
        fillcolor='rgba(0,0,0,0)',
        line=dict(color='rgba(0,0,0,0)'),
        # marker= dict(color="blue"),
        marker=dict(color=GraphConstants.color_palette[0]),
        pointpos=0,
        whiskerwidth=0.2,
        marker_size=6,
        line_width=1,
        showlegend=False,
        #name=key
    ))
    
    if spectral_data is not None:
        ms2_df = spectral_data[spectral_data["MS level"] == 2]
        ms1_df = spectral_data[spectral_data["MS level"] == 1]
        ref_percentiles = stat_dict['metadata']['intensity percentiles']
        
        if scale_ion_injection_time is True:
            # get precursor ion injection time
            prec_inj_time = fetch_precursor_ion_injection_time(
                ms2_df, ms1_df
            )
            sample_trm_loss = transmission_loss(ms2_df["precursor intensity"],
                                                ms2_df["total ion current"],
                                                prec_inj_time,
                                                ms2_df['ion injection time'])
        else:
            sample_trm_loss = transmission_loss(ms2_df["precursor intensity"],
                                                ms2_df["total ion current"])
        yvals, xvals = array_to_percentiles(sample_trm_loss.to_numpy(),
                                            ref_percentiles,
                                            True)
        
        fig.add_trace(go.Scatter(
            y=(1 / np.array(yvals))*100,
            x=xvals,
            mode='markers',
            marker_symbol="x",
            marker_line_width=2, 
            marker_size=15,
            name='Sample',
            showlegend=True
        ))


    fig.update_layout(
        yaxis=dict(title='transmission efficiency (%)', zeroline=False, type="log"),
                   boxmode='group', )

    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
                     gridwidth=GraphConstants.gridwidth,
                     nticks=5)
    fig.update_xaxes(showline=True, linecolor="Black", linewidth=2)
    fig.update_legends(title="Percentile")
    
    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      barmode="overlay")
    #fig.update_yaxes(autorange='reversed')
    
    return fig

    
def ref_score_threshold_barplot(stat_dict: dict,
                                formats: List[str],
                                normalize_psm: bool=False,
                                normalize_rt: bool=False,
                                normalize_fill: bool=False,
                                filter_de_novo_only: bool=False,
                                sample_db_search: MetaPepDbSearch | None=None,
                                sample_de_novo: MetaPepDeNovo | None=None,
                                spectral_metadata: dict | None=None):
    db_search_score_unit = stat_dict['metadata']['db search confidence format']
    de_novo_score_unit = stat_dict['metadata']['de novo confidence format']
    de_novo_ident_group = "de novo only confidence dist" if filter_de_novo_only is True\
        else "de novo confidence dist"
    
    # Define group label dict, each format has a prefix and suffix
    format_to_label = {"db search": [f"{db_search_score_unit}", ""],
                       "de novo": [f"{de_novo_score_unit}", ""],
                       "de novo only": [f"{de_novo_score_unit}", "d-only"]}
    
    # assert that stat dict is compatible with db search and de novo
    if sample_db_search is not None and\
        assert_db_search_compatibility(stat_dict, sample_db_search) is False:
        sample_db_search = None
    if sample_de_novo is not None and\
        assert_de_novo_compatibility(stat_dict, sample_de_novo) is False:
        sample_de_novo = None

    merged_df, sample_order = reference_score_distribution(
        stat_dict,
        formats,
        format_to_label,
        normalize_psm,
        normalize_rt,
        normalize_fill)
    

    if any(x is not None for x in [sample_db_search, sample_de_novo]):
        ncols = 2
    else:
        ncols = 1

    fig = make_subplots(rows=2,
                        cols=ncols,
                        shared_xaxes=True,
                        specs=[[{}]*ncols,
                               [{}]*ncols],
                        vertical_spacing=0.01,
                        horizontal_spacing=0.05,
                        shared_yaxes=True)

    # match category to bar column
    cat_to_col = dict()
    color_generator = chain(GraphConstants.color_palette)
    def add_traces(fig: go.Figure,
                   dataset: pd.DataFrame,
                   col: int,
                   show_legend: int) -> go.Figure:
        """Wrangle sample threshold data into traces for the distribution
        barplot. Colors are assigned manually and subplot positions are
        defined.

        Args:
            fig (go.Figure): Plotly graph object.
            dataset (pd.DataFrame): Value threshold dataset.
            col (int): Sublot column to add trace to.
            show_legend (int): Add trace to the figure legend.

        Returns:
            go.Figure: Updated plotly graph object
        """
        # sort data to obtain best-to-worst distribution in figure
        db_search_data = dataset[dataset["ident method"] == "db search confidence dist"]\
            .sort_values(by="value", ascending=False)
        de_novo_data =  dataset[dataset["ident method"] == de_novo_ident_group]
        

        # add barplot traces for one dataset
        for groupname, groupdata in db_search_data.groupby("x axis"):
            if f"{groupname}" not in cat_to_col:
                cat_to_col[f"{groupname}"] = color_generator.__next__()

            fig.add_trace(
                go.Bar(
                    name=groupname,
                    x=groupdata["sample"],
                    y=groupdata["value"],
                    marker=dict(color=cat_to_col[f"{groupname}"]),
                    legendgroup=f"{groupname}",
                    showlegend=show_legend
                    ),
                row=1,
                col=col)


        for groupname, groupdata in de_novo_data.groupby("x axis"):
            if f"{groupname}" not in cat_to_col:
                cat_to_col[f"{groupname}"] = color_generator.__next__()
            fig.add_trace(
                go.Bar(
                    name=groupname,
                    x=groupdata["sample"],
                    y=groupdata["value"],
                    marker=dict(color=cat_to_col[f"{groupname}"]),
                    legendgroup=f"{groupname}",
                    showlegend=show_legend
                    ),
                row=2,
                col=col)
        return fig

    fig = add_traces(fig, merged_df, ncols, True)

    # add sample data in separate column
    if ncols == 2:
        if spectral_metadata is None:
            ms2_num, total_rt = None, None
        else:
            ms2_num = spectral_metadata.get('MS2 spectrum count')
            total_rt = spectral_metadata.get('total retention time')
        # extract threshold values from sample data
        sample_data = reference_score_dist_peaks(sample_db_search,
                                                 sample_de_novo,
                                                 ms2_num,
                                                 total_rt,
                                                 formats,
                                                 stat_dict['metadata']['db search thresholds'],
                                                 stat_dict['metadata']['de novo thresholds'],
                                                 normalize_psm,
                                                 normalize_rt,
                                                 normalize_fill)
        fig = add_traces(fig, sample_data, 1, False)
        
        # configure sample figure subplot
        fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
            gridwidth=GraphConstants.gridwidth, 
            zerolinecolor="Black",
            domain=[0.3, 1], row=1, col=1)
        fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
            gridwidth=GraphConstants.gridwidth, 
            showline=False,
            autorange="reversed", 
            domain=[0, 0.3], row=2, col=1)
        fig.update_xaxes(row=1, 
                         col=1, 
                         showticklabels=True,
                         side="top", 
                         tickfont=dict(size=16))
        fig.update_xaxes(row=2, 
                         col=1, 
                         showticklabels=False)
        
        # set domain based on number of datapoints
        ref_sample_num = len(merged_df["sample"].unique())
        sample_range = max(1 / (ref_sample_num + 1), 0.02)
        
        fig.update_xaxes(domain=[0, sample_range], row=1, col=1)
        fig.update_xaxes(domain=[0, sample_range], row=2, col=1)
        fig.update_xaxes(domain=[sample_range + 0.02, 1], row=1, col=2)
        fig.update_xaxes(domain=[sample_range + 0.02, 1], row=2, col=2)
        
    
    fig.update_xaxes(gridcolor="rgba(0,0,0,0)",
                     zerolinecolor="Black", col=ncols,
                     showticklabels=False)

    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
        gridwidth=GraphConstants.gridwidth, 
        zerolinecolor="Black",
        domain=[0.3, 1], 
        row=1, 
        col=ncols)
    fig.update_yaxes(gridcolor=GraphConstants.gridcolor, 
        gridwidth=GraphConstants.gridwidth, 
        showline=False,
        autorange="reversed", 
        domain=[0, 0.3],# tickvals=tickvals,
        row=2, 
        col=ncols)
    fig.update_yaxes(title="DB Search", row=1, col=1)
    fig.update_yaxes(title="ALC", row=2, col=1)

    if sample_order is not None:
        fig.update_xaxes(categoryorder='array',
                         categoryarray=sample_order)

    fig.update_layout(margin=dict(l=20, r=20, t=10, b=10),
                      paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      barmode="overlay")
    
    return fig