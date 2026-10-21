#!/usr/bin/env python3
"""Generate Fig. 5 (SWE) and the matching supplemental ASO depth performance.

Use --scores-csv to replay a complete validated 224-row corrected score table
without reopening rasters. Otherwise recompute from unfiltered raw ASO
SWE/depth and saved model snapshots with the unchanged scientific operator.

The ASO validity coverage operator is zero-padded over the complete target
extent, fixing partial-source-edge normalization. Retain C_j>=0.5 separately
for SWE and depth. Scores retain Fig. 5's equal-cell, model-specific finite-pair
definition; fractional watershed weights define eligibility only (w_j>0).
This differs from Fig. 2's fractional-area basin means. No simulations run.

Requires extract_fig2_common_aso.py beside this script, numpy, rasterio, affine,
and matplotlib. --model-directory contains the saved unfiltered D23 model
rasters named YYYYMMDD_{swe,snow_depth}_{isnobal_on_noah,noah_native}.tif.
Raw input rasters and model products are opened read-only.
"""
from __future__ import annotations
import argparse,csv,json,platform,sys
from pathlib import Path
sys.dont_write_bytecode=True
if Path(__file__).with_name('extract_fig2_common_aso.py').is_file():
    sys.path.insert(0,str(Path(__file__).resolve().parent))
else:
    # In the public figure folder this file lives in code/f5/ and shares the
    # validated remapping module at the figure folder's top level.
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
import rasterio
from affine import Affine
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from extract_fig2_common_aso import load_grid,remap,sha256,units_factor
import extract_fig2_common_aso as remapping

MODELS=[('iSnobal','isnobal_on_noah'),('Noah-MP','noah_native')]
VARIABLES=[('SWE','SWE','swe'),('Depth','SD','snow_depth')]
BANDS={'Lower':(None,2896.),'Middle':(2896.,3353.),'Upper':(3353.,None),'All':(None,None)}
FLIGHT_DATES=('2018-03-31','2018-05-24','2019-04-07','2019-06-10',
 '2022-04-21','2022-05-18','2023-04-01','2023-05-23','2024-04-03',
 '2024-05-20','2025-04-07','2025-05-20','2026-03-18','2026-04-05')
SCORE_COLUMNS=('date','water_year','variable','model','elevation_band','n_cells',
 'aso_mean_mm','mb_mm','rmse_mm','r','rmb_percent','rrmse_percent')


def load_scores_csv(path):
    """Validate the frozen fourteen-flight score schema before plotting.

    This checks completeness and internal consistency, not the source-raster
    operator; independently verify the scores from saved fields separately.
    """
    with path.open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.DictReader(stream)
        if tuple(reader.fieldnames or ())!=SCORE_COLUMNS:
            raise ValueError('Scores CSV must have the canonical 12-column schema')
        rows=list(reader)
    lookup={}
    for row in rows:
        key=tuple(row[k] for k in ('date','variable','model','elevation_band'))
        if key in lookup:
            raise ValueError(f'Duplicate score key: {key}')
        try:
            for k in ('water_year','n_cells'):
                row[k]=int(row[k])
            for k in SCORE_COLUMNS[6:]:
                row[k]=None if row[k] in ('','None') else float(row[k])
        except (TypeError,ValueError) as error:
            raise ValueError(f'Invalid numeric score: {key}') from error
        if row['n_cells']<=0 or row['water_year']!=int(row['date'][:4]):
            raise ValueError(f'Invalid support count or water year: {key}')
        for k in ('aso_mean_mm','mb_mm','rmse_mm'):
            if row[k] is None or not np.isfinite(row[k]):
                raise ValueError(f'Nonfinite required score {k}: {key}')
        if row['rmse_mm']+1e-9<abs(row['mb_mm']):
            raise ValueError(f'RMSE must be at least absolute MB: {key}')
        if row['r'] is not None and (not np.isfinite(row['r']) or abs(row['r'])>1+1e-12):
            raise ValueError(f'Invalid correlation: {key}')
        for relative,metric in [('rmb_percent','mb_mm'),('rrmse_percent','rmse_mm')]:
            expected=row[metric]/row['aso_mean_mm']*100 if row['aso_mean_mm'] else None
            value=row[relative]
            if (expected is None)!=(value is None) or (value is not None and
                    (not np.isfinite(value) or not np.isclose(value,expected,rtol=1e-10,atol=1e-8))):
                raise ValueError(f'Inconsistent relative score {relative}: {key}')
        lookup[key]=row
    expected={(date,var,model,band) for date in FLIGHT_DATES
              for var in ('SWE','Depth') for model,_ in MODELS for band in BANDS}
    if set(lookup)!=expected:
        raise ValueError('Scores CSV must contain exactly all 224 canonical flight/model/variable/band rows')
    for date in FLIGHT_DATES:
        for var in ('SWE','Depth'):
            for model,_ in MODELS:
                parts=[lookup[date,var,model,band] for band in ('Lower','Middle','Upper')]
                total=lookup[date,var,model,'All']
                if sum(r['n_cells'] for r in parts)!=total['n_cells']:
                    raise ValueError('Elevation-band support counts must partition All')
    return rows


def scores(observed,model,domain):
    valid=domain & np.isfinite(observed) & np.isfinite(model)
    o,m=observed[valid],model[valid]
    if not o.size:
        raise ValueError('No finite pairs on requested evaluation footprint')
    delta=m-o
    mean=float(np.mean(o));bias=float(np.mean(delta));rmse=float(np.sqrt(np.mean(delta**2)))
    return {'n_cells':int(o.size),'aso_mean_mm':mean,'mb_mm':bias,'rmse_mm':rmse,
            'r':float(np.corrcoef(o,m)[0,1]) if o.size>1 and np.std(o)>0 and np.std(m)>0 else None,
            'rmb_percent':bias/mean*100 if mean else None,
            'rrmse_percent':rmse/mean*100 if mean else None}


def write_csv(path,rows):
    path.parent.mkdir(exist_ok=True,parents=True)
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def raster(path,shape,transform,crs):
    with rasterio.open(path) as src:
        if src.count!=1 or src.shape!=shape or src.crs!=crs or not src.transform.almost_equals(transform):
            raise ValueError(f'Model raster does not align with canonical target grid: {path}')
        return src.read(1,masked=True).astype('float64').filled(np.nan)


def performance_table(rows,variable,path):
    lookup={(r['date'],r['model']):r for r in rows if r['variable']==variable and r['elevation_band']=='All'}
    lines=[r'\begin{table}[p]',r'\centering',r'\begingroup',r'\setstretch{1.0}',
           r'\fontsize{10.5}{12.6}\selectfont',r'\setlength{\tabcolsep}{4pt}',
           r'\renewcommand{\arraystretch}{1.08}',r'\begin{tabular}{lrrrr}',r'\toprule',
           r'Flight date & \shortstack{ASO mean\\(mm)} & \shortstack{Cells\\I / N} & \shortstack{MB (mm)\\I / N} & \shortstack{RMSE (mm)\\I / N} \\',r'\midrule']
    for date in sorted({key[0] for key in lookup}):
        i,n=lookup[date,'iSnobal'],lookup[date,'Noah-MP']
        # Actual retained model supports coincide; refuse silently reporting one
        # reference mean for distinct supports if a future dataset changes this.
        if i['n_cells']!=n['n_cells'] or abs(i['aso_mean_mm']-n['aso_mean_mm'])>1e-8:
            raise ValueError('Table needs separate ASO means for unequal model supports')
        lines.append(f"{date} & {i['aso_mean_mm']:.1f} & {i['n_cells']} / {n['n_cells']} & {i['mb_mm']:+.1f} / {n['mb_mm']:+.1f} & {i['rmse_mm']:.1f} / {n['rmse_mm']:.1f} \\\\")
    noun='SWE' if variable=='SWE' else 'snow-depth'
    lines.extend([r'\bottomrule',r'\end{tabular}',r'\endgroup',
      r'\caption{Flight-specific '+noun+r' mean bias (MB) and RMSE relative to unfiltered ASO. Model pairs are iSnobal / Noah-MP; MB is model minus ASO and RMSE is nonnegative. Values are newly recomputed and displayed to 0.1~mm; full precision is retained in the accompanying CSV. ASO SWE and depth are remapped independently, retaining target cells with at least 50\% valid source coverage over the full target-cell footprint, including source-exterior space as missing. Scores use equally weighted, model-specific finite ASO--model pairs; the two supports happen to coincide on all fourteen flights in these saved model fields. These equal-cell reference means differ from the fractionally weighted Fig.~2 means and the fixed paired-flight WY2026 endpoint means. Model states are at 2300~UTC on nominal flight dates, rather than exact flight-hour matching. ASO SWE incorporates modeled density, including iSnobal in documented production; exact density-production details for the ten later flights remain unresolved.}',
      r'\label{tab:supp_aso_'+('swe' if variable=='SWE' else 'depth')+r'_performance}',r'\end{table}'])
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def plot(rows,figure_directory,font_path,figure_selection='both'):
    """Plot per-flight spatial scores without pooling dates.

    SWE contains MB, RMSE and spatial Pearson r. Depth retains its accepted
    two-panel layout. Select ``swe`` to avoid rewriting a frozen depth asset.
    An unavailable r is never plotted as zero; its row is explicitly marked NA.
    """
    if figure_selection not in ('both','swe','depth'):
        raise ValueError('figure_selection must be both, swe, or depth')
    if font_path:
        font_manager.fontManager.addfont(str(font_path))
        font=font_manager.FontProperties(fname=font_path).get_name()
    else:
        font='DejaVu Sans'
    plt.rcParams.update({'font.family':font,'font.size':8,'axes.labelsize':8,'axes.titlesize':8,
      'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,
      'axes.linewidth':.6,'lines.linewidth':1.1,'savefig.facecolor':'white'})
    lookup={(r['date'],r['variable'],r['model']):r for r in rows if r['elevation_band']=='All'}
    dates=sorted({r['date'] for r in rows})
    colors={'iSnobal':'#0072B2','Noah-MP':'#D55E00'}
    audit=[]
    figure_directory.mkdir(exist_ok=True,parents=True)
    for variable_number,(var,basename,ref) in enumerate([
            ('SWE','JHM_F5_aso_performance',300),
            ('Depth','JHM_S_ASO_depth_performance',1000)]):
        if figure_selection!='both' and var.lower()!=figure_selection:
            continue
        metrics=('mb_mm','rmse_mm','r') if var=='SWE' else ('mb_mm','rmse_mm')
        titles=(f'{var} bias',f'{var} RMSE','Spatial r') if var=='SWE' else (f'{var} bias',f'{var} RMSE')
        fig,axs=plt.subplots(1,len(metrics),figsize=(6.5,5.3),sharey=True)
        fig.subplots_adjust(left=.17,right=.985,top=.91,bottom=.18,wspace=.23)
        for local_panel,(ax,metric,title) in enumerate(zip(axs,metrics,titles)):
            for y,date in enumerate(dates):
                i,n=lookup[date,var,'iSnobal'],lookup[date,var,'Noah-MP']
                if date.startswith('2026'):
                    ax.axhspan(y-.5,y+.5,color='#e9b94a',alpha=.25,lw=0)
                vals=[i[metric],n[metric]]
                if all(v is not None for v in vals):
                    ax.plot(vals,[y,y],color='.65',lw=.7)
                if i['n_cells']!=n['n_cells'] or abs(i['aso_mean_mm']-n['aso_mean_mm'])>1e-8:
                    raise ValueError('Matched marker sizing requires common ASO reference support')
                observed=i['aso_mean_mm'];area=16*observed/ref
                if vals[1] is not None:
                    ax.scatter(vals[1],y,s=area,marker='s',color=colors['Noah-MP'],edgecolor='white',linewidth=.3,zorder=3)
                if vals[0] is not None:
                    ax.scatter(vals[0],y,s=4*area/np.pi,marker='o',facecolor='none',edgecolor=colors['iSnobal'],linewidth=.8,zorder=4)
                if any(v is None for v in vals):
                    missing='/'.join(m for m,v in zip(('I','N'),vals) if v is None)
                    ax.text(.97,y,f'{missing}: NA',transform=ax.get_yaxis_transform(),
                            ha='right',va='center',fontsize=8)
                audit.append({'panel':(0 if var=='SWE' else 3)+local_panel+1,'date':date,'variable':var,'metric':metric,
                              'isnobal_x':vals[0],'noahmp_x':vals[1],'shared_y':y,
                              'aso_reference_mean_mm':observed,'enclosed_marker_area_pt2':area})
            if metric!='r':
                ax.axvline(0,color='.5',lw=.6)
            else:
                values=[lookup[date,var,model]['r'] for date in dates for model,_ in MODELS
                        if lookup[date,var,model]['r'] is not None]
                # Show the observed range, explicitly labelled by axis ticks.
                # Accommodate finite negative and exactly zero correlations.
                lo,hi=(min(values),max(values)) if values else (-1.,1.)
                padding=max(.04,.1*(hi-lo))
                ax.set_xlim(max(-1.,lo-padding),min(1.,hi+padding))
            ax.set_ylim(len(dates)-.4,-.8);ax.set_title(title);ax.set_xlabel('Pearson r' if metric=='r' else 'mm')
            ax.set_yticks(range(len(dates)),dates);ax.tick_params(axis='y',length=0)
            ax.xaxis.set_major_locator(MaxNLocator(3));ax.spines[['top','right','left']].set_visible(False)
            ax.text(.02,.98,f'({chr(97+local_panel)})',transform=ax.transAxes,ha='left',va='top',fontsize=8,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':.5})
        fig.legend(handles=[Line2D([],[],ls='',marker=mk,color=colors[m],mfc='none' if m=='iSnobal' else colors[m],label=m,ms=4)
                            for m,mk in [('iSnobal','o'),('Noah-MP','s')]],ncol=2,loc='upper center',frameon=False)
        h=[Line2D([],[],ls='',marker='o',mfc='none',mec='.4',ms=np.sqrt(4*(16*v/ref)/np.pi),label=str(v))
           for v in ([100,300,600] if var=='SWE' else [300,1000,2000])]
        fig.legend(handles=h,title=f'ASO {var.lower()} (mm)',loc='lower center',bbox_to_anchor=(.58,.002),ncol=3,
                   frameon=False,fontsize=8,title_fontsize=8,handletextpad=.2,columnspacing=.6)
        for extension in ['pdf','png']:
            fig.savefig(figure_directory/f'{basename}.{extension}',dpi=600,facecolor='white')
        plt.close(fig)
    return audit


def generated_figures(directory,figure_selection):
    """List only the assets generated by the requested presentation branch."""
    names={'swe':'JHM_F5_aso_performance','depth':'JHM_S_ASO_depth_performance'}
    selected=names.values() if figure_selection=='both' else (names[figure_selection],)
    return sorted(directory/f'{name}.{extension}' for name in selected
                  for extension in ('pdf','png') if (directory/f'{name}.{extension}').is_file())


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ['aso-directory','model-directory','flights','grid-weights']:
        parser.add_argument('--'+name,type=Path)
    for name in ['output-directory','figure-directory']:
        parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--scores-csv',type=Path,
                        help='Replay the complete corrected 224-row score table; do not read rasters')
    parser.add_argument('--baseline-displayed-csv',type=Path)
    parser.add_argument('--font',type=Path)
    parser.add_argument('--figure-selection',choices=('both','swe','depth'),default='both',
                        help='Plot both figures (default), only main SWE, or only supplemental depth')
    a=parser.parse_args()
    if a.scores_csv:
        if any((a.aso_directory,a.model_directory,a.flights,a.grid_weights,a.baseline_displayed_csv)):
            parser.error('--scores-csv cannot be combined with raw-raster inputs or baseline comparison')
        rows=load_scores_csv(a.scores_csv)
        out=a.output_directory;out.mkdir(exist_ok=True,parents=True)
        # Retain all full-precision rows. Replay changes presentation only.
        write_csv(out/'F5_elevation_metrics_corrected.csv',rows)
        write_csv(out/'F5_all_flight_metrics_corrected.csv',[r for r in rows if r['elevation_band']=='All'])
        audit=plot(rows,a.figure_directory,a.font,a.figure_selection)
        write_csv(out/'F5_marker_placement_audit.csv',audit)
        manifest={'mode':'validated_score_replay','input_scores':{'file':a.scores_csv.name,'sha256':sha256(a.scores_csv)},
          'score_rows':len(rows),'plotted_markers':2*len(audit),
          'presentation':'Fig.5 SWE bias/RMSE/spatial Pearson r; supplemental depth bias/RMSE. Scientific scores unchanged.',
          'figure_selection':a.figure_selection,
          'code_sha256':sha256(Path(__file__)),
          'font':({'file':a.font.name,'sha256':sha256(a.font)} if a.font else
                  {'family':'DejaVu Sans','source':'Matplotlib bundled default'}),
          'versions':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__},
          'figures':[{'file':p.name,'sha256':sha256(p)} for p in generated_figures(a.figure_directory,a.figure_selection)]}
        (out/'F5_split_replay_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        return
    if not all((a.aso_directory,a.model_directory,a.flights,a.grid_weights)):
        parser.error('Raw recomputation requires --aso-directory, --model-directory, --flights and --grid-weights')
    weights,domain,transform,crs=load_grid(a.grid_weights)
    with np.load(a.grid_weights) as z:
        dem=np.asarray(z['dem'],float)
    if not np.isfinite(dem[domain]).all():
        raise ValueError('Canonical DEM must be finite over every watershed target')
    flights=list(csv.DictReader(a.flights.open(encoding='utf-8-sig')))
    if not flights or len({f['date'] for f in flights})!=len(flights):
        raise ValueError('Flight dates must be nonempty and unique')
    rows=[];legacy_rows=[];coverages=[];sources=[];fields=[];obsfields=[];coveragefields=[]
    for flight in sorted(flights,key=lambda f:f['date']):
        day=flight['date'];stamp=day.replace('-','')
        for variable,token,slug in VARIABLES:
            source=a.aso_directory/f'ASO_50M_{token}_USCOGE_{stamp}.tif'
            with rasterio.open(source) as src:
                factor,evidence=units_factor(src,'auto')
                raw=(src.read(1,masked=True).astype(float).filled(np.nan)*src.scales[0]+src.offsets[0])*factor
                observed,C,oldC,unthresholded,padding=remap(raw,src.transform,src.crs,transform,crs,weights.shape)
                sourceinfo={'date':day,'variable':variable,'source_file':source.name,'sha256':sha256(source),
                    'source_shape':[src.height,src.width],'source_transform':list(src.transform),
                    'source_crs_wkt':src.crs.to_wkt(),'source_nodata':str(src.nodata),
                    'scale':src.scales[0],'offset':src.offsets[0],'units':evidence,
                    'mask_flags':[[flag.name for flag in flags] for flags in src.mask_flag_enums],
                    'zero_count':int(np.sum(raw==0)),'negative_count':int(np.sum(raw<0)),
                    'finite_min_mm':float(np.nanmin(raw)),'finite_max_mm':float(np.nanmax(raw)),
                    'padding_pixels':padding,'models':[]}
            oldobs=unthresholded.copy();oldobs[oldC<.5]=np.nan
            oldvalid=domain&np.isfinite(oldobs);newvalid=domain&np.isfinite(observed)
            changed=domain&(oldvalid!=newvalid)
            coverage={'date':day,'variable':variable,'old_target_count':int(oldvalid.sum()),
                'new_target_count':int(newvalid.sum()),'changed_target_count':int(changed.sum()),
                'old_retained_target_area_fraction':float(weights[oldvalid].sum(dtype=float)/weights[domain].sum(dtype=float)),
                'new_retained_target_area_fraction':float(weights[newvalid].sum(dtype=float)/weights[domain].sum(dtype=float)),
                'max_Cj_correction':float(np.max(abs(C[domain]-oldC[domain])))}
            sourceinfo['threshold_changed_targets']=[{'row':int(r),'col':int(c),'dem_m':float(dem[r,c]),
                'watershed_fraction':float(weights[r,c]),'legacy_Cj':float(oldC[r,c]),'corrected_Cj':float(C[r,c]),
                'value_mm':float(unthresholded[r,c])} for r,c in zip(*np.nonzero(changed))]
            model_fields=[]
            for model,modelslug in MODELS:
                path=a.model_directory/f'{stamp}_{slug}_{modelslug}.tif'
                m=raster(path,weights.shape,transform,crs);model_fields.append(m)
                sourceinfo['models'].append({'model':model,'file':path.name,'sha256':sha256(path),
                    'valid_basin_target_count':int(np.sum(domain&np.isfinite(m)))})
                for band,(lo,hi) in BANDS.items():
                    bandmask=domain.copy()
                    if lo is not None:bandmask &= dem>=lo
                    if hi is not None:bandmask &= dem<hi
                    prefix={'date':day,'water_year':int(day[:4]),'variable':variable,'model':model,'elevation_band':band}
                    rows.append({**prefix,**scores(observed,m,bandmask)})
                    legacy_rows.append({**prefix,**scores(oldobs,m,bandmask)})
            sources.append(sourceinfo);coverages.append(coverage);fields.append(model_fields)
            obsfields.append(observed);coveragefields.append(C)
            print(day,variable,'targets',coverage['old_target_count'],'->',coverage['new_target_count'],flush=True)
    out=a.output_directory;out.mkdir(exist_ok=True,parents=True)
    allrows=[r for r in rows if r['elevation_band']=='All']
    write_csv(out/'F5_all_flight_metrics_corrected.csv',allrows)
    write_csv(out/'F5_elevation_metrics_corrected.csv',rows)
    write_csv(out/'F5_metrics_legacy_unpadded.csv',legacy_rows)
    write_csv(out/'F5_coverage_boundary_audit.csv',coverages)
    deltas=[]
    for old,new in zip(legacy_rows,rows):
        delta={key:new[key] for key in ['date','water_year','variable','model','elevation_band']}
        for metric in ['n_cells','aso_mean_mm','mb_mm','rmse_mm','r','rmb_percent','rrmse_percent']:
            delta['old_'+metric]=old[metric];delta['new_'+metric]=new[metric]
            delta['change_'+metric]=new[metric]-old[metric] if new[metric] is not None and old[metric] is not None else None
        if any(abs(delta['change_'+k])>1e-10 for k in ['n_cells','aso_mean_mm','mb_mm','rmse_mm']):
            deltas.append(delta)
    write_csv(out/'F5_changed_metrics_full_precision.csv',deltas)
    if a.baseline_displayed_csv:
        prior=list(csv.DictReader(a.baseline_displayed_csv.open(encoding='utf-8-sig')))
        oldlookup={(r['date'],r['variable'],r['model']):r for r in prior}
        display_changes=[]
        for r in allrows:
            old=oldlookup[r['date'],r['variable'],r['model']]
            for key in ['aso_mean_mm','mb_mm','rmse_mm']:
                if f"{r[key]:.1f}" != f"{float(old[key]):.1f}":
                    display_changes.append({'date':r['date'],'variable':r['variable'],'model':r['model'],
                        'metric':key,'previous_displayed':float(old[key]),'new_displayed':round(r[key],1),
                        'new_full_precision':r[key]})
        write_csv(out/'F5_old_displayed_vs_new.csv',display_changes)
        # Keep source-version/display refresh separate from the mathematical
        # boundary-coverage fix. This comparison still uses legacy unpadded C_j.
        refresh_changes=[]
        for r in legacy_rows:
            if r['elevation_band']!='All':continue
            old=oldlookup[r['date'],r['variable'],r['model']]
            for key in ['aso_mean_mm','mb_mm','rmse_mm']:
                if f"{r[key]:.1f}" != f"{float(old[key]):.1f}":
                    refresh_changes.append({'date':r['date'],'variable':r['variable'],'model':r['model'],
                        'metric':key,'previous_displayed':float(old[key]),
                        'fresh_unpadded_displayed':round(r[key],1),'fresh_unpadded_full_precision':r[key]})
        write_csv(out/'F5_stale_display_vs_fresh_unpadded.csv',refresh_changes)
    for var in ['SWE','Depth']:
        performance_table(rows,var,out/f'JHM_TS_ASO_{var.lower()}_performance_20261019.tex')
    audit=plot(rows,a.figure_directory,a.font,a.figure_selection)
    write_csv(out/'F5_marker_placement_audit.csv',audit)
    counts={}
    for variable in ['SWE','Depth']:
        lookup={(r['date'],r['model']):r for r in allrows if r['variable']==variable}
        dates=sorted({key[0] for key in lookup})
        wins=[date for date in dates if lookup[date,'iSnobal']['rmse_mm']<lookup[date,'Noah-MP']['rmse_mm']]
        counts[variable]={'iSnobal_lower_RMSE':len(wins),'total_flights':len(dates),'iSnobal_win_dates':wins,
            'Noah_MP_lower_RMSE_dates':[date for date in dates if lookup[date,'iSnobal']['rmse_mm']>lookup[date,'Noah-MP']['rmse_mm']]}
    np.savez_compressed(out/'F5_corrected_fields.npz',date=np.array([s['date'] for s in sources]),
        variable=np.array([s['variable'] for s in sources]),ASO_mm=np.array(obsfields),
        model_mm=np.array(fields),model_names=np.array([m[0] for m in MODELS]),
        full_target_valid_coverage=np.array(coveragefields),weights=weights,dem=dem,
        fixed_domain_mask=domain,transform=np.array(tuple(transform)),crs=crs.to_wkt())
    summary={'method':'Raw unfiltered ASO independently averaged to Noah grid for SWE/depth; zero-padded validity C_j>=.5 over full target; equal-cell model-specific finite-pair scores.',
      'score_vs_Fig2':'Fig5 uses equal-cell spatial scores; Fig2 uses fractional watershed-area means. Eligibility is w_j>0 for both. Their reference means and footprints are explicitly different.',
      'coverage_correction':'Only ASO C_j denominator and resulting retained support are corrected; SWE/depth remap values and model fields remain unchanged.',
      'display_refresh_distinction':'F5_stale_display_vs_fresh_unpadded.csv compares prior displayed values with fresh raw source metrics under the OLD unpadded operator. F5_changed_metrics_full_precision.csv isolates ONLY the boundary coverage correction, fresh raw data held fixed.',
      'headline_win_counts':counts,'source_manifest':sources,
      'grid':{'file':a.grid_weights.name,'sha256':sha256(a.grid_weights)},
      'code_sha256':sha256(Path(__file__)),'remap_code_sha256':sha256(Path(remapping.__file__)),
      'baseline_displayed_csv':({'file':a.baseline_displayed_csv.name,'sha256':sha256(a.baseline_displayed_csv)} if a.baseline_displayed_csv else None),
      'versions':{'python':platform.python_version(),'numpy':np.__version__,'rasterio':rasterio.__version__,'GDAL':rasterio.__gdal_version__,'matplotlib':matplotlib.__version__},
      'outputs':[{'file':str(p.relative_to(out)),'sha256':sha256(p)} for p in sorted(out.glob('*')) if p.is_file() and p.suffix in ['.csv','.npz','.tex']],
      'figures':[{'file':p.name,'sha256':sha256(p)} for p in generated_figures(a.figure_directory,a.figure_selection)]}
    (out/'F5_correction_manifest.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(counts,indent=2))

if __name__=='__main__':
    main()
