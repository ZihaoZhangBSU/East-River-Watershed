"""Reduce existing production daily DEPTH fields and corrected ASO depth fields.

No model runs or remapping occur. ASO fields are the frozen, independently
remapped variable-specific F5 fields with true full-target Cj>=0.5. Depth is
never restricted by SWE validity. Inputs are read-only; outputs are separate.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import pandas as pd
from rasterio.crs import CRS

YEARS=tuple(range(2018,2027))
FLIGHTS=('2018-03-31','2018-05-24','2019-04-07','2019-06-10','2022-04-21',
 '2022-05-18','2023-04-01','2023-05-23','2024-04-03','2024-05-20','2025-04-07',
 '2025-05-20','2026-03-18','2026-04-05')

def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def calendar(wy):
    return pd.date_range(f'{wy-1}-10-01',f'{wy}-06-30' if wy==2026 else f'{wy}-09-30')

def weighted_mean(values,weights,domain,require_full=False):
    """Use finite DEPTH and fractional watershed area, including finite zeros.

    ASO means renormalize over each flight's retained depth targets. Daily
    models require the complete fixed domain; missing model depth aborts.
    No Cj multiplier or SWE-based screen enters either calculation.
    """
    values=np.asarray(values,dtype=float);weights=np.asarray(weights,dtype=float)
    domain=np.asarray(domain,dtype=bool)
    if values.shape!=weights.shape or domain.shape!=weights.shape:
        raise ValueError('Depth, weights and domain shapes must agree')
    if not np.isfinite(weights).all() or np.any((weights<0)|(weights>1)):
        raise ValueError('Weights must be finite fractional watershed areas')
    if not np.array_equal(domain,weights>0):
        raise ValueError('The fixed domain must equal the positive watershed weights')
    valid=domain&np.isfinite(values)
    if require_full and not np.array_equal(valid,domain):
        raise ValueError('Incomplete finite model DEPTH on the fixed watershed footprint')
    if not valid.any():raise ValueError('No retained depth target cells')
    mass=float(np.sum(values[valid]*weights[valid],dtype=np.float64))
    retained=float(np.sum(weights[valid],dtype=np.float64))
    full=float(np.sum(weights[domain],dtype=np.float64))
    mean=mass/retained
    independent=math.fsum(float(v)*float(w) for v,w in zip(values[valid],weights[valid]))/math.fsum(map(float,weights[valid]))
    if abs(mean-independent)>1e-9:raise AssertionError('Independent math.fsum depth reduction differs')
    return mean,retained/full,valid,abs(mean-independent)

def depth_from_corrected(values,coverage,weights,domain):
    """Verify the depth field's OWN full-target coverage; never consult SWE."""
    values=np.asarray(values,dtype=float);coverage=np.asarray(coverage,dtype=float)
    if values.shape!=coverage.shape or values.shape!=weights.shape:
        raise ValueError('ASO depth and depth coverage must match the canonical grid')
    if not np.isfinite(coverage).all() or np.any((coverage<0)|(coverage>1)):
        raise ValueError('ASO depth coverage must be finite and in [0,1]')
    if np.any(np.isfinite(values)&(coverage<.5)):
        raise ValueError('A saved finite ASO depth target violates its own Cj>=0.5 threshold')
    return weighted_mean(values,weights,domain)

def write(frame,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    frame.to_csv(path,index=False,date_format='%Y-%m-%d',float_format='%.17g')

def extract(cache_dir,grid_path,corrected_fields,corrected_manifest,output_dir,
            cache_manifest=None,archived_daily=None):
    cache_dir,grid_path,corrected_fields,corrected_manifest,output_dir=map(Path,
        (cache_dir,grid_path,corrected_fields,corrected_manifest,output_dir))
    with np.load(grid_path,allow_pickle=False) as g:
        weights=np.asarray(g['weights'],float);domain=np.asarray(g['fixed_domain_mask'],bool)
        transform=np.asarray(g['transform'],float);crs=CRS.from_wkt(str(g['crs']))
    if weights.shape!=(66,54) or int(domain.sum())!=842 or not np.array_equal(domain,weights>0):
        raise ValueError('Expected the canonical842-target positive fractional-area Noah grid')
    cell_area=abs(transform[0]*transform[4]-transform[1]*transform[3])
    full_area=float(np.sum(weights[domain],dtype=np.float64)*cell_area)
    if not np.isclose(full_area,748679519.6193271,rtol=0,atol=1e-3):
        raise ValueError('Unexpected fixed fractional watershed area')
    if not crs.is_projected or crs.linear_units.lower() not in ('metre','meter'):
        raise ValueError('Expected projected metre target grid')
    provenance=json.loads(corrected_manifest.read_text(encoding='utf-8'))
    expected_outputs={Path(item['file']).name:item['sha256'] for item in provenance['outputs']}
    if sha256(corrected_fields)!=expected_outputs.get(corrected_fields.name):
        raise ValueError('Corrected ASO field archive does not match its original provenance')
    if sha256(grid_path)!=provenance['grid']['sha256']:
        raise ValueError('Canonical grid does not match corrected ASO provenance')
    sources={item['date']:item for item in provenance['source_manifest'] if item['variable']=='Depth'}
    if tuple(sorted(sources))!=FLIGHTS:raise ValueError('Expected all14 depth-specific ASO source records')
    with np.load(corrected_fields,allow_pickle=False) as z:
        for key in ('weights','fixed_domain_mask','transform'):
            expected={'weights':weights,'fixed_domain_mask':domain,'transform':transform}[key]
            if not np.array_equal(z[key],expected):raise ValueError(f'Corrected field grid mismatch: {key}')
        if CRS.from_wkt(str(z['crs']))!=crs:raise ValueError('Corrected field CRS mismatch')
        select=np.where(z['variable']=='Depth')[0]
        dates=z['date'][select].astype(str)
        if tuple(dates)!=FLIGHTS:raise ValueError('Corrected depth fields must contain the14 chronological flights')
        aso_fields=np.asarray(z['ASO_mm'][select],float)
        aso_coverage=np.asarray(z['full_target_valid_coverage'][select],float)
        saved_models=np.asarray(z['model_mm'][select],float)
        if tuple(z['model_names'].astype(str))!=('iSnobal','Noah-MP'):
            raise ValueError('Unexpected corrected-field model order')
    expected_caches={}
    if cache_manifest:
        expected_caches={Path(item['path']).name:item['sha256'] for item in
            json.loads(Path(cache_manifest).read_text())['spatial_cache_inputs']}
    daily=[];cache_receipts=[];coverage_rows=[];field_lookup={};max_error=0.;max_snapshot_difference=0.
    for wy in YEARS:
        fields={}
        for model in ('isnobal','noahmp'):
            path=cache_dir/f'{model}_wy{wy}_noahgrid_daily.npz';digest=sha256(path)
            if cache_manifest and expected_caches.get(path.name)!=digest:
                raise ValueError(f'Production spatial cache differs from accepted provenance: {path.name}')
            with np.load(path,allow_pickle=False) as z:
                date=pd.DatetimeIndex(pd.to_datetime(z['dates'].astype(str)))
                depth=np.asarray(z['snow_depth_mm'],float)
            if not date.equals(calendar(wy)):raise ValueError(f'Incomplete or noncanonical calendar: {path.name}')
            if depth.shape!=(len(date),*weights.shape):raise ValueError('Depth field dimensions disagree with grid/calendar')
            fields[model]=depth
            finite=np.isfinite(depth[:,domain])
            if not finite.all():raise ValueError(f'Incomplete model depth support: {path.name}')
            coverage_rows.append({'wy':wy,'model':model,'variable':'Depth','n_dates':len(date),
                'fixed_target_count':842,'valid_depth_area_min':1.,'valid_depth_area_max':1.,
                'nonfinite_target_cell_days':int((~finite).sum()),
                'finite_zero_target_cell_days':int((depth[:,domain]==0).sum()),
                'negative_target_cell_days':int((depth[:,domain]<0).sum())})
            cache_receipts.append({'file':path.name,'sha256':digest,'bytes':path.stat().st_size,
                'wy':wy,'model':model,'variable':'snow_depth_mm','n_dates':len(date),
                'shape':list(depth.shape),'source_field_units':'mm'})
        for k,day in enumerate(calendar(wy)):
            means=[]
            for model in ('isnobal','noahmp'):
                mean,fraction,mask,error=weighted_mean(fields[model][k],weights,domain,require_full=True)
                means.append(mean);max_error=max(max_error,error)
            daily.append({'date':day,'wy':wy,'variable':'Depth','depth_i_mm':means[0],
                'depth_n_mm':means[1],'fixed_model_area_fraction':1.,'daily_common_valid_area_fraction':1.,
                'fixed_model_cells':842,'state_hour_utc':23})
            day_string=day.date().isoformat()
            if day_string in FLIGHTS:
                field_lookup[day_string]=(fields['isnobal'][k].copy(),fields['noahmp'][k].copy())
    aso_rows=[];matched=[];depth_field_receipts=[]
    for index,day in enumerate(FLIGHTS):
        value=aso_fields[index];coverage=aso_coverage[index]
        mean,fraction,valid,error=depth_from_corrected(value,coverage,weights,domain)
        max_error=max(max_error,error);source=sources[day]
        aso_rows.append({'date':day,'wy':int(day[:4]),'variable':'Depth',
            'aso_watershed_mean_mm':mean,'aso_valid_area_fraction':fraction,
            'retained_target_count':int(valid.sum()),'fixed_target_count':842,
            'retained_target_watershed_area_m2':full_area*fraction,
            'fixed_model_watershed_area_m2':full_area,
            'source_file':source['source_file'],'source_sha256':source['sha256'],
            'finite_zero_retained_targets':int(np.sum(valid&(value==0))),
            'negative_retained_targets':int(np.sum(valid&(value<0))),
            'maximum_retained_depth_mm':float(np.max(value[valid]))})
        means=[]
        for m,model in enumerate(field_lookup[day]):
            max_snapshot_difference=max(max_snapshot_difference,
                float(np.max(np.abs(model[domain]-saved_models[index,m][domain]))))
            numerator=math.fsum(float(v)*float(w) for v,w in zip(model[valid],weights[valid]))
            means.append(numerator/math.fsum(map(float,weights[valid])))
        matched.append({'date':day,'variable':'Depth','aso_depth_mm':mean,
            'isnobal_depth_mm':means[0],'noahmp_depth_mm':means[1],
            'retained_target_count':int(valid.sum()),'retained_area_fraction':fraction,
            'retained_target_watershed_area_m2':full_area*fraction})
        depth_field_receipts.append({**source,'retained_target_count':int(valid.sum()),
            'depth_retained_area_fraction':fraction,'finite_values_with_Cj_below_half':0})
    if max_snapshot_difference>1e-9:raise ValueError('Flight depth snapshots differ from production daily caches')
    frame=pd.DataFrame(daily)
    archived_check=None
    if archived_daily:
        old=pd.read_csv(archived_daily,parse_dates=['date'])
        if not old.date.equals(frame.date):raise ValueError('Archived common-grid depth dates differ')
        archived_check={model:float(np.max(np.abs(frame[f'depth_{model}_mm']-old[f'depth_{model}_mm'])))
                        for model in ('i','n')}
        if max(archived_check.values())>1e-6:raise ValueError('New depth means disagree with archived common-grid means')
    output_dir.mkdir(parents=True,exist_ok=True)
    write(frame,output_dir/'models/watershed_depth_daily.csv')
    write(pd.DataFrame(coverage_rows),output_dir/'models/depth_variable_coverage.csv')
    write(pd.DataFrame(aso_rows),output_dir/'aso/aso_depth_flights.csv')
    write(pd.DataFrame(matched),output_dir/'aso/matched_depth_flight_means.csv')
    np.savez_compressed(output_dir/'aso/aso_depth_retained_fields.npz',
        date=np.asarray(FLIGHTS),variable=np.asarray('Depth'),depth_mm=aso_fields,
        full_target_valid_coverage=aso_coverage,weights=weights,fixed_domain_mask=domain,
        transform=transform,crs=crs.to_wkt())
    inputs=[grid_path,corrected_fields,corrected_manifest]+([Path(cache_manifest)] if cache_manifest else [])
    report={'status':'FRESH_DEPTH_REDUCTION_VALIDATED','variable':'Depth',
      'daily_records':len(frame),'daily_model_scalar_reductions':len(frame)*2,'ASO_flights':14,
      'no_SWE_mask_or_values_used':True,'fixed_depth_targets':842,'fixed_depth_area_m2':full_area,
      'depth_model_coverage':'Complete finite depth on all842 fixed targets on every day; no daily renormalization.',
      'ASO_operator':'Own depth source mask/NoData, finite zeros/negative/extreme retained; independently remapped to Noah grid; Cj>=0.5 true full-target source-valid coverage. Fractional watershed-area mean over each flight retained depth targets; no Cj multiplier.',
      'coverage_denominator':'Retained depth target watershed area / fixed842-target watershed area; not exact native flight polygon.',
      'periods':'Full October-September WY2018-2025; October1-June30 WY2026; D23 UTC states.',
      'operator_inheritance_limit':'Reuses verified production depth spatial caches and frozen corrected ASO fields. It does not independently reread/remap every original50m model day or original ASO raster. Per-day native-subcell coverage arrays were not saved; finite target depth is not100percent native coverage.',
      'max_independent_weighted_mean_error_mm':max_error,
      'daily_vs_F5_depth_snapshot_max_abs_mm':max_snapshot_difference,
      'archived_common_grid_mean_max_abs_mm':archived_check,
      'spatial_cache_inputs':cache_receipts,'ASO_depth_sources_and_coverage':depth_field_receipts,
      'inputs_sha256':{str(p):sha256(p) for p in inputs},'code_sha256':sha256(Path(__file__)),
      'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__},
      'outputs_sha256':{str(p.relative_to(output_dir)):sha256(p) for p in output_dir.rglob('*') if p.is_file()}}
    (output_dir/'depth_extraction_provenance.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:report[key] for key in ('status','daily_records','ASO_flights',
        'max_independent_weighted_mean_error_mm','daily_vs_F5_depth_snapshot_max_abs_mm',
        'archived_common_grid_mean_max_abs_mm')},indent=2))
    return frame,pd.DataFrame(aso_rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('cache-dir','grid','corrected-fields','corrected-manifest','output-dir'):
        p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--cache-manifest',type=Path)
    p.add_argument('--archived-daily',type=Path)
    a=p.parse_args()
    extract(a.cache_dir,a.grid,a.corrected_fields,a.corrected_manifest,a.output_dir,
            a.cache_manifest,a.archived_daily)

if __name__=='__main__':main()
