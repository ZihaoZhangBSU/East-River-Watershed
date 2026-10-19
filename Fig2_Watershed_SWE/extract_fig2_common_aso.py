#!/usr/bin/env python3
"""Remap unfiltered ASO SWE to the Noah-MP grid for Fig. 2.

Inputs are read-only. ASO source masks and finite values define validity: valid
zeros, negative values and extremes are retained. GDAL area-overlap averaging
normalizes SWE over valid contributing source area, as in the Fig. 5 operator.
The validity raster is zero-padded to cover the entire target extent BEFORE
averaging, so C_j denotes valid source area / FULL target-cell area even where a
source raster ends inside a target cell. Cells require C_j >= 0.5.

After remapping, Fig. 2 uses a_j = fractional watershed overlap * target-cell
planimetric area, restricted to the fixed model-domain mask. Each flight's mean
is sum(a_j*SWE_j)/sum(a_j) over retained cells. It does not multiply a_j by C_j.
Thus reported coverage is represented target watershed area / fixed model
watershed area, not the exact native valid flight footprint. This two-stage
operator can include outside-watershed source pixels in a boundary target-cell
mean and is intentionally different from clipping source pixels first.

Example:
 python extract_fig2_common_aso.py --aso-directory ASO --flights flights.csv \
   --grid-weights grid_common.npz --output-directory data/aso

The grid NPZ requires weights, transform (Affine coefficients), crs (WKT) and
fixed_domain_mask. A file containing weights alone is refused. The flights CSV
requires date and optionally source_file (relative to --aso-directory).
Dependencies: numpy, rasterio, affine. No source observations are distributed.
"""
from __future__ import annotations
import argparse
import csv
from datetime import date
import hashlib
import json
from pathlib import Path
import platform
import re
import sys
sys.dont_write_bytecode = True
import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.warp import reproject, Resampling, transform_bounds


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load_grid(path):
    with np.load(path, allow_pickle=False) as data:
        required = {'weights', 'transform', 'crs', 'fixed_domain_mask'}
        if not required.issubset(data.files):
            raise ValueError(f'Grid file must include {sorted(required)}')
        weights = np.array(data['weights'], dtype=float)
        mask = np.array(data['fixed_domain_mask'], dtype=bool)
        transform = Affine(*np.asarray(data['transform'], dtype=float).ravel())
        crs = CRS.from_wkt(str(data['crs'].item()))
        for name, axis in [('x', 1), ('y', 0)]:
            if name in data.files:
                centers = (transform.c + (np.arange(weights.shape[1])+.5)*transform.a
                           if name == 'x' else
                           transform.f + (np.arange(weights.shape[0])+.5)*transform.e)
                if not np.allclose(data[name], centers, atol=1e-6, rtol=0):
                    raise ValueError(f'Grid {name} coordinates disagree with affine')
    if not crs.is_projected or not weights.ndim == 2 or mask.shape != weights.shape:
        raise ValueError('Projected 2-D grid and same-shape fixed mask are required')
    if transform.b != 0 or transform.d != 0 or transform.a == 0 or transform.e == 0:
        raise ValueError('Rotated or degenerate target grids are unsupported')
    if not np.isfinite(weights).all() or np.any((weights < 0) | (weights > 1)):
        raise ValueError('Watershed fractions must lie in [0,1] and be finite')
    if not np.any(mask & (weights > 0)):
        raise ValueError('Fixed model footprint has no watershed area')
    if np.any(mask & (weights <= 0)):
        raise ValueError('Fixed footprint may include only positive watershed fractions')
    # The production grid uses projected metres; refuse undeclared conversions.
    if str(crs.linear_units).lower() not in {'metre', 'meter'}:
        raise ValueError('Target projected coordinate units must be metres')
    return weights, mask, transform, crs


def grid_bounds(transform, shape):
    x0,y0 = transform * (0,0)
    x1,y1 = transform * (shape[1],shape[0])
    return min(x0,x1),min(y0,y1),max(x0,x1),max(y0,y1)


def full_target_validity(valid, source_transform, source_crs,
                         target_transform, target_crs, target_shape):
    """Average a 0/1 validity raster with outside-source space explicitly zero.

    GDAL average without padding normalizes over the source-extent intersection,
    overstating C_j at a partially covered target cell. Padding must extend over
    the full target envelope in source coordinates, rather than one source pixel.
    """
    if source_transform.b or source_transform.d or source_transform.a <= 0 or source_transform.e >= 0:
        raise ValueError('ASO validity requires an unrotated north-up source')
    xmin,ymin,xmax,ymax = transform_bounds(target_crs, source_crs,
        *grid_bounds(target_transform, target_shape), densify_pts=41)
    cols = [(x-source_transform.c)/source_transform.a for x in (xmin,xmax)]
    rows = [(y-source_transform.f)/source_transform.e for y in (ymin,ymax)]
    left=max(0, -int(np.floor(min(cols)))+2)
    right=max(0, int(np.ceil(max(cols)))-valid.shape[1]+2)
    top=max(0, -int(np.floor(min(rows)))+2)
    bottom=max(0, int(np.ceil(max(rows)))-valid.shape[0]+2)
    expanded_shape=(valid.shape[0]+top+bottom,valid.shape[1]+left+right)
    if np.prod(expanded_shape) > 100_000_000:
        raise ValueError('Target envelope requires excessive source padding; check grid CRS and extent')
    padded=np.pad(valid.astype('float32'), ((top,bottom),(left,right)), constant_values=0)
    padded_transform=source_transform*Affine.translation(-left,-top)
    coverage=np.zeros(target_shape,dtype='float32')
    reproject(padded,coverage,src_transform=padded_transform,src_crs=source_crs,
              dst_transform=target_transform,dst_crs=target_crs,src_nodata=None,
              dst_nodata=0.0,resampling=Resampling.average,init_dest_nodata=True)
    return np.clip(coverage,0,1), {'left':left,'right':right,'top':top,'bottom':bottom}


def remap(values, source_transform, source_crs, target_transform, target_crs,
          target_shape, minimum_valid_coverage=.5):
    """Return SWE, true full-target C_j, legacy F5 C_j and padding audit."""
    if not 0 <= minimum_valid_coverage <= 1:
        raise ValueError('Coverage threshold must lie in [0,1]')
    data=np.asarray(values,dtype='float32')
    valid=np.isfinite(data)
    sentinel=np.float32(-3.4e38)
    if np.any(data[valid] == sentinel):
        raise ValueError('Finite data collide with the GDAL NoData sentinel')
    destination=np.full(target_shape,sentinel,dtype='float32')
    legacy=np.zeros(target_shape,dtype='float32')
    shared=dict(src_transform=source_transform,src_crs=source_crs,
                dst_transform=target_transform,dst_crs=target_crs,
                resampling=Resampling.average,init_dest_nodata=True)
    reproject(np.where(valid,data,sentinel),destination,src_nodata=sentinel,
              dst_nodata=sentinel,**shared)
    reproject(valid.astype('float32'),legacy,src_nodata=None,dst_nodata=0.0,**shared)
    coverage,padding=full_target_validity(valid,source_transform,source_crs,
                                        target_transform,target_crs,target_shape)
    unthresholded=destination.astype(float)
    unthresholded[destination == sentinel]=np.nan
    retained=unthresholded.copy()
    retained[coverage < minimum_valid_coverage]=np.nan
    return retained,coverage,legacy,unthresholded,padding


def basin_mean(values, weights, fixed_mask):
    """Fractional-area mean; no equal-cell weighting or C_j multiplier."""
    if values.shape != weights.shape or fixed_mask.shape != weights.shape:
        raise ValueError('All target fields must have the same shape')
    valid=fixed_mask & (weights>0) & np.isfinite(values)
    full=float(np.sum(weights[fixed_mask],dtype='float64'))
    retained=float(np.sum(weights[valid],dtype='float64'))
    if full <= 0 or retained <= 0:
        raise ValueError('No retained positive watershed area')
    return float(np.sum(values[valid]*weights[valid],dtype='float64')/retained),retained/full,valid


def units_factor(src, choice):
    if choice != 'auto':
        return (1000. if choice=='m' else 1.),f'{choice} explicit'
    text=' '.join([str(src.units[0]),str(src.descriptions[0]),*map(str,src.tags().values()),
                   *map(str,src.tags(1).values())]).lower()
    mm=bool(re.search(r'\bmm\b|\bmillimet(?:er|re)s?\b',text))
    metre=bool(re.search(r'\bm\b|\bmet(?:er|re)s?\b',text))
    if mm == metre:
        raise ValueError(f'Source SWE units absent or ambiguous: {src.name}')
    return (1.,'mm source metadata') if mm else (1000.,'m source metadata')


def extract(aso_directory, flights_path, grid_path, output_directory, swe_units='auto',
            native_comparison=None, archived_regrid_directory=None):
    aso_directory,flights_path,grid_path,output_directory=map(Path,
        (aso_directory,flights_path,grid_path,output_directory))
    weights,mask,transform,crs=load_grid(grid_path)
    with flights_path.open(newline='',encoding='utf-8-sig') as stream:
        flights=list(csv.DictReader(stream))
    if not flights or len({f['date'] for f in flights}) != len(flights):
        raise ValueError('Flight dates must be nonempty and unique')
    flights=sorted(flights,key=lambda f:date.fromisoformat(f['date']))
    sources=[aso_directory/(f.get('source_file') or f"ASO_50M_SWE_USCOGE_{f['date'].replace('-','')}.tif") for f in flights]
    protected={p.resolve() for p in [flights_path,grid_path,*sources,Path(__file__)]}
    output_names=['aso_common_grid_flights.csv','aso_common_grid_14_flights.npz','aso_provenance.json']
    if any((output_directory/name).resolve() in protected for name in output_names):
        raise ValueError('An output collides with a source input')
    previous={}
    if native_comparison:
        with Path(native_comparison).open(newline='',encoding='utf-8-sig') as stream:
            previous={f['date']:f for f in csv.DictReader(stream)}
    rows=[];details=[];fields=[];coverages=[];legacy_coverages=[];raw_fields=[]
    cell_area=abs(transform.a*transform.e)
    for flight,source in zip(flights,sources):
        day=date.fromisoformat(flight['date']); digest=sha256(source)
        with rasterio.open(source) as src:
            if src.count != 1 or not src.crs:
                raise ValueError('Single-band source with a declared CRS is required')
            multiplier,unit_evidence=units_factor(src,swe_units)
            values=src.read(1,masked=True).astype('float64').filled(np.nan)
            values=(values*src.scales[0]+src.offsets[0])*multiplier
            value,coverage,legacy,unthresholded,padding=remap(values,src.transform,src.crs,
                                                            transform,crs,weights.shape)
            mean,fraction,retained=basin_mean(value,weights,mask)
            legacy_field=unthresholded.copy();legacy_field[legacy<.5]=np.nan
            legacy_mean,legacy_fraction,legacy_retained=basin_mean(legacy_field,weights,mask)
            changed=mask & (legacy_retained != retained)
            source_details={'date':day.isoformat(),'source_file':source.name,'sha256':digest,
                'crs_wkt':src.crs.to_wkt(),'transform':list(src.transform),'shape':[src.height,src.width],
                'nodata':str(src.nodata),'scale':src.scales[0],'offset':src.offsets[0],
                'mask_flags':[[f.name for f in flags] for flags in src.mask_flag_enums],
                'unit_evidence':unit_evidence,'tags':src.tags(),'padding_pixels':padding,
                'finite_zero_count':int(np.sum(values==0)),
                'finite_negative_count':int(np.sum(values<0)),
                'finite_min_mm':float(np.nanmin(values)),'finite_max_mm':float(np.nanmax(values)),
                'coverage_changed_target_cells':[
                    {'row':int(r),'col':int(c),'watershed_fraction':float(weights[r,c]),
                     'legacy_Cj':float(legacy[r,c]),'full_target_Cj':float(coverage[r,c]),
                     'SWE_mm':float(unthresholded[r,c])}
                    for r,c in zip(*np.nonzero(changed))]}
        row={'date':day.isoformat(),'wy':day.year+(day.month>=10),
             'aso_watershed_mean_mm':mean,'aso_valid_area_fraction':fraction,
             'retained_target_count':int(retained.sum()),'fixed_target_count':int(mask.sum()),
             'retained_target_watershed_area_m2':float(weights[retained].sum(dtype='float64')*cell_area),
             'fixed_model_watershed_area_m2':float(weights[mask].sum(dtype='float64')*cell_area),
             'valid_overlap_weighted_target_proxy_fraction':float(np.sum(weights[mask]*coverage[mask],dtype='float64')/np.sum(weights[mask],dtype='float64')),
             'legacy_unpadded_area_weighted_mean_mm':legacy_mean,
             'legacy_unpadded_retained_area_fraction':legacy_fraction,
             'coverage_threshold_changed_cells':int(changed.sum()),
             'max_Cj_boundary_correction':float(np.max(abs(legacy[mask]-coverage[mask]))),
             'source_file':source.name,'source_sha256':digest}
        if previous:
            old=previous[day.isoformat()]
            if old['source_sha256'] != digest:
                raise ValueError('Prior native mean uses a different source hash')
            row.update(native_footprint_mean_mm=float(old['aso_watershed_mean_mm']),
                       native_footprint_valid_area_fraction=float(old['aso_valid_area_fraction']),
                       common_minus_native_mm=mean-float(old['aso_watershed_mean_mm']))
        if archived_regrid_directory:
            archive=Path(archived_regrid_directory)/f"{day:%Y%m%d}_swe_aso_on_noah.tif"
            with rasterio.open(archive) as src:
                if src.crs != crs or not src.transform.almost_equals(transform) or src.shape != weights.shape:
                    raise ValueError('Archived remap grid does not match verified target grid')
                old=src.read(1,masked=True).filled(np.nan).astype(float)
            paired=mask & np.isfinite(old) & np.isfinite(legacy_field)
            source_details['archived_unfiltered_regrid']={'file':archive.name,'sha256':sha256(archive),
                'paired_target_count':int(paired.sum()),
                'support_mismatch_count':int(np.sum(mask & (np.isfinite(old)!=np.isfinite(legacy_field)))),
                'maximum_absolute_difference_mm':float(np.max(abs(old[paired]-legacy_field[paired])))}
        rows.append(row);details.append(source_details)
        fields.append(value);coverages.append(coverage);legacy_coverages.append(legacy);raw_fields.append(unthresholded)
        print(f'{day}: {mean:.9f} mm; target-area coverage {fraction:.9%}; boundary changes {changed.sum()}',flush=True)
    output_directory.mkdir(parents=True,exist_ok=True)
    with (output_directory/output_names[0]).open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    np.savez_compressed(output_directory/output_names[1],date=np.array([r['date'] for r in rows]),
        aso_swe_mm=np.array(fields),valid_source_target_area_fraction=np.array(coverages),
        legacy_unpadded_coverage=np.array(legacy_coverages),unthresholded_swe_mm=np.array(raw_fields),
        weights=weights,fixed_domain_mask=mask,transform=np.array(tuple(transform)),crs=crs.to_wkt())
    manifest={'method':'GDAL area-overlap valid-source SWE mean then fractional watershed target-area basin mean',
        'minimum_valid_coverage':.5,'coverage_denominator':'Full target cell, including source-exterior space as missing',
        'basin_mean_equation':'sum(w_j*S_j over retained fixed-domain targets)/sum(w_j over same targets)',
        'w_j':'fractional watershed overlap * target planar cell area; no C_j multiplier',
        'reported_coverage':'sum(w_j retained)/sum(w_j fixed model domain); not native-pixel footprint',
        'boundary_semantics':'Target-cell SWE may use source area outside watershed; watershed weighting follows remapping',
        'coverage_proxy':'sum(w_j*C_j)/sum(w_j), reported separately; not exact native polygon coverage',
        'source_validity':'Read declared raster mask/NoData, apply scale and offset; retain all finite SWE including zeros/extremes',
        'area_convention':'Projected target-grid planar area, not geodesic area',
        'grid':{'file':grid_path.name,'sha256':sha256(grid_path),'crs_wkt':crs.to_wkt(),
                'transform':list(transform),'shape':list(weights.shape),'fixed_domain_target_count':int(mask.sum()),
                'fixed_watershed_area_m2':float(weights[mask].sum(dtype='float64')*cell_area)},
        'flights':{'file':flights_path.name,'sha256':sha256(flights_path)},
        'native_comparison':({'file':Path(native_comparison).name,'sha256':sha256(native_comparison)} if native_comparison else None),
        'code_sha256':sha256(Path(__file__)),
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'rasterio':rasterio.__version__,'gdal':rasterio.__gdal_version__},
        'outputs':[{'file':name,'sha256':sha256(output_directory/name)} for name in output_names[:2]],
        'sources':details}
    (output_directory/output_names[2]).write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return rows,manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--aso-directory',required=True,type=Path)
    parser.add_argument('--flights',required=True,type=Path)
    parser.add_argument('--grid-weights',required=True,type=Path)
    parser.add_argument('--output-directory',required=True,type=Path)
    parser.add_argument('--swe-units',choices=['auto','m','mm'],default='auto')
    parser.add_argument('--native-comparison',type=Path)
    parser.add_argument('--archived-regrid-directory',type=Path)
    args=parser.parse_args()
    extract(args.aso_directory,args.flights,args.grid_weights,args.output_directory,args.swe_units,
            args.native_comparison,args.archived_regrid_directory)

if __name__=='__main__':
    main()
