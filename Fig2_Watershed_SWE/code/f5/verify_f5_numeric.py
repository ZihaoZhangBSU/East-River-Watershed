"""Independent F5 reductions from saved fields using Python math.fsum.

Run: python verify_f5_numeric.py --data-directory data/f5_boundary
This does not call the generator's scores() implementation.
"""
import argparse,csv,json,math,hashlib
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-directory',required=True,type=Path)
    a=ap.parse_args();d=a.data_directory
    lookup={}
    with (d/'F5_elevation_metrics_corrected.csv').open() as stream:
        for row in csv.DictReader(stream):lookup[row['date'],row['variable'],row['model'],row['elevation_band']]=row
    with np.load(d/'F5_corrected_fields.npz',allow_pickle=False) as z:
        observed=z['ASO_mm'];modeled=z['model_mm'];dates=z['date'];variables=z['variable']
        names=z['model_names'];domain=z['fixed_domain_mask'];weights=z['weights'];dem=z['dem']
    max_error={key:0. for key in ['aso_mean_mm','mb_mm','rmse_mm','r','rmb_percent','rrmse_percent']}
    count=0;comparisons=[]
    for obs,models,date,var in zip(observed,modeled,dates,variables):
        for model,name in zip(models,names):
            for band in ['Lower','Middle','Upper','All']:
                valid=domain&(weights>0)&np.isfinite(obs)&np.isfinite(model)
                if band=='Lower':valid &= dem<2896.
                if band=='Middle':valid &= (dem>=2896.)&(dem<3353.)
                if band=='Upper':valid &= dem>=3353.
                pairs=list(zip(obs[valid],model[valid]));n=len(pairs)
                mean=math.fsum(float(o) for o,m in pairs)/n
                modelmean=math.fsum(float(m) for o,m in pairs)/n
                bias=math.fsum(float(m)-float(o) for o,m in pairs)/n
                rmse=math.sqrt(math.fsum((float(m)-float(o))**2 for o,m in pairs)/n)
                covariance=math.fsum((float(o)-mean)*(float(m)-modelmean) for o,m in pairs)
                sum_oo=math.fsum((float(o)-mean)**2 for o,m in pairs)
                sum_mm=math.fsum((float(m)-modelmean)**2 for o,m in pairs)
                correlation=covariance/math.sqrt(sum_oo*sum_mm) if sum_oo>0 and sum_mm>0 else None
                relative_bias=bias/mean*100 if mean else None
                relative_rmse=rmse/mean*100 if mean else None
                r=lookup[str(date),str(var),str(name),band]
                assert n==int(r['n_cells'])
                values={'aso_mean_mm':mean,'mb_mm':bias,'rmse_mm':rmse,'r':correlation,
                        'rmb_percent':relative_bias,'rrmse_percent':relative_rmse}
                for key,value in values.items():
                    if value is None:
                        assert r[key] in ['','None']
                        continue
                    error=abs(value-float(r[key]));max_error[key]=max(max_error[key],error)
                    assert error<(1e-12 if key=='r' else 1e-9),(date,var,name,band,key,error)
                count+=1
                comparisons.append({'date':str(date),'variable':str(var),'model':str(name),'elevation_band':band,
                    'n':n,'independent_reference_mean_mm':mean,'independent_MB_mm':bias,'independent_RMSE_mm':rmse,
                    'independent_r':correlation,'independent_RMB_percent':relative_bias,'independent_RRMSE_percent':relative_rmse})
    with (d/'F5_independent_numeric_checks.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
    audit=list(csv.DictReader((d/'F5_marker_placement_audit.csv').open()))
    # New main Fig. 5 adds per-flight spatial r. A replay may intentionally
    # select only SWE or depth, avoiding writes to other accepted figures.
    expected_metrics={'SWE':('mb_mm','rmse_mm','r'),'Depth':('mb_mm','rmse_mm','r')}
    variables_present={point['variable'] for point in audit}
    assert variables_present and variables_present<=set(expected_metrics)
    expected_keys={(str(date),var,metric) for date in dates
                   for var in variables_present for metric in expected_metrics[var]}
    actual_keys=[(p['date'],p['variable'],p['metric']) for p in audit]
    assert len(actual_keys)==len(set(actual_keys)) and set(actual_keys)==expected_keys
    plotted_markers=0
    for point in audit:
        for model,column in [('iSnobal','isnobal_x'),('Noah-MP','noahmp_x')]:
            raw=lookup[point['date'],point['variable'],model,'All'][point['metric']]
            if raw in ('','None'):
                assert point[column] in ('','None')
            else:
                assert float(point[column])==float(raw)
                plotted_markers+=1
    assert count==224
    result={'all_flight_model_variable_scores_checked':56,'all_elevation_rows_including_All_checked':count,
            'paired_cell_counts_checked':count,'plotted_markers_checked':plotted_markers,
            'maximum_absolute_difference':max_error,
            'difference_units':{'aso_mean_mm':'mm','mb_mm':'mm','rmse_mm':'mm','r':'dimensionless',
                                'rmb_percent':'percentage points','rrmse_percent':'percentage points'},'passed':True,
            'method':'Independent Python math.fsum reductions from saved full corrected fields for reference means, MB, RMSE, Pearson correlation, RMB, RRMSE and counts; model-specific finite-pair and independently selected elevation masks; all figure point records checked against recomputed All-band metric inputs.'}
    (d/'F5_independent_numeric_QA.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    # A rerun may leave earlier QA CSVs beside generator outputs. Refresh their
    # final hashes after verification so the manifest describes current files.
    manifest_path=d/'F5_correction_manifest.json'
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
        for item in manifest.get('outputs',[]):
            path=d/item['file']
            if path.is_file():item['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        manifest['independent_numeric_QA']={
            'file':'F5_independent_numeric_QA.json',
            'sha256':hashlib.sha256((d/'F5_independent_numeric_QA.json').read_bytes()).hexdigest(),
            'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'checked_rows':224,'checked_markers':plotted_markers,'passed':True}
        manifest_path.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
