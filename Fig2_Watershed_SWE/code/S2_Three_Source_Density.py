"""Correct S2 density support using independently corrected ASO target fields.

Use --fields for F5_corrected_fields.npz, made by the unfiltered full-target
ASO remapping extractor. SWE/depth are separately remapped before taking their
ratio. Density distributions use equal cell weights, not Fig.2 area weights.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


SOURCES=("ASO","iSnobal","Noah-MP")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def support(fields,threshold=100.):
    required=("ASO_SWE","ASO_depth","iSnobal_SWE","iSnobal_depth","Noah-MP_SWE","Noah-MP_depth")
    mask=fields["basin"].copy()
    for key in required:
        mask &= np.isfinite(fields[key])
    mask &= (fields["ASO_depth"]>threshold)&(fields["iSnobal_depth"]>0)&(fields["Noah-MP_depth"]>0)
    return mask


def independent_quantile(values,p):
    ordered=sorted(float(v) for v in values)
    pos=(len(ordered)-1)*p
    lo=int(math.floor(pos)); hi=int(math.ceil(pos))
    return ordered[lo]*(hi-pos)+ordered[hi]*(pos-lo) if lo!=hi else ordered[lo]


def summarize(fields_path,output_dir,old_fields=None,old_summary=None):
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    with np.load(fields_path,allow_pickle=False) as z:
        data={k:z[k] for k in z.files}
    dates=np.unique(data["date"])
    if len(dates)!=14 or list(data["model_names"])!=["iSnobal","Noah-MP"]:
        raise ValueError("Expected all14 flights and both established models")
    weights=data["weights"].astype(float)
    basin=(weights>0)&np.isfinite(data["dem"])
    if weights.shape!=(66,54) or int(basin.sum())!=842:
        raise ValueError("Wrong analysis grid or watershed support")
    prior=np.load(old_fields,allow_pickle=False) if old_fields else None
    summaries,cells,sensitivity,metrics,impact=[],[],[],[],[]
    all_fields={"weights":weights,"dem":data["dem"]}
    checks=0
    for date in dates:
        a={"basin":basin}
        for var,slug in (("SWE","SWE"),("Depth","depth")):
            inds=np.flatnonzero((data["date"]==date)&(data["variable"]==var))
            if len(inds)!=1:
                raise ValueError(f"Duplicate or missing field {date}/{var}")
            idx=int(inds[0]); raw=data["ASO_mm"][idx]
            # Enforce full-target coverage independently for both components.
            aso=np.where(data["full_target_valid_coverage"][idx]>=.5,raw,np.nan)
            a["ASO_"+slug]=aso
            a["iSnobal_"+slug]=data["model_mm"][idx,0]
            a["Noah-MP_"+slug]=data["model_mm"][idx,1]
        common=support(a)
        fraction=float(weights[common].sum()/weights[basin].sum())
        compact=str(date).replace("-","")
        all_fields[compact+"_common"]=common
        for source in SOURCES:
            all_fields[compact+"_"+source.replace("-","")+"_SWE"]=a[source+"_SWE"]
            all_fields[compact+"_"+source.replace("-","")+"_depth"]=a[source+"_depth"]
        densities={s:1000*a[s+"_SWE"][common]/a[s+"_depth"][common] for s in SOURCES}
        for source,values in densities.items():
            q=np.percentile(values,[5,25,50,75,95])
            mean=float(values.mean())
            if abs(mean-math.fsum(values)/len(values))>1e-10:
                raise AssertionError("Scalar mean check failed")
            checks+=1
            for p,value in zip((.05,.25,.5,.75,.95),q):
                if abs(value-independent_quantile(values,p))>1e-10:
                    raise AssertionError("Independent sorted quantile check failed")
                checks+=1
            summaries.append({"date":str(date),"water_year":int(str(date)[:4]),"source":source,
                "n_cells":int(common.sum()),"watershed_area_fraction":fraction,"mean_kg_m3":mean,
                "p05":q[0],"p25":q[1],"median":q[2],"p75":q[3],"p95":q[4],
                "minimum":float(values.min()),"maximum":float(values.max()),
                "zero_swe_count":int((a[source+"_SWE"][common]==0).sum()),
                "negative_swe_count":int((a[source+"_SWE"][common]<0).sum())})
            for (row,col),density in zip(np.argwhere(common),values):
                cells.append({"date":str(date),"source":source,"row":int(row),"col":int(col),
                    "density_kg_m3":float(density),"swe_mm":float(a[source+"_SWE"][row,col]),
                    "depth_mm":float(a[source+"_depth"][row,col]),"watershed_fraction":float(weights[row,col])})
        for threshold in (0.,100.,200.):
            mask=support(a,threshold)
            for source in SOURCES:
                vals=1000*a[source+"_SWE"][mask]/a[source+"_depth"][mask]
                sensitivity.append({"date":str(date),"source":source,"aso_depth_gt_mm":threshold,
                    "n_cells":int(mask.sum()),"watershed_area_fraction":float(weights[mask].sum()/weights[basin].sum()),
                    "mean_kg_m3":float(vals.mean()),"median_kg_m3":float(np.median(vals))})
        obs=densities["ASO"]
        for source in ("iSnobal","Noah-MP"):
            pred=densities[source];diff=pred-obs
            metrics.append({"date":str(date),"model":source,"n_cells":len(obs),
                "aso_mean_kg_m3":float(obs.mean()),"model_mean_kg_m3":float(pred.mean()),
                "MB_kg_m3":float(diff.mean()),"RMSE_kg_m3":float(np.sqrt(np.mean(diff**2))),
                "r":float(np.corrcoef(obs,pred)[0,1]),
                "relative_MB_percent":float(100*diff.mean()/obs.mean()),
                "relative_RMSE_percent":float(100*np.sqrt(np.mean(diff**2))/obs.mean())})
        if prior is not None:
            oldmask=prior[compact+"_common"].astype(bool)
            removed=np.argwhere(oldmask&~common)
            added=np.argwhere(common&~oldmask)
            rec={"date":str(date),"old_n_cells":int(oldmask.sum()),"new_n_cells":int(common.sum()),
                "removed_cells":removed.tolist(),"added_cells":added.tolist(),"mean_changes_kg_m3":{}}
            for source in SOURCES:
                key=source.replace("-","")
                previous=1000*prior[compact+"_"+key+"_SWE"][oldmask]/prior[compact+"_"+key+"_depth"][oldmask]
                rec["mean_changes_kg_m3"][source]=float(densities[source].mean()-previous.mean())
            impact.append(rec)
    summary=pd.DataFrame(summaries)
    summary.to_csv(output_dir/"S2_bulk_density_summary.csv",index=False,float_format="%.17g")
    pd.DataFrame(cells).to_csv(output_dir/"S2_bulk_density_common_cells.csv",index=False,float_format="%.17g")
    pd.DataFrame(sensitivity).to_csv(output_dir/"S2_depth_screen_sensitivity.csv",index=False,float_format="%.17g")
    pd.DataFrame(metrics).to_csv(output_dir/"S2_density_model_metrics.csv",index=False,float_format="%.17g")
    summary[summary.source.eq("ASO")][["date","water_year","n_cells","watershed_area_fraction"]].to_csv(output_dir/"S2_flight_support.csv",index=False,float_format="%.17g")
    np.savez_compressed(output_dir/"S2_bulk_density_common_fields.npz",**all_fields)
    if old_summary:
        old=pd.read_csv(old_summary)
        joined=summary.merge(old,on=["date","source"],suffixes=("_new","_old"))
        for name in ("n_cells","watershed_area_fraction","mean_kg_m3","p05","p25","median","p75","p95"):
            joined[name+"_change"]=joined[name+"_new"]-joined[name+"_old"]
        joined.to_csv(output_dir/"S2_old_vs_new_distribution_statistics.csv",index=False,float_format="%.17g")
    provenance={"method":"S2 three-source density after separate SWE/depth remapping; full-target ASO validity>=0.5 separately for both fields; all components finite; ASO depth>100mm, model depths>0; equal-cell distribution statistics, no density clipping",
        "n_flights":14,"n_summaries":len(summary),"n_source_cell_samples":len(cells),
        "independent_mean_and_sorted_quantile_checks":checks,
        "support_changes":impact,
        "field_input_sha256":sha256(fields_path),"code_sha256":sha256(__file__),
        "old_fields_sha256":sha256(old_fields) if old_fields else None,
        "old_summary_sha256":sha256(old_summary) if old_summary else None,
        "source_operator_provenance":"../f5_boundary/F5_correction_manifest.json",
        "finite_zeros":"Retained for SWE; ASO depth screen and positive model depth determine ratio admissibility.",
        "limits":"ASO-implied density is not independent density observation. Shared positive-depth selection excludes model-missed snow. Flight supports can change; distributions are spatial variability, not uncertainty intervals."}
    (output_dir/"S2_density_correction_manifest.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    return summary


def render(summary,figure_dir,font=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.lines import Line2D
    figure_dir=Path(figure_dir);figure_dir.mkdir(parents=True,exist_ok=True)
    family="DejaVu Sans"
    if font:
        font_manager.fontManager.addfont(str(font))
        family=font_manager.FontProperties(fname=font).get_name()
    plt.rcParams.update({"font.family":family,"font.size":8,"axes.labelsize":8,
        "xtick.labelsize":7.5,"ytick.labelsize":7.5,"legend.fontsize":7.5,
        "pdf.fonttype":42,"ps.fonttype":42,"axes.linewidth":.65})
    dates=list(summary.date.drop_duplicates());positions=[];cur=0.;prev=None
    for date in dates:
        year=int(date[:4])
        if prev is not None and year!=prev:cur+=.65
        positions.append(cur);cur+=1.;prev=year
    colors={"iSnobal":"#0072B2","Noah-MP":"#D55E00","ASO":"black"}
    markers={"iSnobal":"o","Noah-MP":"s","ASO":"D"}
    fig,ax=plt.subplots(figsize=(6.5,3.5))
    fig.subplots_adjust(left=.105,right=.985,top=.81,bottom=.22)
    for source,offset in (("iSnobal",-.23),("Noah-MP",0),("ASO",.23)):
        v=summary[summary.source.eq(source)].set_index("date").loc[dates]
        xp=np.asarray(positions)+offset
        ax.vlines(xp,v.p05,v.p95,color=colors[source],lw=.7,zorder=2)
        ax.vlines(xp,v.p25,v.p75,color=colors[source],lw=2.8,zorder=3)
        ax.scatter(xp,v["median"],color=colors[source],edgecolor="white",lw=.3,s=10,marker=markers[source],zorder=4)
    ax.set_xticks(positions,[pd.Timestamp(date).strftime("%d %b") for date in dates],rotation=55,ha="right")
    ax.set_xlim(min(positions)-.65,max(positions)+.65)
    ax.set_ylim(0,np.ceil(summary.p95.max()/50)*50+25)
    ax.set_ylabel("Common-grid bulk density (kg m$^{-3}$)")
    ax.spines[["top","right"]].set_visible(False)
    ax.yaxis.grid(True,color=".85",ls="--",lw=.45);ax.set_axisbelow(True)
    trans=ax.get_xaxis_transform()
    for year in sorted({date[:4] for date in dates}):
        inds=[i for i,date in enumerate(dates) if date.startswith(year)]
        ax.text(np.mean(np.asarray(positions)[inds]),1.035,"WY"+year,ha="center",va="bottom",fontsize=7.5,transform=trans)
        ax.plot([positions[inds[0]]-.42,positions[inds[-1]]+.42],[1.015,1.015],color=".45",lw=.5,transform=trans,clip_on=False)
    handles=[Line2D([],[],color=colors[s],marker=markers[s],ls="none",ms=3.5,label=s) for s in ("iSnobal","Noah-MP","ASO")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.54,.995),ncol=3,frameon=False,handletextpad=.4,columnspacing=1.4)
    for ext in ("pdf","png"):
        fig.savefig(figure_dir/("JHM_S2_ASO_three_source_bulk_density."+ext),dpi=600,facecolor="white")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fields",required=True,type=Path)
    p.add_argument("--output-dir",required=True,type=Path)
    p.add_argument("--figure-dir",type=Path)
    p.add_argument("--old-fields",type=Path)
    p.add_argument("--old-summary",type=Path)
    p.add_argument("--font",type=Path)
    a=p.parse_args()
    d=summarize(a.fields,a.output_dir,a.old_fields,a.old_summary)
    if a.figure_dir:
        render(d,a.figure_dir,a.font)
    print(f"Validated{len(d)} summaries and regenerated S2 density support")


if __name__=="__main__":
    main()
