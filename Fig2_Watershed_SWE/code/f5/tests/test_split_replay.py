"""Replay rejects incomplete/ambiguous scientific inputs and preserves artists."""
import copy
import csv
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

source=Path(__file__).resolve().parents[1]/'Fig5_ASO_Performance.py'
spec=importlib.util.spec_from_file_location('f5_split',source)
f5=importlib.util.module_from_spec(spec);spec.loader.exec_module(f5)


class ReplayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.directory=Path(self.temp.name)
        self.rows=[]
        for date in f5.FLIGHT_DATES:
            for variable in ('SWE','Depth'):
                for model,_ in f5.MODELS:
                    for band in f5.BANDS:
                        self.rows.append(dict(zip(f5.SCORE_COLUMNS,
                            (date,int(date[:4]),variable,model,band,
                             6 if band=='All' else 2,100.,10.,20.,.8,10.,20.))))

    def tearDown(self):
        f5.plt.close('all')
        self.temp.cleanup()

    def write(self,rows=None,columns=None):
        path=self.directory/'scores.csv'
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=columns or f5.SCORE_COLUMNS)
            writer.writeheader();writer.writerows(self.rows if rows is None else rows)
        return path

    def test_complete_validated_replay(self):
        rows=f5.load_scores_csv(self.write())
        self.assertEqual(len(rows),224)
        self.assertEqual(rows,self.rows)

    def test_missing_row_rejected(self):
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write(self.rows[:-1]))

    def test_duplicate_key_rejected(self):
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write(self.rows+[self.rows[0]]))

    def test_unknown_flight_rejected(self):
        self.rows[0]['date']='2020-03-31';self.rows[0]['water_year']=2020
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write())

    def test_nonfinite_scientific_value_rejected(self):
        self.rows[0]['mb_mm']='nan'
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write())

    def test_inconsistent_relative_metric_rejected(self):
        self.rows[0]['rmb_percent']=999
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write())

    def test_inconsistent_support_partition_rejected(self):
        self.rows[0]['n_cells']=3
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write())

    def test_negative_rmse_rejected(self):
        self.rows[0]['rmse_mm']=-20.
        with self.assertRaises(ValueError):f5.load_scores_csv(self.write())

    def test_model_specific_reference_cannot_be_silently_shared(self):
        rows=copy.deepcopy(self.rows)
        for row in rows:
            if row['model']=='Noah-MP':row['aso_mean_mm']=200.
        with patch.object(f5.plt.Figure,'savefig'):
            with self.assertRaisesRegex(ValueError,'common ASO reference support'):
                f5.plot(rows,self.directory,None)

    def test_one_combined_figure_has_all_168_markers_and_dates_in_both_rows(self):
        figures=[]
        def capture(fig,*args,**kwargs):
            if fig not in figures:figures.append(fig)
        with patch.object(f5.plt.Figure,'savefig',capture):
            audit=f5.plot(self.rows,self.directory,None)
        self.assertEqual(len(figures),1)
        self.assertEqual(len(audit),84)
        self.assertEqual([len(fig.axes) for fig in figures],[6])
        for fig in figures:
            self.assertEqual(tuple(fig.get_size_inches()),(6.5,7.2))
            self.assertEqual(sum(len(ax.collections) for ax in fig.axes),28*len(fig.axes))
            for index in (0,3):
                self.assertEqual([text.get_text() for text in fig.axes[index].get_yticklabels()],list(f5.FLIGHT_DATES))
            self.assertTrue(all(text.get_fontsize()>=8.05 for ax in fig.axes
                                for text in ax.get_xticklabels()+ax.get_yticklabels()))
        self.assertEqual({r['variable'] for r in audit[:42]},{'SWE'})
        self.assertEqual({r['variable'] for r in audit[42:]},{'Depth'})

    def test_main_only_does_not_save_depth_and_maps_spatial_r(self):
        saved=[];figures=[]
        def capture(fig,path,**kwargs):
            saved.append(Path(path).name)
            if fig not in figures:figures.append(fig)
        with patch.object(f5.plt.Figure,'savefig',capture):
            audit=f5.plot(self.rows,self.directory,None,'swe')
        self.assertEqual(saved,['JHM_F5_swe_audit.pdf','JHM_F5_swe_audit.png'])
        self.assertEqual(len(figures),1)
        self.assertEqual(len(audit),42)
        self.assertEqual([r['metric'] for r in audit[28:]],['r']*14)
        self.assertEqual([r['date'] for r in audit[28:]],list(f5.FLIGHT_DATES))
        self.assertTrue(all(r['isnobal_x']==.8 and r['noahmp_x']==.8 for r in audit[28:]))
        self.assertEqual(figures[0].axes[2].get_title(),'SWE spatial r')
        self.assertEqual(figures[0].axes[2].get_xlabel(),'')

    def test_zero_correlation_plots_but_missing_correlation_is_explicit(self):
        rows=copy.deepcopy(self.rows)
        for row in rows:
            if row['variable']=='SWE' and row['elevation_band']=='All' and row['date']==f5.FLIGHT_DATES[0]:
                row['r']=0. if row['model']=='iSnobal' else None
        figures=[]
        with patch.object(f5.plt.Figure,'savefig',lambda fig,*a,**k:figures.append(fig)):
            audit=f5.plot(rows,self.directory,None,'swe')
        ax=figures[0].axes[2]
        self.assertEqual(len(ax.collections),27)
        self.assertEqual(float(ax.collections[0].get_offsets()[0,0]),0.)
        self.assertTrue(any(text.get_text()=='N: NA' for text in ax.texts))
        self.assertEqual(audit[28]['isnobal_x'],0.)
        self.assertIsNone(audit[28]['noahmp_x'])
        self.assertEqual(f5.plotted_marker_count(audit),83)

    def test_depth_only_audit_has_three_metrics_and_unique_name(self):
        figures=[]
        saved=[]
        def capture(fig,path,**kwargs):
            figures.append(fig);saved.append(Path(path).name)
        with patch.object(f5.plt.Figure,'savefig',capture):
            audit=f5.plot(self.rows,self.directory,None,'depth')
        self.assertEqual(len(audit),42)
        self.assertEqual({r['variable'] for r in audit},{'Depth'})
        self.assertEqual(len(figures[0].axes),3)
        self.assertEqual(saved,['JHM_F5_depth_audit.pdf','JHM_F5_depth_audit.png'])
        self.assertEqual([r['metric'] for r in audit[28:]],['r']*14)
        self.assertEqual([r['panel'] for r in audit[28:]],[6]*14)

    def test_invalid_figure_selection_rejected(self):
        with self.assertRaises(ValueError):f5.plot(self.rows,self.directory,None,'unknown')

    def test_main_only_manifest_excludes_preexisting_depth_sentinel(self):
        for stem in ('JHM_F5_aso_performance','JHM_F5_swe_audit','JHM_S_ASO_depth_performance'):
            for extension in ('pdf','png'):
                (self.directory/f'{stem}.{extension}').write_bytes(b'preexisting sentinel')
        selected=f5.generated_figures(self.directory,'swe')
        self.assertEqual([p.name for p in selected],
            ['JHM_F5_swe_audit.pdf','JHM_F5_swe_audit.png'])
        self.assertEqual((self.directory/'JHM_S_ASO_depth_performance.pdf').read_bytes(),b'preexisting sentinel')

    def test_combined_output_never_writes_or_claims_accepted_depth_asset(self):
        sentinel=self.directory/'JHM_S_ASO_depth_performance.pdf'
        sentinel.write_bytes(b'accepted depth sentinel')
        saved=[]
        def capture(fig,path,**kwargs):
            Path(path).write_bytes(b'combined output');saved.append(Path(path).name)
        with patch.object(f5.plt.Figure,'savefig',capture):
            f5.plot(self.rows,self.directory,None)
        self.assertEqual(saved,['JHM_F5_aso_performance.pdf','JHM_F5_aso_performance.png'])
        self.assertEqual([p.name for p in f5.generated_figures(self.directory,'both')],saved)
        self.assertEqual(sentinel.read_bytes(),b'accepted depth sentinel')

    def test_current_shared_axes_and_safe_future_expansion(self):
        figures=[]
        with patch.object(f5.plt.Figure,'savefig',lambda fig,*a,**k:figures.append(fig)):
            f5.plot(self.rows,self.directory,None)
        for top,bottom,expected in zip(figures[0].axes[:3],figures[0].axes[3:],
                ((-300.,200.),(0.,500.),(.65,1.))):
            self.assertEqual(top.get_xlim(),expected)
            self.assertEqual(bottom.get_xlim(),expected)
        rows=copy.deepcopy(self.rows)
        for row in rows:
            if row['variable']=='Depth' and row['elevation_band']=='All':row['r']=-.8
        limits=f5.performance_limits(rows,{'SWE','Depth'})
        self.assertLessEqual(limits['r'][0],-.8)
        self.assertEqual(limits['r'][1],1.)


if __name__=='__main__':unittest.main()
