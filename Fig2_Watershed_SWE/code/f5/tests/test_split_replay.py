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

    def test_two_figures_contain_all_112_markers_and_14_dates(self):
        figures=[]
        def capture(fig,*args,**kwargs):
            if fig not in figures:figures.append(fig)
        with patch.object(f5.plt.Figure,'savefig',capture):
            audit=f5.plot(self.rows,self.directory,None)
        self.assertEqual(len(figures),2)
        self.assertEqual(len(audit),56)
        self.assertEqual([len(fig.axes) for fig in figures],[2,2])
        for fig in figures:
            self.assertEqual(tuple(fig.get_size_inches()),(6.5,5.3))
            self.assertEqual(sum(len(ax.collections) for ax in fig.axes),56)
            self.assertEqual([text.get_text() for text in fig.axes[0].get_yticklabels()],list(f5.FLIGHT_DATES))
            self.assertTrue(all(text.get_fontsize()>=8 for ax in fig.axes
                                for text in ax.get_xticklabels()+ax.get_yticklabels()))
        self.assertEqual({r['variable'] for r in audit[:28]},{'SWE'})
        self.assertEqual({r['variable'] for r in audit[28:]},{'Depth'})


if __name__=='__main__':unittest.main()
