"""VW TDS utilities script

File: utilities.py
Description: Collection of functions for VW TDS program
Created by: Kyle O'Connor (keoconnor@usgs.gov)
Created: 2020-07-17
Last Modified: 2021-03-26
Developed for: Python 3.6+
"""

# Standard library imports
import csv
import datetime
import os
import statistics

# Third party imports
import numpy as np
import pandas as pd
import requests
from scipy import stats

# remove warnings for overwriting pandas dataframes for future pandas versions
pd.options.mode.chained_assignment = None
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class VolumeWeightedTDS(requests.Session):

    """Class for getting Volume Weighted Total Dissolved Solids data
    from a lake or reservoir. Four constructors need to be defined and are
    shown in the __init__ method."""

    def __init__(self, input_sites: str, stage_storage_relation: str,
                 output_location: str, comp_files_location: str):
        """Constructor for VolumeWeightedTDS class

        :param input_sites: list of site name and site id in a tuple
        :type input_sites: list
        :param stage_storage_relation: csv file of stage and storage relation
        :type stage_storage_relation: str
        :param output_location: directory for output files
        :type output_location: str
        :param comp_files_location: directory for comp files
        :type comp_files_location: str
        """

        super(VolumeWeightedTDS, self).__init__()
        self.input_sites = input_sites
        self.stage_storage_relation = stage_storage_relation
        self.site_names = [site[0] for site in self.input_sites]
        self.site_no = [site[1] for site in self.input_sites]
        # this may need updated after current NWIS deprecation
        self.html_locations = [(
            site[0], f'https://nwis.waterdata.usgs.gov/usa/nwis/qwdata/'
            f'?site_no={site[1]}&agency_cd=USGS&inventory_output=0&rdb_'
            f'inventory_output=file&TZoutput=0&pm_cd_compare=Greater%'
            f'20than&radio_parm_cds=all_parm_cds&qw_attributes=0&format='
            f'rdb&qw_sample_wide=wide&rdb_qw_attributes=0&date_format='
            f'YYYY-MM-DD&rdb_compression=value&submitted_form=brief_list'
        ) for site in self.input_sites]
        self.output_location = output_location
        self.comp_files_location = comp_files_location
        self.stage_to_storage_data = {}
        self.regression_data = {}
        self.ordered_site_list = ()
        self.vw_tds_data = {}
        self.salt_mass = {}

    def csv_to_dict(self):
        """Converts csv to a dict"""

        stage_to_stor_dict = {}
        stor_to_stage_dict = {}
        with open(self.stage_storage_relation, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                stage = row['Stage']
                stor = row['Volume']
                stage_to_stor_dict[stage] = stor
                stor_to_stage_dict[stor] = stage

        self.stage_to_storage_data = stage_to_stor_dict

    def get_raw_data(self):
        """Gets raw data from NWIS and creates a csv file of data"""

        for site, url in self.html_locations:
            request = self.get(url, timeout=120)
            if request.status_code != 200:
                raise Exception(f'ERROR: Unable to get site data:'
                                f'HTTP ERROR: {request.status_code}')
            text_from_url = request.text.splitlines()
            site_file = site + '.csv'
            site_file_location = os.path.join(
                self.comp_files_location, site_file
            )
            with open(site_file_location, 'w') as output:
                for line in text_from_url:
                    if '\t' not in line:
                        continue
                    if '#' not in line:
                        output.write(f'{line}\n')

    def get_regression(self):
        """Gets SC and TDS regression values and creates a regression graph"""

        sc_tds = []
        regression_dict = {}
        for site in self.site_names:
            site_file = site + '.csv'
            site_file_location = os.path.join(
                self.comp_files_location, site_file
            )
            with open(site_file_location, 'r') as file:
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    sc = row.get('p00095', None)
                    if sc is not None:
                        sc = sc.strip()
                    if not sc:
                        continue
                    if not sc.isnumeric():
                        continue
                    tds = row.get('p70300', None)
                    if tds is not None:
                        tds = tds.strip()
                    if not tds:
                        continue
                    if not tds.isnumeric():
                        continue
                    sc_tds.append((float(sc), float(tds)))


        ratios = []
        sc_tds_filtered = []
        for sc, tds in sc_tds:
            ratio = tds/sc
            ratios.append(ratio)
        stdev = statistics.stdev(ratios)
        for sc, tds in sc_tds:
            if tds/sc > statistics.mean(ratios) + (2*stdev):
                continue
            elif tds/sc < statistics.mean(ratios) - (2*stdev):
                continue
            else:
                sc_tds_filtered.append((sc, tds))

        sc = [va[0] for va in sc_tds_filtered]
        tds = [va[1] for va in sc_tds_filtered]
        reg_data = stats.linregress(sc, tds)
        print(reg_data)
        reg_x = np.arange(np.min(sc), np.max(sc), 100)
        reg_y = (reg_data[0] * reg_x) + reg_data[1]
        p_va = reg_data[3]
        if p_va < 0.001:
            print(f'Regression is statisticall significant. '
                  f'p-value = {p_va}')
        else:
            print(f'Regression is not statistically significant. '
                  f'p-value = {p_va}')

        regression_dict['slope'] = reg_data[0]
        regression_dict['intercept'] = reg_data[1]
        self.regression_data = regression_dict
        return {
            'regression_data': regression_dict, 
            'sc': sc, 
            'tds': tds
        }

    def get_deepest_data(self):
        """Determines the deepest site"""

        number_of_sites = len(self.input_sites)
        site_file = self.site_names[0] + '.csv'
        site_file_location = os.path.join(self.comp_files_location, site_file)
        df = pd.read_csv(site_file_location, delimiter='\t')
        df.dropna(how='all')
        df.drop(df.index[0], inplace=True)
        df['sample_dt'] = pd.to_datetime(df['sample_dt'])
        df['p00003'] = df['p00003'].astype(float)  # Sampling depth
        datetime_list = pd.to_datetime(df['sample_dt'].unique().tolist())
        set_date = datetime_list[0]
        depth_data = [[] for _ in range(number_of_sites)]
        for site in range(len(self.site_names)):
            site_file = self.site_names[site] + '.csv'
            site_file_location = os.path.join(self.comp_files_location,
                                              site_file)
            df = pd.read_csv(
                site_file_location, delimiter='\t'
            )
            df.dropna(how='all')
            df.drop(df.index[0], inplace=True)
            df.fillna(-1, inplace=True)
            df['p70300'] = df['p70300'].astype(float)
            df = df[df['p70300'] < 0]
            df['sample_dt'] = pd.to_datetime(df['sample_dt'])
            df['p00003'] = df['p00003'].astype(float)   # Sampling depth
            df = df.loc[df['sample_dt'] == set_date]
            depth = df['p00003'].tolist()
            depth_data[site].append(depth)
        depth_dict = {
            max(depth_data[0][0]): self.site_names[0], max(depth_data[1][0]):
                self.site_names[1], max(depth_data[2][0]): self.site_names[2]
        }
        depth_dict_reversed = {
            self.site_names[0]: depth_data[0][0], self.site_names[1]:
                depth_data[1][0], self.site_names[2]: depth_data[2][0]
        }
        depths_from_dict = []
        sites_from_dict = []
        for key, value in depth_dict.items():
            depths_from_dict.append(key)
            sites_from_dict.append(self.site_names)
        sites_list = sorted(
            depth_dict_reversed, key=depth_dict_reversed.get
        )
        x=1

        self.ordered_site_list = tuple(sites_list)

    def get_volume_weighted_tds(self):
        """Calculates the volume-weighted TDS for a lake"""

        mgl2kgacft = 1.23348
        bottom = float(
            list(self.stage_to_storage_data.keys())[0]
        )
        vw_tds_dict = {}
        stage_list = []
        storage_list = []
        slope = self.regression_data['slope']
        intercept = self.regression_data['intercept']
        stage2stor = self.stage_to_storage_data
        # use dates from deepest site
        site_file = self.ordered_site_list[0] + '.csv'
        site_file_location = os.path.join(self.comp_files_location, site_file)
        df = pd.read_csv(site_file_location, delimiter='\t')
        df.dropna(how='all')
        df.drop(df.index[0], inplace=True)
        sample_dt = list(set(df['sample_dt']))
        dates = sorted(
            [datetime.datetime.strptime(dt, '%Y-%m-%d') for dt in sample_dt]
        )
        for date in dates:
            test_dict = df[df['p00062'].astype(float) > 0]
            test_date = date.strftime('%Y-%m-%d')
            test_dict = test_dict[test_dict['sample_dt'] == test_date]
            if test_dict.empty:
                if date == dates[-1]:
                    stage_last = input(f'Enter stage datum for {date}: ')
                    storage_last = input(f'Enter storage datum for {date}: ')
                    last_date = date.strftime('%Y-%m-%d')
                else:
                    stage_last = None
                    storage_last = None
                    last_date = None

        for i, site in enumerate(self.ordered_site_list):
            if i == 0:
                site_file = site + '.csv'
                site_file_location = os.path.join(
                    self.comp_files_location, site_file
                )
                df = pd.read_csv(site_file_location, delimiter='\t')
                df.dropna(how='all')
                df.drop(df.index[0], inplace=True)
                df['Datetime'] = [
                    datetime.datetime.strptime(
                        dt, '%Y-%m-%d'
                    ) for dt in list(df['sample_dt'])
                ]
                if stage_last:
                    df.loc[
                        df['sample_dt'] == last_date, 'p00062'
                    ] = stage_last
                    df.loc[
                        df['sample_dt'] == last_date, 'p00054'
                    ] = storage_last
                df_stage_storage = df[df['p00062'].astype(float) > 0]
                df.fillna(-1, inplace=True)
                vw_tds_dict[site] = {}
                for date in dates:
                    stage_storage = df_stage_storage.loc[
                        df_stage_storage['Datetime'] == date
                    ]
                    stage = set(
                        round(stage_storage['p00062'].astype(float), 1)
                    )
                    storage = set(
                        round(stage_storage['p00054'].astype(float), 1)
                    )
                    if bool(stage) is False and bool(storage) is False:
                        continue
                    stage = list(stage)[0]
                    stage_list.append(stage)
                    storage = list(storage)[0]
                    storage_list.append(storage)
                    df2 = df.loc[df['Datetime'] == date]
                    df3 = df2[df2['p70300'].astype(float) < 0]
                    df3['Depth'] = df3['p00003'].astype(float)
                    df3 = df3.sort_values('Depth')
                    df3.dropna(subset=['p00095'], how='all', inplace=True)
                    depths = list(round(df3['Depth'], 1))
                    stage_m_depth = []
                    for depth in depths:
                        if stage - depth > bottom:
                            stage_m_depth.append(stage - depth)
                        else:
                            stage_m_depth.append(bottom)
                    stage_m_depth = [
                        round(depth, 1) for depth in stage_m_depth
                    ]
                    storages = [
                        stage2stor[
                            str(stage_m_depth_value)
                        ] for stage_m_depth_value in stage_m_depth
                    ]
                    storages = [float(i) for i in storages]
                    storages.insert(0, storage)
                    storages.append(0)
                    sliver = [l - m for l, m in zip(storages, storages[1:])]
                    sc_i = df3['p00095'].astype(float)
                    tds_mgl = slope * sc_i + intercept
                    tds_kgacft = tds_mgl * mgl2kgacft
                    vwtds_sliver = [
                        s * t for s, t in zip(sliver, tds_kgacft)
                    ]
                    vwtds_i = sum(vwtds_sliver)/storage
                    vwtds_mgl = vwtds_i/mgl2kgacft
                    vw_tds_dict[site][date] = vwtds_mgl
            else:
                site_file = site + '.csv'
                site_file_location = os.path.join(
                    self.comp_files_location, site_file
                )
                df = pd.read_csv(site_file_location, delimiter='\t')
                df.dropna(how='all')
                df.drop(df.index[0], inplace=True)
                df['Datetime'] = [
                    datetime.datetime.strptime(
                        dt, '%Y-%m-%d'
                    ) for dt in list(df['sample_dt'])
                ]
                if stage_last:
                    df.loc[
                        df['sample_dt'] == last_date, 'p00062'
                    ] = stage_last
                    df.loc[
                        df['sample_dt'] == last_date, 'p00054'
                    ] = storage_last
                df_stage_storage = df[df['p00062'].astype(float) > 0]
                df.fillna(-1, inplace=True)
                vw_tds_dict[site] = {}

                deepest_site_file = self.ordered_site_list[0] + '.csv'
                deepest_site_file_location = os.path.join(
                    self.comp_files_location, deepest_site_file
                )
                df_deepest = pd.read_csv(
                    deepest_site_file_location, delimiter='\t'
                )
                df_deepest.dropna(how='all')
                df_deepest.drop(df_deepest.index[0], inplace=True)
                df_deepest['Datetime'] = [
                    datetime.datetime.strptime(
                        dt, '%Y-%m-%d'
                    ) for dt in list(df_deepest['sample_dt'])
                ]
                df_deepest.fillna(-1, inplace=True)
                for date in dates:
                    stage_storage = df_stage_storage.loc[
                        df_stage_storage['Datetime'] == date
                    ]
                    stage = set(
                        round(stage_storage['p00062'].astype(float), 1)
                    )
                    storage = set(
                        round(stage_storage['p00054'].astype(float), 1)
                    )
                    if bool(stage) is False and bool(storage) is False:
                        continue
                    stage = list(stage)[0]
                    storage = list(storage)[0]
                    df2 = df.loc[df['Datetime'] == date]
                    df3 = df2[df2['p70300'].astype(float) < 0]
                    df3['Depth'] = df3['p00003'].astype(float)
                    df3 = df3.sort_values('Depth')
                    df3.dropna(subset=['p00095'], how='all', inplace=True)
                    df_deep = df_deepest.loc[
                        df_deepest['Datetime'] == date
                    ]
                    df3_depths = list(df3['Depth'].astype(float))
                    bottom_depth = df3_depths[-1]
                    df_deep = df_deep[df_deep['p70300'].astype(float) < 0]
                    df_deep['Depth'] = df_deep['p00003'].astype(float)
                    df_deep = df_deep.loc[df_deep['Depth'] > bottom_depth]
                    df_deep = df_deep.sort_values('Depth')
                    df4 = df3.append(df_deep, sort=False)
                    depths = list(round(df4['Depth'], 1))
                    stage_m_depth = []
                    for depth in depths:
                        if stage - depth > bottom:
                            stage_m_depth.append(stage - depth)
                        else:
                            stage_m_depth.append(bottom)
                    stage_m_depth = [
                        round(depth, 1) for depth in stage_m_depth
                    ]
                    storages = [
                        stage2stor[
                            str(stage_m_depth_value)
                        ] for stage_m_depth_value in stage_m_depth
                    ]
                    storages = [float(i) for i in storages]
                    storages.insert(0, storage)
                    storages.append(0)
                    storages = [
                        storage for storage in storages
                        if str(storage) != 'nan'
                    ]
                    sliver = [l - m for l, m in zip(storages, storages[1:])]
                    sc_i = df4['p00095'].astype(float)
                    tds_mgL = slope * sc_i + intercept
                    tds_kgacft = tds_mgL * mgl2kgacft
                    vwtds_sliver = [
                        s * t for s, t in zip(sliver, tds_kgacft)
                    ]
                    vwtds_i = sum(vwtds_sliver)/storage
                    vwtds_mgl = vwtds_i/mgl2kgacft
                    if vwtds_mgl < 0:
                        vwtds_mgl = statistics.mean(
                            [sc for sc in list(sc_i) if sc > 0]
                        ) * slope + intercept
                    vw_tds_dict[site][date] = vwtds_mgl
        final_df = pd.DataFrame(vw_tds_dict)
        final_df['Mean'] = final_df.mean(axis=1)
        final_df['Stage (ft)'] = stage_list
        final_df['Storage (ac-ft)'] = storage_list
        final_df.to_csv(self.output_location + '\\VW-TDS_Results.csv')

        self.vw_tds_data = final_df.to_dict()
        return final_df.to_dict()

    def get_salt_mass(self):
        """Calculates Salt Mass and exports a csv to output folder"""

        vw_tds = pd.DataFrame(self.vw_tds_data)
        mean_vw_tds = np.array(vw_tds['Mean'])
        mean_vw_tds_kg_acft = mean_vw_tds * 1.23348
        storage = np.array(vw_tds['Storage (ac-ft)'])
        salt_mass = mean_vw_tds_kg_acft * storage
        salt_mass_mt = salt_mass/907184740
        salt_mass_df = pd.DataFrame({
            'Date': vw_tds.index.values, 'Salt Mass (kg)': salt_mass,
            'Salt Mass (Millions of tons)': salt_mass_mt

        })
        salt_mass_df.to_csv(self.output_location + '\\Salt_Mass.csv')

        self.salt_mass = salt_mass_df.to_dict()


