import datetime
import inspect
import json
import os
import sys

# assert correct versions
PY_VERSION_INFO = sys.version_info
MAJOR, MINOR = (PY_VERSION_INFO[0], PY_VERSION_INFO[1])
assert MAJOR >= 3 and MINOR >= 6

# Local imports
SCRIPT_PATH = os.path.realpath(inspect.stack()[0][1])
SCRIPT_LOC = os.path.dirname(SCRIPT_PATH)
sys.path.insert(0, os.path.dirname(SCRIPT_LOC))
sys.path.insert(0, SCRIPT_LOC)
from utilities import VolumeWeightedTDS

# remove pandas warnings - not updating pandas
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def lambda_handler(event, context):

    start_time = datetime.datetime.now()
        # Get site data
    sites = []
    sites_location = os.path.join(SCRIPT_LOC, 'sites.txt')
    with open(sites_location, 'r') as f:
        for line in f.readlines():
            site_name, site_no = line.strip().split('\t')
            site_id = (site_name, site_no)
            sites.append(site_id)

    # Get comp files location
    comp_loc = 'comp_files'
    comp_files_location = os.path.join(SCRIPT_LOC, comp_loc)
    if not os.path.exists(comp_files_location):
        os.mkdir(comp_files_location)

    # Get Stage-Storage Relation
    stage_storage_file = 'LakeStageVolRelation.csv'
    stage_storage_relation = os.path.join(
        SCRIPT_LOC, stage_storage_file
    )

    # Create output directory
    output_dir = os.path.join(SCRIPT_LOC, 'output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output = f'output_{start_time.strftime("%Y%M%d_%H%M%S")}'
    output_location = os.path.join(output_dir, output)
    os.mkdir(output_location)

    # initiate VolumeWeightedTDS object
    print('Initializing Volume Weighted TDS Calc')
    VolumeWeightedObject = VolumeWeightedTDS(
        sites, stage_storage_relation, output_location, comp_files_location
    )
    # get data from NWIS
    print('Getting data from NWIS...')
    try:
        VolumeWeightedObject.get_raw_data()
        print('Data obtained')
    except Exception as e:
        print(f'Error in getting data from NWIS -- Error: {e}')
        sys.exit()

    print('Getting regression data...')
    # try:
    regression_data = VolumeWeightedObject.get_regression()
    print('regression obtained')
    # except Exception as e:
    #     print(f'Error in getting regression -- Error: {e}')
    #     sys.exit()

    print('Getting order list of deepest to shallowest sites...')
    try:
        VolumeWeightedObject.get_deepest_data()
        print('list obtained')
    except Exception as e:
        print(f'Error in getting ordered list -- Error: {e}')
        sys.exit()

    print('Getting stage-storage relation...')
    try:
        VolumeWeightedObject.csv_to_dict()
        print('Stage-storage relation obtained')
    except Exception as e:
        print(f'Error in getting stage-storage relation -- Error: {e}')
        sys.exit()

    print('Getting volume-weighted TDS...')
    try:
        vwtds_data = VolumeWeightedObject.get_volume_weighted_tds()
        for site in vwtds_data.copy():
            site_dict = vwtds_data[site]
            for entry in site_dict.copy():
                data = site_dict[entry]
                try:
                    time_string = entry.strftime('%Y-%m-%d %X')
                except:
                    time_string = entry
                site_dict[time_string] = data
                del site_dict[entry]
        x=1
        print('Volume Weighted TDS obtained')
    except Exception as e:
        print(f'Error in obtaining VWTDS -- Error: {e}')
        sys.exit

    print('Getting salt mass...')
    try:
        VolumeWeightedObject.get_salt_mass()
        print('Salt mass obtained')
        print('VW-TDS Successfully Calculated!')
    except Exception as e:
        print(f'Error in obtaining salt mass -- Error: {e}')
    result = {**vwtds_data, **regression_data}
    x=1
    return json.dumps(result)
