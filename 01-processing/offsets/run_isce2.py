import os
import shutil
import sys
import subprocess
import glob
import time
import re 
import getpass
import asf_search as asf
import numpy as np
import xarray as xr
import logging
import json
from pathlib import Path

# Module level variables
cfg = {}

########### LOGGING ###########################################################

def get_logger(logger=None):
    return logger or logging.getLogger("topsapp")


def setup_logger(log_path):
    logger = logging.getLogger("topsapp")
    logger.setLevel(logging.DEBUG)
    
    file_handler = logging.FileHandler(log_path)
    file_handler.setLevel(logging.DEBUG)
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)  # suppresses DEBUG messages

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


########### DATA ACQUISITION ##################################################

def get_SAFE():
    """Download SAFE Files"""
    n_files = int(input('Number of files to download: '))

    if n_files == 0:
        return None, None
    
    print('Earth Access login credentials')
    username = input('  Username: ')
    password = getpass.getpass('  Password: ')    
    
    # CreateASF session for downloads
    session = asf.ASFSession()
    try:
        user_pass_session = asf.ASFSession().auth_with_creds(username, password)
    except asf.ASFAuthenticationError as e:
        print(f'Auth failed: {e}')
    else:
        print('ASF Login Successful!')
        
    # Shirase AIO
    aoi = cfg['scene']['aoi']    
    opts = {
        'platform':'S1',
        'start':str(input('Start Date (YYYY-MM-DD): ')),
        'processingLevel':'SLC',
        'frame':cfg['scene']['frames'],
    }
    results = asf.search(intersectsWith=aoi, **opts)[-n_files:]

    print('Downloading SAFE files...')
    results.download(path='SAFE/', session=user_pass_session)
    print('SAFE files downloaded.\n')

    return username, password


def get_orbits():
    safe_dir = '../SAFE/'
    safe_files = glob.glob(f'{safe_dir}*.zip')

    if len(safe_files) < 2:
        raise ValueError(f'You need at least 2 SAFE files, only found {len(safe_files)}')

    for file in safe_files:
        os.system(f'../fetchOrbit.py -i {file[8:-5]}')

    
def get_aux(ED_username, ED_password):
    os.system('echo "machine urs.earthdata.nasa.gov login {ED_username} password {ED_password}" > ~/.netrc')
    os.system('chmod 0600 ~/.netrc')

    bash_script = '''
    #!/bin/bash
    URL=https://s1qc.asf.alaska.edu/aux_cal
    cd ../aux
    wget -r -l2 -nc -nd -np -nH -A SAFE $URL
    '''
    # Run the script
    subprocess.run(bash_script, shell=True, executable='/bin/bash', check=True)


########### HELPERS ###########################################################

def date_from_safe(file):
    """Extract acquisition date from SAFE file"""
    match = re.search(r'(\d{8})T\d{6}', file)
    return match.group(1)


def compute_n_passes(start_dt, end_dt):
    """Calculate number of 12-day Sentinel-1 passes between two YYYYMMDD date strings."""
    to_dt64 = lambda s: np.datetime64(f"{s[:4]}-{s[4:6]}-{s[6:]}")
    return int(np.ceil((to_dt64(end_dt) - to_dt64(start_dt)).astype("int64") / 12))


def validate_dem(dem_file, logger=None):
    """Warn if expected DEM sidecar files are missing."""
    log = get_logger(logger)
    dem_file = cfg['scene']['dem_file']

    if not dem_file.endswith(".bil"):
        return

    log.info(f"DEM file: {dem_file}")
    for ext, label in [(".xml", "ISCE XML"), (".vrt", "VRT"), ('.aux.xml', 'XML')]:
        path = dem_file + ext
        status = "exists" if os.path.exists(path) else "missing"
        log.info(f"  DEM {label} file: {status}")
    
    if not os.path.exists(dem_file + ".xml"):
        log.warning("DEM XML file missing. This may cause processing to fail.")
        log.warning("Consider using a .bil file with proper metadata.")


def load_config(config_path="config.json"):
    """Load and validate user configuration."""
    global cfg
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(
            f"Config file '{config_path}' not found. "
            f"Copy config.example.json to config.json and fill in your values."
        )
    with open(config_path) as f:
        cfg = json.load(f)

    # . . Validate top-level keys
    required_top = ["scene", "insar", "ampcor", "crop", "compute"]
    missing_top = [k for k in required_top if k not in cfg]
    if missing_top:
        raise KeyError(f"Missing required config keys: {missing_top}")

    # . . Validate scene block
    required_scene = ["aoi", "frames", "dem_file", "swaths", "region_of_interest", "polarization"]
    missing_scene = [k for k in required_scene if k not in cfg["scene"]]
    if missing_scene:
        raise KeyError(f"Missing required 'scene' config keys: {missing_scene}")
    # Coerce dem_file to absolute path so it's correct everywhere
    cfg["scene"]["dem_file"] = os.path.abspath(cfg["scene"]["dem_file"])

    # . . Validate insar block
    required_insar = ["do_interferogram", "do_esd", "do_unwrap", "do_unwrap_2stage", "do_ionosphere", "geocode_list"]
    missing_insar = [k for k in required_insar if k not in cfg["insar"]]
    if missing_insar:
        raise KeyError(f"Missing required 'insar' config keys: {missing_insar}")

    # . . Validate ampcor block
    required_ampcor = ["do_denseoffsets", "window_width", "window_height", "skip_width", "skip_height",
                       "search_width_per_pass", "search_height_per_pass"]
    missing_ampcor = [k for k in required_ampcor if k not in cfg["ampcor"]]
    if missing_ampcor:
        raise KeyError(f"Missing required 'ampcor' config keys: {missing_ampcor}")

    # . . Validate crop block
    required_crop = ["xmin", "xmax", "ymin", "ymax", "output_epsg"]
    missing_crop = [k for k in required_crop if k not in cfg["crop"]]
    if missing_crop:
        raise KeyError(f"Missing required 'crop' config keys: {missing_crop}")

    # . . Validate compute block
    required_compute = ["omp_num_threads"]
    missing_compute = [k for k in required_compute if k not in cfg["compute"]]
    if missing_compute:
        raise KeyError(f"Missing required 'compute' config keys: {missing_compute}")

    return cfg


########### XML GENERATION ####################################################

def make_scene_xml(role, safe_file, orbit_dir='../orbits'):
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<component name="{role}">
    <property name="safe">['../{safe_file}']</property>
    <property name="output directory">{role}</property>
    <property name="orbit directory">{orbit_dir}</property>
    <property name="polarization">{cfg['scene']['polarization']}</property>
</component>"""


def make_topsapp_xml(n_passes):
    amp   = cfg['ampcor']
    insar = cfg['insar']
    scene   = cfg['scene']
    return f"""<?xml version="1.0" encoding="UTF-8"?>
    <topsApp>
        <component name="topsinsar">
            <property name="Sensor name">SENTINEL1</property>
    
            <!-- Scene XML files -->
            <component name="reference">
                <catalog>reference.xml</catalog>
            </component>
            <component name="secondary">
                <catalog>secondary.xml</catalog>
            </component>
    
            <property name="swaths">{scene['swaths']}</property>
            <property name="region of interest">{scene['region_of_interest']}</property>
            <property name="demFilename">{scene["dem_file"]}</property>
    
            <!-- InSAR processing steps -->
            <property name="do interferogram">{insar['do_interferogram']}</property>
            <property name="do ESD">{insar['do_esd']}</property>
            <property name="do unwrap">{insar['do_unwrap']}</property>
            <property name="do unwrap 2 stage">{insar['do_unwrap_2stage']}</property>
            <property name="do ionosphere correction">{insar['do_ionosphere']}</property>
            <property name="geocode list">{insar['geocode_list']}</property>
            
            <!-- Parameters for dense offsets -->
            <property name="do denseoffsets">{amp['do_denseoffsets']}</property>
            <property name="Ampcor window width">{amp['window_width']}</property>
            <property name="Ampcor window height">{amp['window_height']}</property>
            <property name="Ampcor search window width">{amp['search_width_per_pass'] * n_passes}</property>
            <property name="Ampcor search window height">{amp['search_height_per_pass'] * n_passes}</property>
            <property name="Ampcor skip width">{amp['skip_width']}</property>
            <property name="Ampcor skip height">{amp['skip_height']}</property>
            
        </component>
    </topsApp>"""


def write_xml_files(ref_safe, sec_safe, n_passes):
    """Write reference.xml, secondary.xml, and topsApp.xml to the current directory."""
    with open("reference.xml", "w") as f:
        f.write(make_scene_xml("reference", ref_safe))
    with open("secondary.xml", "w") as f:
        f.write(make_scene_xml("secondary", sec_safe))
    with open("topsApp.xml", "w") as f:
        f.write(make_topsapp_xml(n_passes))


########### PROCESSING ####################################################

def run_topsapp_cmd(logger=None):
    """Run topsApp.py, stream output to logger, and return the exit code."""
    # . . Run topsApp to generate offsets
    log = get_logger(logger)
    cmd = ['topsApp.py', 'topsApp.xml']
    
    process = subprocess.Popen(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    
    for line in iter(process.stdout.readline, ''):
        log.debug(line.strip())
    
    process.stdout.close()
    return process.wait()


def process_pair(ref_safe, sec_safe, parent_dir, logger=None):
    log = get_logger(logger)
    
    dates = [date_from_safe(sf) for sf in (ref_safe, sec_safe)]
    pair_dir = os.path.join(parent_dir, "-".join(dates))
    os.makedirs(pair_dir, exist_ok=True)
    os.chdir(pair_dir)
    log.info(f"\tWorking in: {pair_dir}")

    n_passes = compute_n_passes(*dates)
    write_xml_files(ref_safe, sec_safe, n_passes)
    log.info("\tCreated reference.xml, secondary.xml, topsApp.xml")
    
    return run_topsapp_cmd(log)


def run_topsapp():
    # . . Setup steps
    os.environ['OMP_NUM_THREADS'] = str(cfg['compute']['omp_num_threads'])
    log = setup_logger('./topsapp_processing.log')
    parent_dir = os.getcwd()

    log.info(f"Working directory: {parent_dir}\n")

    safe_files = sorted(glob.glob('SAFE/*.zip'))
    if len(safe_files) < 2:
        log.error(f"Need at least 2 SAFE files, found {len(safe_files)}")
        return False

    validate_dem(log)

    # Number of offsets to compute (1 less than number of images)
    num_pairs = len(safe_files) - 1
    log.info(f"Total pairs to process: {num_pairs}")

    # . . Begin processing image pairs
    return_codes = []
    start_time = time.time()
    for i, (ref, sec) in enumerate(zip(safe_files, safe_files[1:]), start=1):
        msg = f"=== Pair {i}/{num_pairs} started at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
        log.info(f"{msg:=^80}")
        log.info(f"\tReference: {ref}")
        log.info(f"\tSecondary: {sec}")

        # . . Run topsApp command
        try:
            rc = process_pair(ref, sec, parent_dir, log)
        except Exception as e:
            log.error(f"Pair {i} raised an exception: {e}")
            return False
            
        return_codes.append(rc)
        
        if rc == 0:
            log.info("Cleaning Directory ...")
            clean_dir(os.getcwd(), log)
            remove_safe(os.path.join(parent_dir, ref))
            log.info(f"\tRemoved reference file: {ref}")

            msg = f"=== Pair {i} complete at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
            log.info(f"{msg:=^80}")

            if i < num_pairs: 
                log.info(f"Continuing to pair {i+1}")
            
        else:
            msg = f"=== Pair {i} failed with return code: {rc} ==="
            log.warning(f"{msg:^80}")
            if i < num_pairs: 
                log.info(f"Continuing to pair {i+1}")
                log.info(f"Press Ctrl + C to cancel")

        os.chdir(parent_dir)

    if all(rc == 0 for rc in return_codes):
        remove_safe(os.path.join(parent_dir, safe_files[-1]))
        log.info(f'Removed SAFE file: {safe_files[-1]}')
    
        log.info(f"Total Time: {round(time.time()-start_time)}s")
        return True
    
    else:
        return False


########### POSTPROCESSING ####################################################

def crop_denseoff(file):
    crop = cfg['crop']
    xmin, xmax = crop['xmin'], crop['xmax']
    ymin, ymax = crop['ymin'], crop['ymax']
    da = xr.open_dataarray(file, engine='rasterio').rio.reproject(crop['output_epsg'], nodata=np.nan)
    return da.sel(x=slice(xmin, xmax), y=slice(ymax, ymin))


def clean_dir(processing_dir, logger=None):
    log = get_logger(logger)
    keep_files = ['isce.log', 'reference.xml', 'secondary.xml', 'topsApp.xml']
    keep_dirs  = ['merged']    

    processing_dir = Path(os.getcwd())
    
    # . . Remove unwanted files and directories
    for entry in processing_dir.iterdir():
        if entry.is_dir():
            if entry.name not in keep_dirs:
                log.info(f"Deleting directory: {entry}")
                shutil.rmtree(entry)
        else:
            if entry.name not in keep_files:
                log.info(f"Deleting file: {entry}")
                entry.unlink()
    
    # . . Remove large intermediates from merged/
    merged = processing_dir / 'merged'
    for pattern in ['dem*', 'reference*', 'secondary*', '*rdr*', 'filt_dense_offsets.bil*']:
        for f in merged.glob(pattern):
            f.unlink(missing_ok=True)

    # . . Crop dense offsets, SNR, and COV; save as netCDF and remove .bil* files
    stems  = ['dense_offsets', 'dense_offsets_snr', 'dense_offsets_cov']
    
    for stem in stems:
        log.info(f"Cropping {stem}...")
        crop_denseoff(merged / f'{stem}.bil.geo').to_netcdf(merged / f'{stem}.nc')
        for f in merged.glob(f'{stem}.bil*'):
            f.unlink()
        log.info(f"  Done.")

        
def remove_safe(safe_path):
    os.remove(safe_path)

########### ENTRY POINT #######################################################

def main():
    print("="*80)
    print(f"{'ISCE2 Run script':^80}")
    print("="*80)

    cfg = load_config()
    print('='*80)
    print(f"{'Get SAFE files':^80}")
    print('='*80)
    username, password = get_SAFE()

    print('='*80)
    print(f"{'Get Orbit files':^80}")
    print('='*80)
    os.chdir('orbits')
    get_orbits()
    os.chdir('..')

    # # . . Uncomment to enable AUX file download:
    # print('='*80)
    # print(f"{'Get AUX files':^80}")
    # print('='*80)
    # os.chdir('../aux/')
    # get_aux(username, password)

    print('='*80)
    print(f"{'Begin Offset Computation':^80}")
    print('='*80)
    success = run_topsapp()

    if success:
        print("\nProcessing completed successfully!")
        return 0
    else:
        print("\nProcessing failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())