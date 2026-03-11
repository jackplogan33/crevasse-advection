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

########### LOGGING ###########################################################

def get_logger(logger=None):
    return logger or logging.getLogger("topsapp")


def setup_logger(log_path="topsapp_processing.log"):
    logger = logging.getLogger("topsapp")
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.FileHandler(log_path))
    logger.addHandler(logging.StreamHandler())
    return logger


########### DATA ACQUISITION ##################################################

def get_SAFE(username, password):
    """Download SAFE Files"""
    # CreateASF session for downloads
    session = asf.ASFSession()
    try:
        user_pass_session = asf.ASFSession().auth_with_creds(username, password)
    except asf.ASFAuthenticationError as e:
        print(f'Auth failed: {e}')
    else:
        print('ASF Login Successful!')
        
    # Shirase AIO
    aoi = 'POLYGON((38.0336 -69.7358,38.0336 -70.4952,39.6985 -70.4952,39.6985 -69.7358,38.0336 -69.7358))'
    
    n_files = int(input('Number of files to download: '))
    if n_files:
        opts = {
            'platform':'S1',
            'start':str(input('Start Date (YYYY-MM-DD): ')),  # Input start date
            'processingLevel':'SLC',
            'frame':[830, 834, 936, 938, 939],
        }
        results = asf.search(intersectsWith=aoi, **opts)[-n_files:]
    
        print('Downloading SAFE files...')
        results.download(path='SAFE/', session=user_pass_session)
        print('SAFE files downloaded.\n')


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


########### ENVIRONMENT #######################################################

def get_envs_root():
    # Ask conda for its config info
    result = subprocess.run(
        ['conda', 'info', '--json'],
        capture_output=True, text=True
    )
    info = json.loads(result.stdout)
    return info['envs_dirs'][0]  # First/default envs directory

    
def setup_environment():
    """Setup ISCE2 environment"""
    # . . Automatically pull environment root
    envs_root = get_envs_root()

    isce_home = f"{envs_root}/envs/isce2/lib/python3.11/site-packages/isce"
    isce_stack = f"{envs_root}/envs/isce2/share/isce2"
    isce_lib_path = f"{envs_root}/envs/isce2/lib"
    
    os.environ['LD_LIBRARY_PATH'] = f"{isce_lib_path}:{os.environ.get('LD_LIBRARY_PATH', '')}"
    os.environ['ISCE_HOME'] = isce_home
    os.environ['ISCE_STACK'] = isce_stack
    os.environ['ISCE_ROOT'] = f"{isce_lib_path}/python3.11/site-packages"
    
    path_components = [
        f"{isce_home}/bin",
        f"{isce_home}/applications",
        f"{isce_stack}/topsStack",
        os.environ.get('PATH', '')
    ]
    os.environ['PATH'] = ':'.join(filter(None, path_components))
    os.environ['OMP_NUM_THREADS'] = '8'


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
    
    if not dem_file.endswith(".bil"):
        return
    for ext, label in [(".xml", "ISCE XML"), (".vrt", "VRT"), ('.aux.xml', 'XML')]:
        path = dem_file + ext
        status = "exists" if os.path.exists(path) else "missing"
        log.info(f"  DEM {label} file: {path} — {status}")
    if not os.path.exists(dem_file + ".xml"):
        log.warning("DEM XML file missing. This may cause processing to fail.")
        log.warning("Consider using a .dem.wgs84 file with proper metadata.")


########### XML GENERATION ####################################################

def make_scene_xml(role, safe_file, orbit_dir='../orbits'):
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<component name="{role}">
    <property name="safe">['../{safe_file}']</property>
    <property name="output directory">{role}</property>
    <property name="orbit directory">{orbit_dir}</property>
    <property name="polarization">hh</property>
</component>"""


def make_topsapp_xml(dem_file, n_passes):
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
    
            <!-- The swaths to process -->
            <property name="swaths">[2]</property>

            <!-- The region of interest -->
            <property name="region of interest">[-70.45441464, -69.88490745, 38.29553816,  39.83817435]</property>

            <!-- DEM for processing -->
            <property name="demFilename">{dem_file}</property>
    
            <!-- Unset all InSAR processing steps -->
            <property name="do interferogram">False</property>
            <property name="do ESD">False</property>
            <property name="do unwrap">False</property>
            <property name="do unwrap 2 stage">False</property>
            <property name="do ionosphere correction">False</property>
            <property name="geocode list">[]</property>
    
            <!-- Parameters for dense offsets -->
            <property name="do denseoffsets">True</property>
            <property name="Ampcor window width">256</property>
            <property name="Ampcor window height">64</property>
            <property name="Ampcor search window width">{30*n_passes}</property>
            <property name="Ampcor search window height">{10*n_passes}</property>
            <property name="Ampcor skip width">44</property>
            <property name="Ampcor skip height">8</property>

        </component>
    </topsApp>"""


def write_xml_files(ref_safe, sec_safe, dem_file, n_passes):
    """Write reference.xml, secondary.xml, and topsApp.xml to the current directory."""
    with open("reference.xml", "w") as f:
        f.write(make_scene_xml("reference", ref_safe))
    with open("secondary.xml", "w") as f:
        f.write(make_scene_xml("secondary", sec_safe))
    with open("topsApp.xml", "w") as f:
        f.write(make_topsapp_xml(dem_file, n_passes))


########### PROCESSING ####################################################

def run_topsapp_cmd(logger=None):
    """Run topsApp.py, stream output to logger, and return the exit code."""
    # . . Run topsApp to generate offsets
    log = get_logger(logger)
    cmd = ['conda', 'run', 'topsApp.py', 'topsApp.xml']
    
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

    
def process_pair(ref_safe, sec_safe, dem_file, parent_dir, logger=None):
    log = get_logger(logger)
    
    dates = [date_from_safe(sf) for sf in (ref_safe, sec_safe)]
    pair_dir = os.path.join(parent_dir, "-".join(dates))
    os.makedirs(pair_dir, exist_ok=True)
    os.chdir(pair_dir)
    log.info(f"  Working in: {pair_dir}")

    n_passes = compute_n_passes(*dates)
    write_xml_files(ref_safe, sec_safe, dem_file, n_passes)
    log.info("  Created reference.xml, secondary.xml, topsApp.xml")
    
    return run_topsapp_cmd(log)


def run_topsapp():
    # . . Setup steps
    setup_environment()
    log = setup_logger()
    parent_dir = os.getcwd()

    log.info(f"Working directory: {parent_dir}")
    log.info(f"ISCE_HOME:  {os.environ.get('ISCE_HOME')}")
    log.info(f"ISCE_STACK: {os.environ.get('ISCE_STACK')}")

    safe_files = sorted(glob.glob('SAFE/*.zip'))
    if len(safe_files) < 2:
        log.error(f"Need at least 2 SAFE files, found {len(safe_files)}")
        return False

    # Fixed DEM path
    dem_file = "../dem/REMA_10m_shirase.bil"  # Update to your DEM filename 
    log.info(f"DEM file: {dem_file}")
    validate_dem(dem_file, log)

    # Number of offsets to compute (1 less than number of images)
    num_pairs = len(safe_files) - 1
    log.info(f"Total pairs to process: {num_pairs}")

    # . . Begin processing image pairs
    return_codes = []
    start_time = time.time()
    for i, (ref, sec) in enumerate(zip(safe_files, safe_files[1:]), start=1):
        msg = f"=== Pair {i} started at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
        log.info(f"{msg:=^80}")
        log.info(f"Pair {i}/{num_pairs}")
        log.info(f"\tReference: {ref}")
        log.info(f"\tSecondary: {sec}")

        # . . Run topsApp command
        try:
            rc = process_pair(ref, sec, dem_file, parent_dir, log)
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
    xmin, xmax = 1.345e6, 1.43e6
    ymin, ymax = 1.646e6, 1.76e6
    da = xr.open_dataarray(file, engine='rasterio').rio.reproject(3031, nodata=np.nan)
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

    # Get EarthAccess credentials
    print('Earth Access login credentials')
    username = input('  Username: ')
    password = getpass.getpass('  Password: ')

    print('='*80)
    print(f"{'Get SAFE files':^80}")
    print('='*80)
    get_SAFE(username, password)

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


if __name__ == "__main__":
    sys.exit(main())