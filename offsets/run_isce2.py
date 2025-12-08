import os
import shutil
import sys
import subprocess
import glob
import time
import re 
import getpass
from zipfile import ZipFile
import asf_search as asf
import numpy as np
import xarray as xr

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

    # # Extract ZIP files 
    # zips = glob.glob('SAFE/*.zip')
    # if zips:
    #     print('Extracting ZIP files')
    #     for zip_file in zips:
    #         print(f"  Extracting {zip_file}")
    #         with ZipFile(zip_file, 'r') as zf:
    #             zf.extractall('SAFE/')
            
    #         # Remove zip after extraction
    #         print(f"  Removing {zip_file}")
    #         os.remove(zip_file)
    
    #     print('All ZIP files extracted')

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

    bash_script = """
    #!/bin/bash
    URL=https://s1qc.asf.alaska.edu/aux_cal
    cd ../aux
    wget -r -l2 -nc -nd -np -nH -A SAFE $URL
    """
    # Run the script
    subprocess.run(bash_script, shell=True, executable="/bin/bash", check=True)

def setup_environment():
    """Setup ISCE2 environment"""
    conda_path = '/home/jovyan'
    isce_env = 'envs/isce2'
    
    isce_home = f'{conda_path}/{isce_env}/lib/python3.8/site-packages/isce'
    isce_stack = f'{conda_path}/{isce_env}/share/isce2'
    isce_lib_path = f"{conda_path}/{isce_env}/lib"
    
    os.environ["LD_LIBRARY_PATH"] = f"{isce_lib_path}:{os.environ.get('LD_LIBRARY_PATH', '')}"
    os.environ['ISCE_HOME'] = isce_home
    os.environ['ISCE_STACK'] = isce_stack
    os.environ['ISCE_ROOT'] = f'{conda_path}/{isce_env}/lib/python3.8/site-packages'
    
    # Update PATH and PYTHONPATH
    path_components = [
        f"{isce_home}/bin",
        f"{isce_home}/applications",
        f"{isce_stack}/topsStack",
        os.environ.get('PATH', '')
    ]
    os.environ['PATH'] = ':'.join(filter(None, path_components))
    
    pythonpath_components = [
        f'{conda_path}/{isce_env}/lib/python3.8/site-packages',
        isce_home,
        isce_stack,
        f"{isce_home}/applications",
        f"{isce_home}/components",
        os.environ.get('PYTHONPATH', '')
    ]
    os.environ['PYTHONPATH'] = ':'.join(filter(None, pythonpath_components))
    
    os.environ['OMP_NUM_THREADS'] = '8'
    
    return isce_home, isce_stack

# def find_safe_files():
#     """Find SAFE files in the safe directory"""
#     safe_dirs = ["SAFE", "SAFE", "SAFE"]
#     safe_dir = None
    
#     for dirname in safe_dirs:
#         if os.path.exists(dirname):
#             safe_dir = dirname
#             break
    
#     if not safe_dir:
#         return []
    
#     safe_files = []
#     for item in os.listdir(safe_dir):
#         if (item.endswith('.SAFE') or item.endswith('.zip')) and os.path.isdir(os.path.join(safe_dir, item)):
#             safe_files.append(os.path.join(safe_dir, item))
    
#     return sorted(safe_files)

def date_from_safe(file):
    """Extract acquisition date from SAFE file"""
    match = re.search(r'(\d{8})T\d{6}', file)
    return match.group(1)

def find_dem_file():
    """Find DEM file - use fixed path since DEM doesn't change"""
    # Fixed DEM path - update this to your actual DEM location
    fixed_dem_path = "/home/jovyan/crevasse-advection/offsets/dem/fixed_dem.wgs84"

    if os.path.exists(fixed_dem_path):
        return fixed_dem_path
    
    # Fallback to local search if fixed path doesn't exist
    dem_patterns = ["*.dem.wgs84", "*.dem", "*.tif", "*.tiff"]
    dem_files = []
    for pattern in dem_patterns:
        dem_files.extend(glob.glob(pattern))
    
    if dem_files:
        # Prefer .dem.wgs84 files as they usually work better with ISCE
        wgs84_dems = [f for f in dem_files if f.endswith('.dem.wgs84')]
        if wgs84_dems:
            return wgs84_dems[0]
        return dem_files[0]
    return None

def run_topsapp():
    # run environment setup
    isce_home, isce_stack = setup_environment()

    # Get parent dir
    parent_dir = os.getcwd()

    # Log file setup
    log_file_path = "topsapp_processing.log"

    with open(log_file_path, "w") as log_file:
        def log_and_print(message):
            print(message)
            log_file.write(message + "\n")
            log_file.flush()

        def log(message):
            log_file.write(message + "\n")
            log_file.flush()

        log_and_print(f"\nWorking directory: {os.getcwd()}")
        log_and_print(f"ISCE_HOME: {os.environ.get('ISCE_HOME')}")
        log_and_print(f"ISCE_STACK: {os.environ.get('ISCE_STACK')}")

        # Find SAFE files
        safe_files = sorted(glob.glob('SAFE/*.zip'))
        if len(safe_files) < 2:
            log_and_print(f"ERROR: Need at least 2 SAFE files, found {len(safe_files)}")
            return False

        # Find DEM
        dem_file = find_dem_file()
        if not dem_file:
            log_and_print("ERROR: No DEM file found")
            return False

        log_and_print(f"DEM file: {dem_file}")

        # Check if DEM has proper metadata files
        if dem_file.endswith('.dem'):
            xml_file = dem_file + '.xml'
            vrt_file = dem_file + '.vrt'
            log_and_print(f"Checking for DEM metadata files...")
            log_and_print(f"  XML file: {xml_file} - {'exists' if os.path.exists(xml_file) else 'missing'}")
            log_and_print(f"  VRT file: {vrt_file} - {'exists' if os.path.exists(vrt_file) else 'missing'}")
            
            # If no XML file exists, try to create a basic one for the DEM
            if not os.path.exists(xml_file):
                log_and_print("WARNING: DEM XML file missing. This may cause processing to fail.")
                log_and_print("Consider using a .dem.wgs84 file with proper metadata.")

        # Number of offset paits to compute (1 less than SAFE files)
        num_pairs = len(safe_files) - 1
        log_and_print(f"\nTotal number of offsets to compute: {num_pairs}")
        return_codes = []

        start_time = time.time()
        for i in range(num_pairs):
            log_and_print("\n")
            msg = f"=== Preprocessing for Offset {i+1} started at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
            log_and_print(f"{msg:=^80}")
            log_and_print(f"\nReference SAFE: {safe_files[i]}")
            log_and_print(f"Secondary SAFE: {safe_files[i+1]}")
            log_and_print(f"DEM file: {dem_file}")

            # ... Make directories for pairs ...
            dates = [date_from_safe(sf) for sf in safe_files[i:i+2]]
            pair_dir = f"{parent_dir}/{'-'.join(dates)}"
            os.makedirs(pair_dir, exist_ok=True)
            log_and_print(f"\n✓ Created directory: {pair_dir}")

            # ... change into that directory ...
            os.chdir(pair_dir)
            log_and_print(f"✓ Moved to directory: {pair_dir}")

            # Create reference catalog XML
            reference_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
    <component name="reference">
        <property name="safe">['../{safe_files[i]}']</property>
        <property name="output directory">reference</property>
        <property name="orbit directory">../orbits</property>
        <property name="polarization">hh</property>
    </component>'''
            
            # Create secondary catalog XML  
            secondary_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
    <component name="secondary">
        <property name="safe">['../{safe_files[i+1]}']</property>
        <property name="output directory">secondary</property>
        <property name="orbit directory">../orbits</property>
        <property name="polarization">hh</property>
    </component>'''

            # ... Get number of 12-day intervals between images
            start_dt, end_dt = dates
            start = np.datetime64(f'{start_dt[:4]}-{start_dt[4:6]}-{start_dt[6:]}')
            end = np.datetime64(f'{end_dt[:4]}-{end_dt[4:6]}-{end_dt[6:]}')
            n_passes = np.ceil((end - start).astype('int64') / 12)
            
            # Create main topsApp XML configuration for dense offsets
            xml_config = f'''<?xml version="1.0" encoding="UTF-8"?>
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
    </topsApp>'''
            
            # Write all XML files
            with open("reference.xml", "w") as f:
                f.write(reference_xml)
            log_and_print("✓ Created reference.xml")
            
            with open("secondary.xml", "w") as f:
                f.write(secondary_xml)  
            log_and_print("✓ Created secondary.xml")
            
            config_file = "topsApp.xml"
            with open(config_file, "w") as f:
                f.write(xml_config)
            log_and_print(f"✓ Created topsApp configuration: {config_file}\n")
            
            # Start TOPSAPP
            msg = f"=== topsApp for Pair {i+1} started at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
            log_and_print(f"{msg:=^80}")
            log_and_print(f"\nWorking Directory: {os.getcwd()}")
            log_and_print(f"Parent Directory:  {parent_dir}")
            
            # Run ISCE2 command
            cmd = f"conda run -n isce2 topsApp.py {config_file}"

            log_and_print(f"\nStarting topsApp with commmand: {cmd}")
            log_and_print("="*80 + '\n')

            try:
                process = subprocess.Popen(
                    cmd, 
                    shell=True, 
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    env=os.environ,
                    cwd=os.getcwd(),
                )

                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        message = output.strip()
                        log(message)

                return_code = process.poll()

                log_and_print("="*80)
                log_and_print(f"Processing finished at {time.strftime('%Y-%m-%d %H:%M:%S')}")
                log_and_print(f"Return code: {return_code}")
                
                if return_code == 0:
                    log_and_print("Cleaning Directory ...")
                    clean_dir(os.getcwd())
                    
                    log_and_print('\n')
                    log_and_print('Removing reference SAFE file:')
                    remove_safe(f'{parent_dir}/{safe_files[i]}')
                    log_and_print(f'Removed SAFE file: {safe_files[i]}\n')

                    msg = f"=== topsApp for Pair {i+1} complete at {time.strftime('%Y-%m-%d %H:%M:%S')} ==="
                    log_and_print(f"{msg:=^80}\n")
                    log_and_print(f"Offsets computed: {i+1} of {num_pairs}")
                    if i+1 < num_pairs: 
                        log_and_print(f"Continuing to pair {i+2}")
                    
                    return_codes.append(return_code)
                else:
                    log_and_print(f"TopsApp processing pair {i} failed with return code: {return_code}")
                    if i+1 < num_pairs: log_and_print(f"Continuing to pair {i+2}")
                    log_and_print(f"Press Ctrl + C to cancel")
                    return_codes.append(return_code)

            except Exception as e:
                log_and_print(f"Error during processing: {e}")
                return False

        if all((return_code == 0) for return_code in return_codes):
            log_and_print('Removing final SAFE file:')
            remove_safe(f'{parent_dir}/{safe_files[-1]}')
            log_and_print(f'Removed SAFE file: {safe_files[-1]}\n')

            log_and_print(f"Total Time: {round(time.time()-start_time)} seconds")
            return True

        else:
            return False

def crop_denseoff(file):
    da = xr.open_dataarray(file, engine='rasterio').rio.reproject(3031, nodata=np.nan)
    da = da.sel(x=slice(1.345e6, 1.43e6), y=slice(1.76e6, 1.646e6))
    return da

def clean_dir(processing_dir):
    keep_files = {
        "isce.log",
        "reference.xml",
        "secondary.xml",
        "topsApp.xml"
    }
    keep_dirs = {"merged"}

    # Clean rest of processing directory
    for entry in os.listdir():
        full_path = os.path.join(processing_dir, entry)

        if os.path.isdir(full_path):
            if entry not in keep_dirs:
                print(f"Deleting directory: {full_path}")
                shutil.rmtree(full_path)

        else:
            if entry not in keep_files:
                print(f"Deleting file: {full_path}")
                os.remove(full_path)

    # Remove DEM, reference and secondary SLCs from merged
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/dem*')]
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/reference*')]
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/secondary*')]
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/*rdr*')]
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/filt_dense_offsets.bil*')]

    print("Cropping Offset output files...")
    # crop dense offsets, save as netCDF, remove .bil*
    crop_denseoff(
        os.path.join(processing_dir, 'merged', 'dense_offsets.bil.geo')
    ).to_netcdf(os.path.join(processing_dir, 'merged', 'dense_offsets.nc'))
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/dense_offsets.bil*')]
    print("Offsets Cropped")

    print("Cropping SNR files...")
    # crop SNR, save as nc, remove .bil*
    crop_denseoff(
        os.path.join(processing_dir, 'merged', 'dense_offsets_snr.bil.geo')
    ).to_netcdf(os.path.join(processing_dir, 'merged', 'dense_offsets_snr.nc'))
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/dense_offsets_snr.bil*')]
    print("SNR Cropped")

    print("Cropping COV files")
    # crop COV, save as nc, remove .bil*
    crop_denseoff(
        os.path.join(processing_dir, 'merged', 'dense_offsets_cov.bil.geo')
    ).to_netcdf(os.path.join(processing_dir, 'merged', 'dense_offsets_cov.nc'))
    [os.remove(f) for f in glob.glob(f'{processing_dir}/merged/dense_offsets_cov.bil*')]
    print('COV Cropped')

def remove_safe(safe_path):
    if os.path.isdir(safe_path):
        shutil.rmtree(safe_path)

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
