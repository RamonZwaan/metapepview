import io

import requests
import gzip

from typing import Tuple
from pathlib import Path

from constants import *

from backend.utils import check_file_presence


def check_gtdb_taxonomy_present(
    dir_loc: str | Path=GlobalConstants.gtdb_taxonomy_dir) -> Tuple[bool, bool]:
    """Check if gtdb taxonomy database is found in the specified location.

    Args:
        dir_loc (str, optional): Directory containing gtdb db files.
            Defaults to GlobalConstants.gtdb_taxonomy_dir.

    Returns:
        bool: True if files located in expected directory
    """
    gtdb_files = GlobalConstants.gtdb_taxonomy_files
    return check_file_presence(dir_loc, gtdb_files)


def download_gtdb_taxonomy(
    dir_loc: str | Path=GlobalConstants.gtdb_taxonomy_dir,
    remote_url: str=GlobalConstants.gtdb_taxonomy_url,
    overwrite: bool=False,
    create_parent_dirs: bool=False,
    fetch_zlib_archive: bool=True) -> Tuple[bool, str | None]:
    """Download latest version of gtdb taxonomy database into specified
    directory.

    Args:
        dir_loc (str | Path, optional): Parent dir for files to be stored in.
            Defaults to GlobalConstants.ncbi_taxonomy_dir.
        remote_url (str, optional): URL location to fetch remote database from.
            Defaults to GlobalConstants.ncbi_taxonomy_url.
        overwrite (bool, optional): Replace existing file with latest version.
            Defaults to False.
        create_parent_dirs (bool, optional): If write path location does not
            exist, create path. Defaults to False.
        fetch_zlib_archive (bool, optional): Download zlib compressed archive
            from server to be decompressed in memory. Defaults to True.
    
    Returns:
        Tuple[bool, str | None]: Success status, error message
    """
    dir_loc = Path(dir_loc)
    
    # if directory does not exist, create it
    if not (dir_loc.exists() and dir_loc.is_dir()):
        if create_parent_dirs is False:
            msg = "Specified Path does not exist. To create path, select 'Create parent dir'."
            return (False, msg)
        print(dir_loc.resolve().as_posix())
        dir_loc.mkdir(parents=True)

    # check if dir contains existing ncbi taxonomy
    gtdb_exist = check_gtdb_taxonomy_present(dir_loc)[1]
    
    if gtdb_exist is True and overwrite is False:
        msg = "(part of) GTDB Taxonomy dataset exists in dir location, select 'Overwrite' to continue"
        return (False, msg)
    
    # download tar into memory and extract
    for file_name in GlobalConstants.gtdb_taxonomy_files:
        # name file based on archive option
        remote_file = file_name + ".gz" if fetch_zlib_archive is True else file_name
        
        # fetch data from server
        response = requests.get(remote_url + f"/{remote_file}")
        if response.status_code != 200:
            msg = f"Failed to fetch remote data: Status code: {response.status_code}"
            return (False, msg)
        
        # pass in-memory archive into gzip file object for extraction
        if fetch_zlib_archive is True:
            extr_archive = gzip.decompress(response.content)
        else:
            extr_archive = response.content
            
        # store extracted file data
        if Path(dir_loc, file_name).exists() and overwrite is True:
            Path(dir_loc, file_name).unlink()
        try:
            extr_archive = extr_archive.decode()
            write_file = open(Path(dir_loc, file_name).resolve().as_posix(),
                                "wt")
            write_file.write(extr_archive)
            write_file.close()
        except:
            msg = "Failed to write data onto disk..."
            return (False, msg)
    
    # check that files are present in directory
    if not check_gtdb_taxonomy_present(dir_loc)[0]:
        msg = "failed to download gtdb taxonomy database..."
        return (False, msg)
    
    return (True, None)
