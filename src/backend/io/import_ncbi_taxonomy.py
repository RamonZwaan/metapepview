import io

from ftplib import FTP
import tarfile
from zipfile import ZipFile

from typing import Tuple
from pathlib import Path

from constants import *

from backend.utils import check_file_presence


def check_ncbi_taxonomy_present(
    dir_loc: str | Path=GlobalConstants.ncbi_taxonomy_dir) -> Tuple[bool, bool]:
    """Check if ncbi taxonomy database is found in the specified location.

    Args:
        dir_loc (str, optional): Directory containing ncbi db files.
            Defaults to GlobalConstants.ncbi_taxonomy_dir.

    Returns:
        bool: True if files located in expected directory
    """
    ncbi_files = GlobalConstants.ncbi_taxonomy_files
    return check_file_presence(dir_loc, ncbi_files)


def download_ncbi_taxonomy(
    dir_loc: str | Path=GlobalConstants.ncbi_taxonomy_dir,
    ftp_url: str=GlobalConstants.ncbi_taxonomy_url,
    overwrite: bool=False,
    create_parent_dirs: bool=False) -> Tuple[bool, str | None]:
    """Download latest version of ncbi taxonomy database into specified
    directory.

    Args:
        dir_loc (str | Path, optional): Parent dir for files to be stored in.
            Defaults to GlobalConstants.ncbi_taxonomy_dir.
        ftp_url (str, optional): URL location to fetch remote database from.
            Defaults to GlobalConstants.ncbi_taxonomy_url.
        overwrite (bool, optional): Replace existing file with latest version.
            Defaults to False.
        create_parent_dirs (bool, optional): If write path location does not
            exist, create path. Defaults to False.
    
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
    ncbi_exist = check_ncbi_taxonomy_present(dir_loc)[1]
    
    if ncbi_exist is True and overwrite is False:
        msg = "(part of) NCBI Taxonomy dataset exists in dir location, select 'Overwrite' to continue"
        return (False, msg)
    
    # configure ftp connection
    split_url = ftp_url.split("/")
    host = split_url[0]
    path = "/".join(split_url[1:-1])
    archive_name = split_url[-1]
    
    try:
        ftp = FTP(host)
        ftp.login()
        ftp.cwd(path)
    except:
        msg = "Failed to connect to server..."
        return (False, msg)
    
    # download tar into memory and extract
    with io.BytesIO() as local_archive:
        try:
            ftp.retrbinary("RETR {}".format(archive_name), local_archive.write)
        except:
            msg = "Failed to fetch remote data..."
            return (False, msg)
        
        # after writing archive into object, set position to zero
        local_archive.seek(0)
        
        # pass in-memory archive into tarfile object for extraction
        if ".tar" in Path(archive_name).suffixes:
            archive_obj = tarfile.open(fileobj=local_archive)
        elif Path(archive_name).suffix == ".zip":
            archive_obj = ZipFile(local_archive)
        else:
            msg = "Failed to extract data. Select tar of zip archive of database."
            return (False, msg)
        
        # extract only required files, remove existing files first
        for file_name in GlobalConstants.ncbi_taxonomy_files:
            if Path(dir_loc, file_name).exists() and overwrite is True:
                Path(dir_loc, file_name).unlink()
            try:
                archive_obj.extract(file_name, dir_loc.resolve().as_posix())
            except:
                msg = "Failed to locate NCBI Taxonomy files in archive"
                return (False, msg)
        
        archive_obj.close()
    
    ftp.quit()
    
    # check that files are present in directory
    if not check_ncbi_taxonomy_present(dir_loc)[0]:
        msg = "failed to download ncbi taxonomy database"
        return (False, msg)
    
    return (True, None)