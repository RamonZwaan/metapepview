import io

from ftplib import FTP
import tarfile
import zipfile
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
    ncbi_archive = GlobalConstants.ncbi_taxonomy_archive
    
    # check presence of either files or archive
    return check_file_presence(dir_loc, [ncbi_archive]) # or \
        # check_file_presence(dir_loc, ncbi_archive)


def download_ncbi_taxonomy(
    dir_loc: str | Path=GlobalConstants.ncbi_taxonomy_dir,
    ftp_url: str=GlobalConstants.ncbi_taxonomy_url,
    overwrite: bool=False,
    create_parent_dirs: bool=False,
    extract_contents: bool=False) -> Tuple[bool, str | None]:
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
        print("Store NCBI taxonomy datasets in new directory: " + 
              dir_loc.resolve().as_posix())
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
    
    # if extract_contents is False:
    #     file_name = Path(dir_loc, archive_name)
    #     if file_name.exists() and overwrite is False:
    #         msg = "Output file exists... Provide other location or allow overwrite"
    #         return (False, msg)
        
    #     with open(Path(dir_loc, archive_name).resolve().as_posix(), "w+b") as out_file:
    #         ftp.retrbinary("RETR {}".format(archive_name), out_file.write)
            
    # download tar into memory and extract
    # else:
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
            msg = "Failed to extract data. Select tar or zip archive of database."
            return (False, msg)
        
        # extract only required files, remove existing files first'
        if extract_contents is True:
            for file_name in GlobalConstants.ncbi_taxonomy_files:
                if Path(dir_loc, file_name).exists() and overwrite is True:
                    Path(dir_loc, file_name).unlink(True)
                try:
                    archive_obj.extract(file_name, dir_loc.resolve().as_posix())
                except:
                    msg = "Failed to locate NCBI Taxonomy files in archive"
                    return (False, msg)
        # compress files back into archive for smaller storage
        else:
            # setup archive
            archive_loc = Path(dir_loc, GlobalConstants.ncbi_taxonomy_archive)

            if archive_loc.exists() and overwrite is True:
                Path(dir_loc, archive_loc).unlink(True)
            elif archive_loc.exists():
                raise FileExistsError(f"Archive '{archive_loc}' exists already...")
            
            try:
                with ZipFile(archive_loc, 'w') as out_archive:
                    # mem_buffers = [io.BytesIO() for x in range(3)]
                    for file_name in GlobalConstants.ncbi_taxonomy_files:
                        if isinstance(archive_obj, ZipFile):
                            bytes = archive_obj.read(file_name)
                        else:
                            io_bytes = archive_obj.extractfile(file_name)
                            if io_bytes is None: raise ValueError
                            bytes = io_bytes.read()
                        
                        # write contents to archive
                        out_archive.writestr(file_name, 
                                             bytes,
                                             compress_type=zipfile.ZIP_DEFLATED, 
                                             compresslevel=6)
            except:
                msg = "Failed to locate NCBI Taxonomy files in archive"
                # delete archive if write failed
                archive_loc.unlink(True)
                return (False, msg)
        
        archive_obj.close()
    
    ftp.quit()
    
    # check that files are present in directory
    if not check_ncbi_taxonomy_present(dir_loc)[0]:
        msg = "failed to download ncbi taxonomy database"
        return (False, msg)
    
    return (True, None)
