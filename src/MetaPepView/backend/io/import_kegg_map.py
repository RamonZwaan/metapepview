from typing import Tuple
from pathlib import Path
import json
import pandas as pd

from metapepview.constants import *
from metapepview.backend.utils import request_to_file, request_kegg_link_dbs, contents_to_file


def check_kegg_mapping_present() -> bool:
    """Check if mapping between kegg ko and orthology info is present in
    specified location.

    Args:
        dir_loc (str | Path, optional): Directory to check dataset presence.
            Defaults to GlobalConstants.kegg_map_dir.

    Returns:
        bool: True if file present, else False.
    """
    dir_loc = Path(GlobalConstants.kegg_map_dir)
    
    # check if directory exists
    if not dir_loc.is_dir():
        return False
    
    kegg_file_names = [
        GlobalConstants.kegg_brite_ko_file_name,
        GlobalConstants.kegg_brite_mod_file_name,
        GlobalConstants.kegg_pathway_list_file_name,
        GlobalConstants.kegg_module_list_file_name,
        GlobalConstants.kegg_ko_list_file_name,
        GlobalConstants.kegg_brite_list_file_name,
        GlobalConstants.kegg_mod_ko_link_name,
        GlobalConstants.kegg_path_ko_link_name,
        GlobalConstants.kegg_brite_ko_link_name,
        GlobalConstants.kegg_path_mod_link_name,
        GlobalConstants.kegg_ko_ec_link_name
    ]
    
    # name of file to check
    dir_files = [i.name for i in dir_loc.iterdir()]
    
    # return true if all files present
    return all(i in dir_files for i in kegg_file_names)
    

def download_kegg_ko_map(overwrite: bool=False,
                         create_parent_dirs: bool=False) -> Tuple[bool, str | None]:
    """Fetch ko to info map dataset from KEGG servers through use of its
    public api.

    Args:
        dir_loc (str | Path, optional): Directory to store file in.
            Defaults to GlobalConstants.kegg_map_dir.
        kegg_url (str | Path, optional): URL location to fetch remote db from.
            Defaults to GlobalConstants.kegg_map_url
        overwrite (bool, optional): Replace existing file with latest version.
            Defaults to False.
        create_parent_dirs (bool, optional): If write path location does not
            exist, create path. Defaults to False.
            
    Returns:
        Tuple[bool, str | None]: Success status, error message
    """
    dir_loc = Path(GlobalConstants.kegg_map_dir)
    
    if not dir_loc.is_dir() and create_parent_dirs is False:
        msg = "Specified Path does not exist. To create path, select 'Create parent dir'."
        return (False, msg)
        
    # fetch data from kegg server with api calls
    api_locs = [
        (GlobalConstants.kegg_brite_ko_url,
         GlobalConstants.kegg_brite_ko_file_name,
         _validate_brite_ko),
        (GlobalConstants.kegg_brite_mod_url,
         GlobalConstants.kegg_brite_mod_file_name,
         _validate_brite_module),
        (GlobalConstants.kegg_pathway_list_url,
         GlobalConstants.kegg_pathway_list_file_name,
         _validate_pathway_list),
        (GlobalConstants.kegg_module_list_url,
         GlobalConstants.kegg_module_list_file_name,
         _validate_module_list),
        (GlobalConstants.kegg_ko_list_url,
         GlobalConstants.kegg_ko_list_file_name,
         _validate_ko_list),
        (GlobalConstants.kegg_brite_list_url,
         GlobalConstants.kegg_brite_list_file_name,
         _validate_brite_list),
    ]
    
    # Fetch data from api location, store into file and validate data format
    for (url, file_name, validate_func) in api_locs:
        file_loc = Path(dir_loc, file_name)
        (success, msg) = request_to_file(url,
                                         file_loc,
                                         create_parent_dirs,
                                         overwrite)
        if success is False:
            return (success, msg)
        
        # if file is invalid format, remove file.
        (success, msg) = validate_func(file_loc)
        if success is False:
            file_loc.unlink(True)
            return (success, msg)

    # Fetch mappings between databases:
    link_files = [
        (GlobalConstants.kegg_mod_ko_link_name, 'ko', 'module'),
        (GlobalConstants.kegg_path_ko_link_name, 'ko', 'pathway'),
        (GlobalConstants.kegg_brite_ko_link_name, 'ko', 'brite'),
        (GlobalConstants.kegg_path_mod_link_name, 'module', 'pathway'),
        (GlobalConstants.kegg_ko_ec_link_name, 'enzyme', 'ko'),
    ]
    for file_name, target_db, source_db in link_files:
        link_file_loc = Path(dir_loc, file_name)
        (success, msg) = download_link_data(link_file_loc,
                                            target_db, #type:ignore
                                            source_db, #type:ignore
                                            create_parent_dirs,
                                            overwrite)
    
    return (True, None)


def download_link_data(file_loc: Path,
                       target_db: KeggDbs,
                       source_db: KeggDbs,
                       create_parent_dirs: bool,
                       overwrite: bool) -> Tuple[bool, str | None]:
    fetch_data = request_kegg_link_dbs(target_db, source_db)
    (success, msg) = contents_to_file(fetch_data, 
                     file_loc,
                     create_parent_dirs,
                     overwrite)
    if success is False:
        return (success, msg)
    
    # validate link file
    (success, msg) = _validate_mod_ko_link(file_loc)
    if success is False:
        file_loc.unlink(True)
        return (success, msg)
    
    return (True, None)
    

def _validate_brite_module(file: Path) -> Tuple[bool, str | None]:
    """Valida brite module dataset. It checks the format of the brite data 
    to be json and checks that the root name is as expected.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    with open(file) as contents:
        try:
            json_data = json.load(contents)
        except:
            return False, "Brite file not in json format..."
        brite_name = json_data.get("name")
        if brite_name is None:
            return False, "File data not a brite dataset..."
        elif brite_name != "ko00002":
            return False, "File data is not kegg module brite dataset..."
        else:
            return True, None


def _validate_brite_ko(file: Path) -> Tuple[bool, str | None]:
    """Valida brite ko dataset. It checks the format of the brite data to be
    json and checks that the root name is as expected.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    with open(file) as contents:
        try:
            json_data = json.load(contents)
        except:
            return False, "Brite file not in json format..."
        brite_name = json_data.get("name")
        if brite_name is None:
            return False, "File data not a brite dataset..."
        elif brite_name != "ko00001":
            return False, "File data is not kegg orthology brite dataset..."
        else:
            return True, None
    

def _validate_kegg_list_format(file: Path) -> Tuple[bool, str | None]:
    """Validata kegg list file format. It checks that file data corresponds
    to the expected output from a 'https://rest.kegg.jp/list/<db>' call.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    try:
        file_data = pd.read_csv(file,
                                sep="\t",
                                engine="python")
        if file_data.shape[0] == 0 or file_data.shape[1] != 2:
            return False, "KEGG list dataset of invalid shape..."
        return True, None
    except:
        return False, "Unable to parse KEGG list file"
    

def _validate_pathway_list(file: Path) -> Tuple[bool, str | None]:
    """Validate kegg pathway list file format.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    res, msg = _validate_kegg_list_format(file)
    
    if res == False:
        return False, "KEGG pathway dataset of invalid shape..."
    else:
        return True, None


def _validate_module_list(file: Path) -> Tuple[bool, str | None]:
    """Validate kegg module list file format.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    res, msg = _validate_kegg_list_format(file)
    
    if res == False:
        return False, "KEGG module list dataset of invalid shape..."
    else:
        return True, None


def _validate_ko_list(file: Path) -> Tuple[bool, str | None]:
    """Validate kegg ko list file format.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    res, msg = _validate_kegg_list_format(file)
    
    if res == False:
        return False, "KEGG ko list dataset of invalid shape..."
    else:
        return True, None


def _validate_brite_list(file: Path) -> Tuple[bool, str | None]:
    """Validate kegg brite list file format.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    res, msg = _validate_kegg_list_format(file)
    
    if res == False:
        return False, "KEGG brite list dataset of invalid shape..."
    else:
        return True, None


def _validate_mod_ko_link(file: Path) -> Tuple[bool, str | None]:
    """Validate kegg module to ko link file format.

    Args:
        file (Path): File location

    Returns:
        Tuple[bool, str | None]: True if valid format. Potential error
            message if validation failed.
    """
    res, msg = _validate_kegg_list_format(file)
    
    if res == False:
        return False, "KEGG link file in invalid format..."
    else:
        return True, None    
