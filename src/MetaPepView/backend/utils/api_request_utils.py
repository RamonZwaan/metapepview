from time import sleep
from pathlib import Path
from typing import List, Tuple, Dict, Any
import requests

import pandas as pd

from metapepview.constants import *


def fetch_request(url: str,
                  error_msg: str | None = None) -> requests.Response:
    req_get = requests.get(url)

    if error_msg is None:
        error_msg = "Request failed"
    error_msg += ": statuscode {}".format(req_get.status_code)
    
    if req_get.status_code != 200:
        raise RuntimeError(error_msg)
    return req_get


def request_unipept_pept_to_lca(pept_df: pd.DataFrame,
                                seq_col: str) -> pd.DataFrame:
    """From a dataset of peptides, fetch lca taxonomy and rank from unipept
    database.

    Args:
        pept_df (pd.DataFrame): Peptide dataset.
        seq_col (str): Column with peptide sequences.

    Returns:
        pd.DataFrame: Peptide sequences with LCA taxonomy and LCA rank
    """
    # perform global peptide search through unipept
    base_url = "http://api.unipept.ugent.be/api/v1/pept2lca.json?equate_il=true"

    # batch_size defines the amount of peptides to send to unipept in one request
    batch_size = 100
    seq_series = "&input[]=" + pept_df[seq_col].drop_duplicates()
    
    # initialize dataframe to store lca
    lca_df = pd.DataFrame(columns=[seq_col, "Global LCA", "Global LCA Rank"], dtype=object)

    # set index counter
    x = 0
    while True:
        if x+batch_size >= seq_series.size:
            proteins = [_ for _ in seq_series[x:]]
        else:
            proteins = [_ for _ in seq_series[x:x+batch_size]]

        # create url and send request to unipept
        req_str = "".join([base_url, *proteins])
        response = fetch_request(
            req_str,
            "Unipept request unsuccessfull"
            ).json()
    	
        lca_df = pd.concat(
            [
                lca_df,
                pd.DataFrame(
                    [(elem["peptide"], elem['taxon_id'], elem['taxon_rank']) for elem in response],
                    columns=[seq_col, "Global LCA", "Global LCA Rank"])
            ])
        x += batch_size

        # stop when end is reached
        if x >= seq_series.size:
            break
    return lca_df


def request_to_file(url: str,
                    file_loc: Path,
                    create_parent_dir: bool = False,
                    overwrite: bool = False) -> Tuple[bool, str | None]:
    """Fetch api response data and write contents to file.

    Args:
        url (str): Api URL.
        file_loc (Path): Path to file to write contents
        overwrite (bool, optional): Overwrite current file if exists.
            Defaults to False.

    Returns:
        Tuple[bool, str | None]: Success status and potential error message
    """
    try:
        response = fetch_request(url)
    except RuntimeError as err:
        return (False, repr(err))
    return contents_to_file(response.text, file_loc, create_parent_dir, overwrite)


def contents_to_file(contents: str,
                     file_loc: Path,
                     create_parent_dir: bool = False,
                     overwrite: bool = False) -> Tuple[bool, str | None]:
    """Write contents to file while checking for potential conflicts within
    the local file system, like overwrites or absent directories.

    Args:
        contents (str): File contents
        file_loc (Path): File location
        create_parent_dir (bool, optional): Create parent directories if
            they do not exist. Defaults to False.
        overwrite (bool, optional): Overwrite current file if exists.
            Defaults to False.

    Returns:
        Tuple[bool, str | None]: Success status and potential error message
    """
    # if file exist, either delete or stop function
    if file_loc.exists() and overwrite is True:
        file_loc.unlink()
    elif file_loc.exists() and overwrite is False:
        msg = "File exists, abort data retrieval. Select 'Overwrite' to continue"
        return (False, msg)

    # create directory if it does not exist (at this point, it is clear that
    # dirs can be created)
    if not file_loc.parent.is_dir() and create_parent_dir == True:
        file_loc.parent.mkdir(parents=True)
    elif not file_loc.parent.is_dir() and create_parent_dir == False:
        msg = "Parent directory does not exist, abort data retrieval."
        return (False, msg)
        
    # store response content in file
    with file_loc.open("w") as output_file:
        output_file.write(contents)

    return (True, None)


def request_kegg_link_entries(target_db: str,
                            db_entry: str | List[str]) -> str:
    """Fetch data from kegg using the link api method. Select the target
    database to retrieve data fron and specify entries to search for links
    in the target database. Results are returned as db_entry, target_db pairs.

    Args:
        target_db (str): Target_db to fetch links from
        db_entry (str | List[str]): Entries to match links from target_db

    Returns:
        str: List of entry, target_db link pairs.
    """
    # KEGG link url
    base_url = "https://rest.kegg.jp/link/"
    
    if isinstance(db_entry, str):
        db_entry = [db_entry]
    
    # run in batches, kegg limits to 100
    # pointer in query list
    response = ""
    batch_size = 50
    x = 0
    while True:
        if x + batch_size >= len(db_entry):
            query_elems = [_ for _ in db_entry[x:]]
        else:
            query_elems = [_ for _ in db_entry[x:x+batch_size]]
        
        # configure db entry string:
        query_elems[1:] = ["+" + x for x in query_elems[1:]]
        
        # create url and send request to unipept
        req_str = "".join([base_url, target_db, '/', *query_elems])
        response += fetch_request(
            req_str,
            "Failed to fetch KEGG link data"
            ).text
        
        # wait a bit before repeating api call
        sleep(0.2)
        
        x += batch_size

        # stop when end is reached
        if x >= len(db_entry):
            break
    
    # remove module and ko prefixes from link dataset
    return response.replace("md:", "").replace("ko:", "")
        
    
def request_kegg_link_dbs(target_db: KeggDbs,
                          source_db: KeggDbs) -> str:
    """Fetch data from kegg using the link api method. Select the target
    database to retrieve data fron and specify entries to search for links
    in the target database. Results are returned as db_entry, target_db pairs.

    Args:
        target_db (str): Target_db to fetch links from
        db_entry (str | List[str]): Entries to match links from target_db

    Returns:
        str: List of entry, target_db link pairs.
    """
    # KEGG link url
    url = f"https://rest.kegg.jp/link/{target_db}/{source_db}"
    
    response = fetch_request(
        url,
        "Failed to fetch KEGG link data"
        ).text
    
    # remove prefix
    t_pref = GlobalConstants.kegg_db_prefix_map[target_db]
    s_pref = GlobalConstants.kegg_db_prefix_map[source_db]
    
    # remove prefixes from link dataset
    for prefix in t_pref + s_pref:
        response = response.replace(prefix, "")
    return response
    
    

