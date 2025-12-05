import re
from pathlib import Path
from collections import defaultdict
from xml.etree.ElementTree import Element
import xml.etree.ElementTree as ET
from typing import Any, Callable, Dict, IO, Tuple

import pandas as pd
import numpy as np

from metapepview.constants import StyleConstants
from metapepview.html_templates import import_single_file
from metapepview.backend.utils import determine_archive_format, \
    memory_to_file_like, \
    compress_string


featurexml_types: Dict[str, Callable] = {
    'int': int,
    'float': float,
    'string': str
}

_SCAN_NUM_PATTERN = re.compile(r"(?<=scan=)[0-9]+")


def import_features(content: str, 
                    filename: str, 
                    mzml_metadata: Dict[str, Any] | None) -> Tuple[str, 
                                                                   Dict[str, Any],
                                                                   bool] | \
                                                             Tuple[None,
                                                                   None,
                                                                   bool]:
    """Import featurexml data uploaded to dashboard.

    Args:
        content (str): Feature data content.
        filename (str): Features file name.
        mzml_metadata (Dict[str, Any] | None): Metadata from mzml spectral file.

    Returns:
        Tupple[...]: Features datasets with valid indicator.
    """
    # only process features file after mzml file is imported
    if mzml_metadata is None:
        qa_box_style["background-color"] = StyleConstants.import_failed_color
        return (None, None)
    
    archive_format = determine_archive_format(filename)
    
    print("Process feature data...")
    try:
        data, metadata = featurexml_to_df(memory_to_file_like(content, archive_format), 
                                          None)
    except Exception as err:
        print(err)
        return (None, None, False)
    
    print("Finished feature processing...")
    return (compress_string(data.to_json()), metadata, True)


def featurexml_to_df(file: str | Path | IO[bytes],
                     mzml_name: str | None = None) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Extract feature data from featureXML and return dataset as dataframe.
    Optionally checks if mzml name in featurexml corresponds to user specified\
    mzml name.

    Args:
        file (str | Path | IO[bytes]): Feature data in featureXML format.
        mzml_name (str | None): mzml file name that featurexml should correspond
            to. Defaults to None.

    Raises:
        ValueError: Cannot parse the mzml file name.
        ValueError: Featurexml does not belong to specified mzml file.

    Returns:
        Tuple[pd.DataFrame, Dict[str, Any]]: Features dataset and metadata
    """
    if isinstance(file, str):
        file = Path(file)
    
    # parse xml file
    tree = ET.parse(file)
    root = tree.getroot()

    # get mzml_name
    feature_mzml_name = None
    for metadata_param in root.findall('{*}dataProcessing/{*}UserParam'):
        if metadata_param.get('name') == "parameter: in":
            feature_mzml_name = Path(metadata_param.get('value')).name
            break

    # check if mzml from features is same as imported mzml
    # if feature_mzml_name is None and mzml_name is not None:
    #     raise ValueError("Cannot resolve spectral data source from featureXML file...")
    # elif feature_mzml_name != mzml_name:
    #     raise ValueError(f"mzML source of features '{feature_mzml_name}' is different from imported mzML '{mzml_name}'...")

    # get number of features
    feature_num = int(root.find('{*}featureList').get('count'))

    # iterate over features and parse data
    feature_dict = defaultdict(list)
    for feature in root.findall('{*}featureList/{*}feature'):
        feature_dict = parse_feature(feature, feature_dict)
    
    dataset = pd.DataFrame(feature_dict)
    dataset = dataset.sort_values(by="retention time").reset_index(drop=True)

    total_feature_intensity = dataset['intensity'].sum()

    return (dataset,
            {"mzml file name": feature_mzml_name,
             "feature count": feature_num,
             "combined intensity": total_feature_intensity})


def parse_feature(feature: Element,
                  feature_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Parse xml feature element and return extracted data as dictionary.

    Args:
        feature (Element): Feature element.

    Returns:
        Dict[str, Any]: Extracted feature data.
    """
    # parse userparams: map names and values, and convert value to type specified from featureXML
    params = dict()
    for userparam in feature.findall('UserParam'):
        name = userparam.get('name')
        type_name = userparam.get('type')
        type_cls = featurexml_types.get(type_name)
        value = type_cls(userparam.get('value'))
        
        params[name] = value

    # get feature retention time and mono isotopic mass
    pos_data = [float(x.text) for x in feature.findall('position')]

    # add all parameters to feature dataset
    feature_dict["feature id"].append(feature.get('id'))
    feature_dict["feature label"].append(params.get('label'))
    feature_dict["scan number"].append(get_scan_number(params.get("spectrum_native_id")))
    feature_dict["quality"].append(float(feature.find('overallquality').text))
    feature_dict["retention time"].append(pos_data[0] / 60) # convert sec to min
    feature_dict["monoisotopic m/z"].append(pos_data[1])
    feature_dict["intensity"].append(float(feature.find("intensity").text))
    feature_dict["charge"].append(int(feature.find('charge').text))
    feature_dict["isotope count"].append(len(feature.findall("convexhull")))
    feature_dict["peak width"].append(params.get("FWHM") / 60) # convert sec to min
    feature_dict["score fit"].append(params.get("score_fit"))
    feature_dict["score correlation"].append(params.get("score_correlation"))

    return feature_dict


def get_scan_number(spectrum_id: str) -> int | float:
    """Return the scan number from a spectrum id string.
    It uses a regular expression to extract the scan number from the id string.
    If no match is found, None is returned.

    Args:
        spectrum_id (str): Spectrum id string.

    Returns:
        int | None: Scan number if matched from string.
    """
    match = re.search(_SCAN_NUM_PATTERN, spectrum_id)
    if match is not None:
        return int(match.group(0))
    else:
        return np.nan
