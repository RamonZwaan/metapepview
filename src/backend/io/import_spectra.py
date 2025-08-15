import re
import base64
import zlib
import struct
from pathlib import Path
from datetime import datetime
from typing import IO, Sequence, Callable, Tuple, TypeVar, Dict, List, Any, overload, Literal
from xml.etree.ElementTree import Element
import xml.etree.ElementTree as ET

import pandas as pd
import numpy as np



# regex patterns
_SCAN_NUM_PATTERN = re.compile(r"(?<=scan=)[0-9]+")
_PEAKS_COUNT_PATTERN = re.compile(r"(?<=defaultArrayLength=)[0-9]+")


def mzxml_to_df(file: str | Path | IO[bytes],
                fields: Sequence[str] | None = None) -> pd.DataFrame:
    if isinstance(file, str):
        file = Path(file)
    
    # parse xml file
    tree = ET.parse(file)
    root = tree.getroot()
    
    # output data
    available_fields = [
        'num',
        'scanType',
        'centroided',
        'msLevel',
        'peaksCount',
        'polarity',
        'retentionTime',
        'collisionEnergy',
        'lowMz',
        'highMz',
        'basePeakMz',
        'basePeakIntensity',
        'totIonCurrent',
        'msInstrumentID',
        'peaks',
        'compressionType',
        'compressedLen',
        'byteOrder',
        'precision',
        'contentType',
        'precursorScanNum',
        'precursorIntensity',
        'precursorCharge',
        'activationMethod',
        'windowWideness',
        'precursorMz',
    ]
    
    ms2_only_cols = [
        'collisionEnergy',
        'precursorScanNum',
        'precursorIntensity',
        'precursorCharge',
        'activationMethod',
        'windowWideness',
        'precursorMz'
    ]
    
    # if a list of fields is supplied, create output dict with only those fields, else, create one with all fields
    # Note:
    #   Some fields are required, and therefore will always be present in the output.
    if fields is not None:
        if not all(i in available_fields for i in fields):
            raise ValueError("Invalid field name given. The following fields are valid:\n\n{}".format("\n".join(available_fields)))
        scans = {k: [] for k in fields}
        relevant_fields = set(fields)
    else:
        scans = {k: [] for k in available_fields}
        relevant_fields = set(available_fields)
    
    
    def check_add_field(key: str, element: ET.Element | str | None):
        if key not in relevant_fields:
            return
        elif not isinstance(element, ET.Element): 
            scans[key].append(element)
        else:
            scans[key].append(element.get(key))
    
    
    # parse scans
    for element in root.findall('{*}msRun/{*}scan'):
        if element.get('msLevel') == '1':
            peaks = element[0]
            for key in ms2_only_cols:
                check_add_field(key, None)
            
        else:
            precursor = element[0]
            peaks = element[1]
            
            for key in precursor.keys():
                check_add_field(key, precursor)
            check_add_field('precursorMz', precursor.text)
        
        for key in element.keys():
            # retention time is stored in ISO 8601 format ('PT...S'), currently, the prefix and suffix
            # are simply removed
            if key == 'retentionTime' and 'retentionTime' in relevant_fields:
                value = element.get(key)
                if value is not None:
                    value = value[2:-1]
                scans[key].append(value)
            else:
                check_add_field(key, element)

        for key in peaks.keys():
            check_add_field(key, peaks)
        check_add_field('peaks', peaks.text)

    return pd.DataFrame(scans)
    




newtype = TypeVar('newtype')

def mzml_to_df(file: str | Path | IO[bytes],
               fields: Sequence[str] | None = None) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Parse mzML spectral data file and fetch spectral data for storage into tabular
    format.

    Args:
        file (str | Path | IO[bytes]): mzML file data.
        fields (Sequence[str] | None, optional): fields to return from data.
            Defaults to None.

    Raises:
        ValueError: Incorrect data format encountered

    Returns:
        Tuple[pd.DataFrame, Dict[str, Any]]: Spectral data stored into dataframe,
            and metadata stored in separate dictionary.
    """
    # Overview of field names
    run_fields = (
        'timestamp',
        'spectrum count',
        'raw file name'
    )
    
    spectrum_fields = (
        'scan number',
        'centroided',
        'MS level',
        'peaks count',
        'polarity',
        'retention time',
        'collision energy',
        'low m/z',
        'high m/z',
        'base peak m/z',
        'base peak intensity',
        'ion injection time',
        'total ion current',
        'precursor m/z',
        'precursor intensity',
        'precursor charge',
        'precursor scan number',
        # 'precursor window wideness',
        'm/z array',
        'intensity array'
    )

    ms2_only_fields = (
        'collision energy',
        'precursor m/z',
        'precursor intensity',
        'precursor charge',
        'precursor scan num',
        # 'precursor window wideness',
    )

    available_fields = run_fields + spectrum_fields

    # data extractor function for each field
    exp_data_extractors: Dict[str, Callable[[ET.Element], Any]] = {
        'timestamp': lambda x: _fetch_timestamp(x),
        'spectrum count': lambda x: int(x.find('./{*}spectrumList').get('count')), #type:ignore
        'raw file name': lambda x: x.get('id'),
    }

    # assumes single elements in scan, precursor, etc.
    spectral_data_extractors: Dict[str, Callable[[Any], Any]] = {
        'scan number': lambda x: _get_scan_number(x['spectrum dict']),
        'centroided': lambda x: 'centroid spectrum' in x['spectrum dict'].keys(),
        'MS level': lambda x: _fetch_convert_data(x['spectrum dict'], 'ms level', int),
        'peaks count': lambda x: _fetch_convert_data(x['spectrum dict'], 'defaultArrayLength', int),
        'polarity': lambda x: _get_scan_polarity(x['spectrum dict']),
        'retention time': lambda x: _fetch_convert_data(x['scan dict'], 'scan start time', float),
        'collision energy': lambda x: _fetch_convert_data(x['precursor dict'], 'collision energy', float),
        'low m/z': lambda x: _fetch_convert_data(x['spectrum dict'], 'lowest observed m/z', float),
        'high m/z': lambda x: _fetch_convert_data(x['spectrum dict'], 'highest observed m/z', float),
        'base peak m/z': lambda x: _fetch_convert_data(x['spectrum dict'], 'base peak m/z', float),
        'base peak intensity': lambda x: _fetch_convert_data(x['spectrum dict'], 'base peak intensity', float),
        'ion injection time': lambda x: _fetch_convert_data(x['scan dict'], 'ion injection time', float),
        'total ion current': lambda x: _fetch_convert_data(x['spectrum dict'], 'total ion current', float),
        'precursor m/z': lambda x: _fetch_convert_data(x['precursor dict'], 'selected ion m/z', float),
        'precursor intensity': lambda x: _fetch_convert_data(x['precursor dict'], 'peak intensity', float),
        'precursor charge': lambda x: _fetch_convert_data(x['precursor dict'], 'charge state', int),
        'precursor scan number': lambda x: _get_precursor_scan_number(x['precursor dict']),
        # 'precursor window wideness': ...,
        'm/z array': lambda x: _fetch_convert_data(x['peaks dict']['m/z array'], 'binary'),
        'm/z array compression': lambda x: _fetch_convert_data(x['peaks dict']['m/z array'], 'compression type'),
        'm/z array binary type': lambda x: _fetch_convert_data(x['peaks dict']['m/z array'], 'binary type'),
        'intensity array': lambda x: _fetch_convert_data(x['peaks dict']['intensity array'], 'binary'),
        'intensity array compression': lambda x: _fetch_convert_data(x['peaks dict']['intensity array'], 'compression type'),
        'intensity array binary type': lambda x: _fetch_convert_data(x['peaks dict']['intensity array'], 'binary type')
        } # type:ignore

    if isinstance(file, str):
        file = Path(file)
    
    # parse xml file
    tree = ET.parse(file)
    root = tree.getroot()
    
    # if a list of fields is supplied, create output dict with only those fields, else, create one with all fields
    # Note:
    if fields is None:
        fields = spectrum_fields
    elif not all(i in available_fields for i in fields):
        raise ValueError("Invalid field name given. The following fields are valid:\n\n{}".format("\n".join(available_fields)))
    
    # setup dictionary to store all field data per spectrum
    metadata = dict()
    scans = {k: [] for k in fields}
    relevant_fields = set(fields)
    
    # fetch experiment metadata
    run_item = root.find('{*}mzML/{*}run')
    for field in exp_data_extractors.keys():
        metadata[field] = exp_data_extractors[field](run_item)

    # parse spectra
    update_metadata = False
    for element in root.findall('{*}mzML/{*}run/{*}spectrumList/{*}spectrum'):
        # setup data structure to store xml fetched data
        spec_data = dict()

        # fetch spectrum attributes and cvParam data
        spec_data["spectrum dict"] = _parse_element_data(element)
        spec_data["scan dict"] = _get_scan_data(element)
        spec_data["peaks dict"] = _get_peaks_data(element)

        # only assign precursor dict if ms2 scan given
        ms_level = spectral_data_extractors['MS level'](spec_data)
        if ms_level == 2:
            prec_dict = _get_precursor_data(element)
            spec_data["precursor dict"] = prec_dict
        else:
            spec_data["precursor dict"] = dict()

        # fetch data from dict and append to output dataset
        for field in fields:
            res = spectral_data_extractors[field](spec_data)
            scans[field].append(res)

        # some fields inside spectra are valid for complete experiment,
        # these are added to metadata only from the first spectrum
        if update_metadata is False:
            metadata['compression type'] = spectral_data_extractors['intensity array compression'](spec_data)
            metadata['binary type'] = spectral_data_extractors['intensity array binary type'](spec_data)
            metadata['byte order'] = 'little endian'
            update_metadata = True
    
    # format spectrum data into dataframe
    scan_data = pd.DataFrame(scans)
    
    # fetch total retention time by taking the last scan RT
    metadata['total retention time'] = scan_data.iloc[-1]['retention time']
    metadata['MS2 spectrum count'] = (scan_data['MS level'] == 2).sum()
    metadata['MS1 spectrum count'] = (scan_data['MS level'] == 1).sum()
    metadata['combined MS1 tic'] = scan_data[scan_data['MS level'] == 1]['total ion current'].sum()
    
    return (scan_data, metadata)


def _get_scan_polarity(spectrum: Dict) -> str | None:
    """Get polarity of scan from a given spectrum element

    Args:
        spectrum (Dict): Spectrum element

    Returns:
        str | None: 'positive' or 'negative' polarity. If no polarity field
            found, return None
    """
    if 'positive scan' in spectrum.keys():
        return 'positive'
    elif 'negative scan' in spectrum.keys():
        return 'negative'
    else:
        return None

def _get_scan_number(spectrum_elem: Dict | None) -> int | None:
    """Get scan number from a spectrum id element.

    Args:
        spectrum_id (Dict | None): Spectrum id

    Returns:
        int | None: scan number.
    """
    if spectrum_elem is None:
        return None
    
    spectrum_id = spectrum_elem['id']

    match = re.search(_SCAN_NUM_PATTERN, spectrum_id)
    if match is None:
        return None
    else:
        return int(match.group(0))

def _get_precursor_scan_number(precursor_elem: Dict | None) -> int | None:
    """Get scan number from precursor

    Args:
        precursor_elem (Dict | None): precursor element

    Returns:
        int | None: scan number.
    """
    if precursor_elem is None:
        return None
    
    precursor_id = precursor_elem.get('spectrumRef', None)
    if precursor_id is None:
        return None

    match = re.search(_SCAN_NUM_PATTERN, precursor_id)
    if match is None:
        return None
    else:
        return int(match.group(0))


def _fetch_timestamp(run_element: ET.Element) -> float | None:
    """Fetch experiment start timestamp from mzML run element

    Args:
        run_element (ET.Element): Run element block of mzML.

    Returns:
        datetime | None: Timestamp in datetime format if found.
    """
    timestamp_attr = run_element.get('startTimeStamp')
    if timestamp_attr is not None:
        timestamp_format = r'%Y-%m-%dT%H:%M:%SZ'
        # format timestamp: drop Z, then drop optional decimals, finally add Z again
        timestamp_attr = timestamp_attr[:-1].split(r".")[0] + 'Z'
        datetime_obj = datetime.strptime(timestamp_attr, timestamp_format)
        return datetime_obj.timestamp()
    return None


def _fetch_convert_data(dict_element: Dict,
                        field: str,
                        new_type: Callable[[str], newtype] | None = None) -> newtype | None:
    """Fetch data from a dict if present and convert into desired type.
    If no data found, return None.

    Args:
        dict_element (Dict): Element to fetch data from
        field (str): Dict key to fetch data from
        new_type (Callable[[str], newtype] | None): Type to convert data to.

    Returns:
        newtype | None: _description_
    """
    val = dict_element.get(field, None)
    if val is None or new_type is None:
        return val
    else:
        return new_type(val)


def _get_scan_data(spectrum: ET.Element,
                   allow_multiple: bool = False) -> Dict[str, str] | List[Dict[str, str]] | None:
    """Fetch data from scan elements of given spectrum element.

    Args:
        spectrum (Dict): Spectrum element
        allow_multiple (bool, Optional): If multiple scan elements in one spectrum,
            allow fetching data from each scan. Defaults to False.

    Returns:
        Dict[str, str] | None: Fetched scan data, if no scans stored, return None.
    """
    # can have only one or zero scanList elements
    scan_list = spectrum.find("./{*}scanList")
    if scan_list is None:
        return None
    
    scan_count = int(scan_list.get("count"))

    if allow_multiple is False and scan_count != 1:
        raise ValueError("Multiple scans encountered in Spectrum element during parsing, only one expected...")
    
    # fetch data from scan elements
    scan_data = []
    for scan in scan_list.findall("./{*}scan"):
        # rt = scan['scan start time']
        scan_params = _parse_element_data(scan)
        
        if allow_multiple is False:
            return scan_params
        else:
            scan_data.append(scan_params)
    
    return scan_data


def _get_precursor_data(spectrum: ET.Element,
                        allow_multiple: bool = False) -> Dict[str, str] | List[Dict[str, str]] | None:
    """Fetch data from precursor elements of given spectrum element.

    Args:
        spectrum (Dict): Spectrum element
        allow_multiple (bool, Optional): If multiple precursor elements in one spectrum,
            allow fetching data from each precursor. Defaults to False.

    Returns:
        Dict[str, str] | None: Fetched precursor data, if no precursors stored, return None.
    """
    # can have only one or zero precursorList elements
    precursor_list = spectrum.find("./{*}precursorList")
    if precursor_list is None:
        return None
    
    precursor_count = int(precursor_list.get("count"))

    if allow_multiple is False and precursor_count != 1:
        raise ValueError("Multiple precursors encountered in Spectrum element during parsing, only one expected...")
    
    # fetch data from precursor elements
    precursor_data = []
    for precursor in precursor_list.findall("./{*}precursor"):
        # parse attributes and cvParams
        precursor_params = _parse_element_data(precursor)
        
        # extend with params from sub elements, like isolation witdow and activation data
        precursor_params.update(
            _parse_element_data(precursor.find("./{*}activation"))
        )
        precursor_params.update(
             _parse_element_data(precursor.find("./{*}isolationWindow"))
        )
        
        if allow_multiple is False:
            precursor_params.update(_get_ion_data(precursor, allow_multiple))
            return precursor_params
        else:
            precursor_params['selected ions'] = _get_ion_data(precursor, allow_multiple)
            precursor_data.append(precursor_params)
        
    
    return precursor_data


def _get_ion_data(precursor: ET.Element,
                  allow_multiple: bool = False) -> Dict[str, str] | List[Dict[str, str]] | None:
    """Fetch data from selected ion elements of given precursor element.

    Args:
        precursor (Dict): Precursor element
        allow_multiple (bool, Optional): If multiple ion elements in one precursor,
            allow fetching data from each selected ion. Defaults to False.

    Returns:
        Dict[str, str] | None: Fetched ion data, if no ions stored, return None.
    """
    # can have only one or zero precursorList elements
    ion_list = precursor.find("./{*}selectedIonList")
    if ion_list is None:
        return None
    
    ion_count = int(ion_list.get("count"))

    if allow_multiple is False and ion_count != 1:
        raise ValueError("Multiple selected ions encountered in precursor element during parsing, only one expected...")
    
    # fetch data from precursor elements
    ion_data = []
    for ion_elem in ion_list.findall("./{*}selectedIon"):
        ion_params = _parse_element_data(ion_elem)
        
        if allow_multiple is False:
            return ion_params
        else:
            ion_data.append(ion_params)
    
    return ion_data


def _get_peaks_data(spectrum: ET.Element) -> Dict[str, str] | None:
    """Fetch peak mz and intensity data of given spectrum element.

    Args:
        spectrum (Dict): Spectrum element

    Returns:
        Dict[str, str] | None: Fetched binary. If no binary array, return None.
    """
    # can have only one or zero binaryDataArrayList elements
    binary_array_list = spectrum.find("./{*}binaryDataArrayList")
    if binary_array_list is None:
        return None
    
    # fetch data from binary elements
    binary_data = dict()
    for bin_array in binary_array_list.findall("./{*}binaryDataArray"):
        bin_params = _parse_element_data(bin_array)

        compression_options = ['no compression', 'zlib compression']
        binary_type_options = ['32-bit float', '64-bit float']

        if 'm/z array' in bin_params.keys():
            dict_key = 'm/z array'
        elif 'intensity array' in bin_params.keys():
            dict_key = 'intensity array'
        else:
            continue

        binary_data[dict_key] = {
            'binary': bin_array.find("./{*}binary").text,  
            'compression type': _return_if_present(compression_options, bin_params.keys()),
            'binary type': _return_if_present(binary_type_options, bin_params.keys())
        }

    return binary_data


def _return_if_present(options: Sequence,
                       key_array: Sequence) -> str | None:
    """Given a series of strings, return the string that is present in a
    given key array. The first match is returned

    Args:
        options (Sequence): Series of strings to match.
        key_array (Sequence): Array to match options to

    Returns:
        str: First option string that matched to the key array.
    """
    for val in options:
        if val in key_array:
            return val
    return None


def _parse_element_data(element: ET.Element) -> Dict[str, str]:
    """Fetch parameter values from mzML element, including sub elements.

    Args:
        element (ET.Element): ElementTree element object.

    Returns:
        Dict[str, str]: Dictionary of field names and values fetched.
    """
    # fetch attribute info
    elem_data = dict(element.items())

    # extend data with cvparam elements
    for param in element.findall("./{*}cvParam"):
        elem_data[param.get('name')] = param.get('value')

    return elem_data


def decode_mzxml_peaks(content: str | None,
                       peak_number: int,
                       compression_type: str='zlib',
                       byteorder: str="network",
                       precision: str="64") -> Tuple[np.ndarray, np.ndarray]:
    """Decode compressed binary array from mzxml into m/z and intensity
    numpy arrays.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Tuple of peaks.
            m/z and intensity array.
    """
    # if no value is present, return empty lists
    if content is None:
        return (np.array([]), np.array([]))
    # decode string to binary data
    binary = base64.b64decode(content.encode('ascii'))
    
    # decompress binary string
    if compression_type == 'zlib':
        binary = zlib.decompress(binary)
    
    # Define format of bytestring
    if byteorder == "network":
        b_order = "!"
    elif byteorder == "little endian":
        b_order = "<"
    else:
        b_order = "@"
    
    if precision == '64':
        p_precision = "d"
    else:
        p_precision = "f"    
    byte_format = "{o}{n}{p}".format(o=b_order, n=int(peak_number*2), p=p_precision)
    
    # unpack mz and int values from buffer
    peaks_tuple = struct.unpack(byte_format, binary)
    
    return (np.array(peaks_tuple[::2]), np.array(peaks_tuple[1::2]))


def decode_mzml_peaks(content: str | None,
                      peak_number: int,
                      compression_type: str='zlib compression',
                      byteorder: str="little endian",
                      precision: str="64-bit float") -> np.ndarray:
    """Decode compressed binary array from mzml into m/z and intensity
    numpy arrays.

    Returns:
        np.ndarray: Numpy array of peaks data.
    """
    # if no value is present, return empty lists
    if content is None:
        return np.array([])
    # decode string to binary data
    binary = base64.b64decode(content.encode('ascii'))
    
    # decompress binary string
    if compression_type == 'zlib compression':
        binary = zlib.decompress(binary)
    
    # Define format of bytestring
    if byteorder == "network":
        b_order = "!"
    elif byteorder == "little endian":
        b_order = "<"
    else:
        b_order = "@"
    
    if precision == '64-bit float':
        p_precision = "d"
    else:
        p_precision = "f"    

    byte_format = "{o}{n}{p}".format(o=b_order, n=int(peak_number), p=p_precision)
    
    # unpack mz and int values from buffer
    peaks_unpacked = struct.unpack(byte_format, binary)
    
    return np.array(peaks_unpacked)

    
def fetch_mzml_peaks_data(mzml_row: pd.Series,
                          peaks_dict: Dict[str, Dict[str, str]],
                          peaks_compression: str,
                          peaks_precision: str,
                          include_mz: bool = True,
                          include_int: bool = True) -> np.ndarray | Tuple[np.ndarray, np.ndarray]:
    """Return spectral peak arrays from mzml row and spectral peaks dict.

    Args:
        mzml_row (pd.Series): Row from mzml dataframe
        peaks_dict (Dict[str, Dict[str, str]]): Dictionary of peaks data.
        peaks_compression (str): Compression method for peaks data.
        peaks_precision (str): Precision of compression.
        include_mz (bool, optional): Return mz data from peaks dataset.
            Defaults to True.
        include_int (bool, optional): Return intensity data from peaks dataset.
            Defaults to True.

    Returns:
        np.ndarray | Tuple[np.ndarray, np.ndarray]: if both mz and intensity are
            fetched, a tuple containing mz data and int data is returned, else
            an array of mz data or intensity data is returned.
    """
    spectrum_peaks = peaks_dict[str(mzml_row.name)]
    
    if include_int is False and include_mz is False:
        raise ValueError("Neither intensity or mz data requested.")
    
    mz_arr, int_arr = np.ndarray([]), np.ndarray([])
    if include_mz is True:
        mz_arr = decode_mzml_peaks(spectrum_peaks["m/z array"],
                                   mzml_row['peaks count'],
                                   compression_type=peaks_compression,
                                   precision=peaks_precision)
        if include_int is False:
            return mz_arr
    if include_int is True:
        int_arr = decode_mzml_peaks(spectrum_peaks["intensity array"],
                                    mzml_row['peaks count'],
                                    compression_type=peaks_compression,
                                    precision=peaks_precision)
        if include_mz is False:
            return int_arr
        
    return (mz_arr, int_arr)
