"""
This module contains functions to deal with file data loaded into memory.

MetaPepView provides imported data as binary string data. To process import data
(which may be compressed archive data), several file-like operations need to be
done. The functions here allow quick passing of raw data provided by the 
dcc.Upload component to provide file-like or text-like data for further processing
"""

import base64
import pandas as pd
import io
import csv
from collections import defaultdict
from pathlib import Path
from typing import Type, Dict , Tuple, List, Any, IO, Sequence, Literal, TypeVar, overload
import tarfile
from zipfile import ZipFile
import zlib


def csv_to_dict(file_loc: Path | IO[str],
                keys_col: int = 0,
                values_col: int = 1,
                sep: str = ",",
                concat_values: bool = True) -> Dict[str, str] | Dict[str, List[str]]:
    """Convert csv or csv like file as dictionary, choosing one column as keys
    and another as values.

    Args:
        file_loc (Path | IO[str]): File location or contents
        keys_col (int, optional): Column index to take as keys. Defaults to 0.
        values_col (int, optional): Column index to take as values. Defaults to 1.
        sep (str: optional): Delimiter character. Defaults to ",".
        concat_values (bool, optional): append values as list to store
            multiple values under single key. Defaults to True. 

    Returns:
        Dict[str, str]: Dictionary.
    """
    if isinstance(file_loc, Path):
        file = open(file_loc, mode='r', newline='')
    else:
        file = file_loc

    reader = csv.reader(file, delimiter=sep)
    if concat_values is False:
        return {row[keys_col]: row[values_col] for row in reader}
    else:
        out_dict: Dict[str, List[str]] = defaultdict(list)
        for row in reader:
            out_dict[row[keys_col]].append(row[values_col])
        return out_dict 


def upload_contents_to_bytes(upload_contents: str) -> bytes:
    """Convert dash upload content into bytes string. In dash,
    content is uploaded as base64 encoded binary string.

    Args:
        upload_contents (str): Content uploaded through dash upload component.

    Returns:
        bytes: bytes array of data
    """
    
    # remove content type from the upload data
    content_type, content_string = upload_contents.split(',')
    return base64.b64decode(content_string)


def memory_to_file_like(upload_contents: str,
                        archive_format: str | None = None,
                        filename: str | None = None) -> IO[bytes]:
    """Convert in-memory encoded string into a file-like
    object for processing with functions that expect a file.

    Args:
        upload_contents (str): Content uploaded through dash upload component.

    Returns:
        io.BytesIO: file-like object
    """
    content = upload_contents_to_bytes(upload_contents)
    
    # either return bytes buffer object, or extract archive and return
    if archive_format is None:
        return io.BytesIO(content)
    else:
        return extract_in_memory_archive(io.BytesIO(content),
                                         archive_format,
                                         filename)


def memory_to_str(upload_contents: str,
                  archive_format: str | None = None,
                  filename: str | None = None) -> str:
    """Convert in-memory encoded string into string object for functions that
    require strings.

    Args:
        upload_contents (str): Content uploaded through dash upload component.

    Returns:
        str: Decoded data as string.
    """
    content = upload_contents_to_bytes(upload_contents)
    
    # either return bytes buffer object, or extract archive and return
    if archive_format is None:
        # convert bytes array to StringIO
        return content.decode('utf-8')
    else:
        extracted_content = extract_in_memory_archive(io.BytesIO(content),
                                                      archive_format,
                                                      filename)
        # convert BytesIO to StringIO-like object
        return io.TextIOWrapper(extracted_content, encoding='utf-8').read()


def compress_string(content: str,
                    compression_level: int = -1) -> str:
    """Perform zlib compression on string and return base64 encoded string
    representation of bytes object.

    Args:
        content (str): Input data
        compression_level (int, optional): Level of compression as managed in
            zlib.compress. -1 is the fastest compression, 0 is no compression,
            9 is the most compressed, but the slowest. Defaults to -1.

    Returns:
        bytes: compressed string as bytes object
    """
    content_bytes = content.encode(encoding='utf-8')
    compressed_bytes = zlib.compress(content_bytes, level=compression_level)
    return base64.b64encode(compressed_bytes).decode(encoding='utf-8')


def decompress_string(content: str) -> str:
    """Decompress base64 encoded zlib object into original data string

    Args:
        content (str): Zlib compressed and base64 encoded data.

    Returns:
        str: Decompressed string data.
    """
    bytes_data = base64.b64decode(content)
    decompressed = zlib.decompress(bytes_data)
    return decompressed.decode(encoding='utf-8')


def memory_to_stringio(upload_contents: str,
                       archive_format: str | None = None,
                       filename: str | None = None) -> IO[str]:
    """Convert in-memory encoded string into string buffer object for
    parsing of text data.

    Args:
        upload_contents (str): Content uploaded through dash upload component.

    Returns:
        io.StringIO: String buffer object
    """
    content = upload_contents_to_bytes(upload_contents)
    
    # either return bytes buffer object, or extract archive and return
    if archive_format is None:
        # convert bytes array to StringIO
        return io.StringIO(content.decode('utf-8'))
    else:
        extracted_content = extract_in_memory_archive(io.BytesIO(content),
                                                      archive_format,
                                                      filename)
        # convert BytesIO to StringIO-like object
        return io.TextIOWrapper(extracted_content, encoding='utf-8')


def import_csv_file(upload_contents: str | IO[str]) -> pd.DataFrame:
    """Convert in-memory encoded string into pandas dataframe.
    Requires data to be a csv type object

    Args:
        upload_contents (str): Content uploaded through dash upload component

    Returns:
        pd.DataFrame: Data stored in pandas DataFrame
    """
    # if data has been decoded, directly return csv DataFrame
    if isinstance(upload_contents, io.TextIOBase):
        return pd.read_csv(upload_contents)
    elif isinstance(upload_contents, str):
        return pd.read_csv(memory_to_stringio(upload_contents))
    else:
        raise TypeError("invalid content type supplied...")




def extract_in_memory_archive(file: IO[bytes] | Path,
                              archive_format: str,
                              filename: str | None = None) -> IO[bytes]:
    """Extract a file, or file-like object, archive into
    a string file-like object, to parse through.

    Args:
        file (io.BytesIO | Path): Path to archive or file-like
            binary of archive.
        archive_format (str): Archive format. Supported: {'.tar', '.zip'}.
        filename (str, optional): file inside archive to extract. If only one file
            in archive, no filename is required. In that case the single file
            is extracted. Defaults to "".
    Returns:
        io.StringIO: Extracted file as string file-like object.
    """
    
    # pass archive into tar- or zipfile object for extraction
    if archive_format == ".tar":
        if isinstance(file, Path):
            archive_obj = tarfile.open(name=file)
        elif isinstance(file, io.BytesIO):
            archive_obj = tarfile.open(fileobj=file)
        
        # if no file name given, the function can extract the single file
        # in the archive
        if filename is None:
            members = archive_obj.getmembers()
            if len(members) != 1:
                raise ValueError("Supplied archive contains multiple files, give filename to extract.")
            member = members[0]
        else:
            member = filename
        
        extracted_file = archive_obj.extractfile(member)
    
    elif archive_format == ".zip":
        archive_obj = ZipFile(file=file)
        
        # if no file name given, the function can extract the single file
        # in the archive
        if filename is None:
            members = archive_obj.infolist()
            if len(members) != 1:
                raise ValueError("Supplied archive contains multiple files, give filename to extract.")
            member = members[0]
        else:
            member = filename
        
        extracted_file =  archive_obj.open(member)
    else:
        raise TypeError("Supplied incorrect file type.")
    
    # for tarfiles, supplying directory as member gives empty file buffer.
    if extracted_file is None:
        raise TypeError("Supplied incorrect file type.")
    
    # return io.TextIOWrapper(extracted_file, encoding='utf-8')
    return extracted_file


def determine_archive_format(file_name: str) -> Literal['.tar', '.zip'] | None:
    """Function that returns the archive format for a given file name.
    supports '.tar', '.tar.gz', '.zip'

    Args:
        file_name (str): name of file

    Returns:
        str: archive format {'.tar', '.zip', None}
    """
    if file_name.endswith(".tar.gz") or \
        file_name.endswith("tar"):
        return ".tar"
    elif file_name.endswith(".zip"):
        return ".zip"
    else:
        return None


def archive_to_file_list(archive: str,
                         archive_format: Literal['.zip', '.tar']) -> Tuple[Sequence[IO[str]], List[str]]:
    """Extract list of files and names from archive file.

    Args:
        archive (str): Dash dcc.Upload contents of archive file
        archive_format (str): Format of archive file, determine with
            'determine_archive_format' function in utils.

    Returns:
        Tuple[List[IO[str]], List[str]]: List of files as string buffers,
            list of file names.
    """
    # decode data, do not extract, as an object is created separately
    archive_data = memory_to_file_like(archive)

    file_list, file_names = extract_all_in_memory_archive(archive_data,
                                                          archive_format)
    # convert files to StringIO
    file_list = [io.TextIOWrapper(file, encoding='utf-8') for file in file_list]

    return (file_list, file_names)

    
def extract_all_in_memory_archive(archive_file: IO[bytes] | Path,
                                  archive_format: Literal['.tar', '.zip']) -> Tuple[Sequence[IO[bytes]], List[str]]:
    """Given an archive file in '.zip' or '.tar.gz' format. Extract all files
    in the archive and return a list of file-like objects and a list of 
    file names. Folders within the archive are not supported, in that case
    a ValueError is given.

    Args:
        archive_file (IO[bytes] | Path): In-memory archive, or path to archive.
        archive_format (str): Format or archive file. Valid: {'.tar', '.zip'}.

    Raises:
        TypeError: Invalid archive format supplied.

    Returns:
        Tuple[List[IO[bytes] | None], List[str]]: Two lists, one list of 
            file-like objects. One list of file names.
    """
    # pass archive into tar- or zipfile object for extraction
    if archive_format == ".tar":
        if isinstance(archive_file, Path):
            archive_obj = tarfile.open(name=archive_file)
        elif isinstance(archive_file, io.BytesIO):
            archive_obj = tarfile.open(fileobj=archive_file)
        
        # get all members from tarfile
        members = archive_obj.getmembers()
        
        
        # extract each file in the archive
        extracted_files = [archive_obj.extractfile(mem) for mem in members]
        # filter None out (for type hint as it should be filtered already)
        extracted_files = [file for file in extracted_files if file is not None]
        
        file_names = archive_obj.getnames()

        # check presence of directories, should result in incompatible list lengths
        if len(extracted_files) != len(file_names):
            raise TypeError("Invalid archive content, directories are not supported")
    
    elif archive_format == ".zip":
        archive_obj = ZipFile(file=archive_file)
        
        # get all members from zipfile
        members = archive_obj.infolist()
        
        # check for presence of directories
        if any(mem.is_dir() for mem in members):
            raise TypeError("Directories in archive are not supported.")

        # extract files in the archive
        extracted_files =  [archive_obj.open(mem) for mem in members]
        file_names = [mem.filename for mem in members]
    else:
        raise TypeError("Supplied incorrect file type...")
    
    # return io.TextIOWrapper(extracted_file, encoding='utf-8')
    return (extracted_files, file_names)    


def check_file_presence(
    parent_dir: str | Path,
    file_list: Sequence[str]) -> Tuple[bool, bool]:
    """Check if all or any files from a given file list is present within
    a supplied parent directory.

    Args:
        parent_dir (str | Path): Parent directory.
        file_list (Sequence[str]): List of file names (No Path's).

    Returns:
        Tuple[bool, bool]: presence all files, presence any file
    """
    parent_dir = Path(parent_dir)

    if not parent_dir.exists():
        return (False, False)
    
    # get files for directory and match to expected file names
    dir_files = [i.name for i in parent_dir.iterdir()]

    
    # check presence all files
    if any(file not in dir_files for file in file_list):
        all_present = False
    else:
        all_present = True
        
    # check presence any file (required to know if new version is downloaded)
    if any(file in dir_files for file in file_list):
        any_present = True
    else:
        any_present = False
        
    return (all_present, any_present)
