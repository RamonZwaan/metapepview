from pathlib import Path
from metapepview.constants import GlobalConstants
import json


def import_ref_statistics() -> str:
    statistics_dir = Path(GlobalConstants.reference_dataset_loc)
    
    if not statistics_dir.exists():
        raise ValueError("Ref stats directory location does not exist")

    combined_stats = dict()
    
    
    # iterate over all json files
    for stats_file in statistics_dir.glob("*.json"):
        filename = stats_file.stem
        
        if filename in combined_stats.keys():
            print("duplicate names encountered in reference statistics, skip...")
            continue
        
        # load contents into dict and add to combined dict
        contents = json.loads(stats_file.read_text())
        if check_valid_ref_statistics(contents) is False:
            print("invalid statistics file encountered...")
            continue
        combined_stats[filename] = contents
    
    # return combined dataset as string
    return json.dumps(combined_stats)


# function to check validity of performance ref file, not yet implemented
# TODO: Implement ref stats validity function
def check_valid_ref_statistics(file_contents: dict) -> bool:
    return True
