from pathlib import Path
from collections import defaultdict
import json
import io
import re
from typing import Self, Dict, Any, List, Tuple, Callable, Type, overload

import pandas as pd

from backend.utils import *
from backend.io.import_kegg_map import _validate_brite_ko,\
    _validate_brite_module,\
    _validate_pathway_list

class KeggDatabase:
    """KEGG database object to import and interact with kegg datasets
    """
    
    MODULE_NAME = re.compile(r"[a-zA-z\d]+")
    BRITE_PATHWAY_LIST = re.compile(r"(?<=\[).+(?=\])")
    BRITE_PATHWAY_ELEMS = re.compile(r"\S+")

    BRITE_EC_LIST = BRITE_PATHWAY_LIST
    BRITE_EC_ELEMS = re.compile(r"[\d]+\.[\d]+\.[\d]+\.\d")

    KO_NAME = MODULE_NAME
    BRITE_KO_PATHWAY = re.compile(r"(?<=\[PATH:).*(?=\])")
    BRITE_KO_BRITE = re.compile(r"(?<=\[BR:).*(?=\])")

    GENE_SYMBOL_LIST = re.compile(r"^(?:[^,\s]+(?:, )?)+(?=;)")
    GENE_SYMBOL_ELEMS = re.compile(r"[^,\s]+")
    

    def __init__(self,
                 pathway_dict: Dict[str, str],
                 module_dict: Dict[str, str],
                 ko_dict: Dict[str, str],
                 brite_dict: Dict[str, str],
                 pathway_to_module: Dict[str, List[str]],
                 pathway_to_ko: Dict[str, List[str]],
                 brite_to_ko: Dict[str, List[str]],
                 module_to_ko_dict: Dict[str, List[str]],
                 ko_to_module_dict: Dict[str, List[str]],
                 ko_to_ec_dict: Dict[str, List[str]],
                 ko_to_symbol_dict: Dict[str, str]):
        # Id description objects
        self.pathway_dict = pathway_dict
        self.module_dict = module_dict
        self.ko_dict = ko_dict
        self.brite_dict = brite_dict
        
        # mapping objects
        self.pathway_to_module_dict = pathway_to_module
        self.pathway_to_ko_dict = pathway_to_ko
        self.brite_to_ko_dict = brite_to_ko
        self.module_to_ko_dict = module_to_ko_dict
        self.ko_to_module_dict = ko_to_module_dict
        self.ko_to_ec_dict = ko_to_ec_dict
        self.ko_to_symbol_dict = ko_to_symbol_dict
    
    
    def list_pathways(self) -> List[str]:
        """Fetch list of all pathways.

        Returns:
            List[str]: list of pathways
        """
        return list(self.pathway_dict.keys())


    def list_modules(self, pathway: str | None = None) -> List[str]:
        """Fetch list of modules, with option to filter modules by pathway.

        Args:
            pathway (str | None): Pathway to fetch module from.

        Returns:
            List[str]: List of modules.
        """
        if pathway is None:
            return list(self.module_dict.keys())
        else:
            return self.pathway_to_module_dict.get(pathway, [])


    def list_brite(self) -> List[str]:
        """Fetch list of brite objects, with option to filter modules by pathway.

        Args:
            pathway (str | None): Pathway to fetch module from.

        Returns:
            List[str]: List of modules.
        """
        return list(self.brite_dict.keys())


    def list_ko(self,
                pathway: str | None = None,
                module: str | None = None,
                brite: str | None = None) -> List[str]:
        """Fetch list of orthologous groups, with option to filter modules by
        pathway or module.

        Args:
            pathway (str | None): Pathway to fetch module from.

        Returns:
            List[str]: List of modules.
        """
        if pathway is None and module is None and brite is None:
            return list(self.ko_dict.keys())
        elif module is not None:
            return self.module_to_ko_dict.get(module, [])
        elif pathway is not None:
            return self.pathway_to_ko_dict.get(pathway, [])
        elif brite is not None:
            return self.brite_to_ko_dict.get(brite, [])
        return []


    def module_to_ko(self, modules: str | List[str] | None) -> List[List[str]] | None:
        """Fetch KEGG Orthology id's from one or multiple module id's.

        Args:
            modules (str | List[str]): Module id's to search.

        Returns:
            List[str] | List[List[str]]: KEGG Orthology id's under modules.
        """
        if modules is None:
            return None
        elif isinstance(modules, str):
            modules = [modules]
        return [self.list_ko(module=mod) for mod in modules]
    
    def ko_to_module(self, ko: str) -> List[str]:
        """Fetch modules to which kegg orthology belongs under.

        Args:
            ko (str): keg orthology id.

        Returns:
            List[str] | None: module id's matched to ko
        """
        return self.ko_to_module_dict.get(ko, [])
        
    
    
    def ko_to_symbol(self, ko: str) -> str | None:
        """Return symbol for given kegg orthology.

        Args:
            ko (str): Kegg orthology id.

        Returns:
            str | None: Symbol name if present in dataset.
        """
        return self.ko_to_symbol_dict.get(ko, None)
   
    
    def ko_to_ec(self, ko: str) -> List[str]:
        """Return enzyme ec for given kegg orthology.

        Args:
            ko (str): Kegg orthology id.

        Returns:
            List[str] | None: ec codes if present in dataset.
        """
        return self.ko_to_ec_dict.get(ko, [])
    
    
    def to_json(self) -> str:
        """Write KeggDatabase object to json format string.
        """
        return json.dumps(
            {
                "pathway df": json.dumps(self.pathway_dict),
                "module df": json.dumps(self.module_dict),
                "ko df": json.dumps(self.ko_dict),
                "brite df": json.dumps(self.brite_dict),
                "pathway to module": json.dumps(self.pathway_to_module_dict),
                "pathway to ko": json.dumps(self.pathway_to_ko_dict),
                "brite to ko": json.dumps(self.brite_to_ko_dict),
                "module to ko": json.dumps(self.module_to_ko_dict),
                "ko to module": json.dumps(self.ko_to_module_dict),
                "ko to ec": json.dumps(self.ko_to_ec_dict),
                "ko to symbol": json.dumps(self.ko_to_symbol_dict)
            }
        )
    
    @classmethod
    def read_json(cls, json_data: str) -> Self:
        """Read json string representation of KeggDatabase object into class
        instance.

        Args:
            json_data (str): string containing class data in json format

        Returns:
            KeggDatabase: Class instance of KeggDatabase.
        """
        class_dict = json.loads(json_data)
        
        elems = [
            json.loads(class_dict["pathway df"]),
            json.loads(class_dict["module df"]),
            json.loads(class_dict["ko df"]), 
            json.loads(class_dict["brite df"]), 
            json.loads(class_dict["pathway to module"]),
            json.loads(class_dict["pathway to ko"]),
            json.loads(class_dict["brite to ko"]),
            json.loads(class_dict["module to ko"]),
            json.loads(class_dict["ko to module"]),
            json.loads(class_dict["ko to ec"]),
            json.loads(class_dict["ko to symbol"])
        ]
        
        return cls(*elems)
        

    @classmethod
    def from_files(cls,
                   pathways_file: Path,
                   module_file: Path,
                   ko_file: Path,
                   brite_file: Path,
                   module_ko_link_file: Path,
                   pathway_ko_link_file: Path,
                   brite_ko_link_file: Path,
                   pathway_module_link_file: Path,
                   ko_ec_link_file: Path) -> Self:
        # import pathways tsv file
        pathway_dict = csv_to_dict(pathways_file,
                                   sep="\t",
                                   concat_values=False)
        module_dict = csv_to_dict(module_file,
                                  sep="\t",
                                  concat_values=False)
        ko_dict = csv_to_dict(ko_file,
                              sep="\t",
                              concat_values=False)
        brite_dict = csv_to_dict(brite_file,
                                 sep="\t",
                                 concat_values=False)
        
        module_ko_dict = csv_to_dict(module_ko_link_file,
                                     sep="\t")
        ko_module_dict = csv_to_dict(module_ko_link_file,
                                     keys_col=1,
                                     values_col=0,
                                     sep="\t")
        
        pathway_ko_dict = csv_to_dict(pathway_ko_link_file,
                                      sep="\t")
        brite_ko_dict = csv_to_dict(brite_ko_link_file,
                                    sep="\t")
        pathway_module_dict = csv_to_dict(pathway_module_link_file,
                                          sep="\t")
        ko_ec_dict = csv_to_dict(ko_ec_link_file,
                                 sep="\t")
        ko_symbol_dict = cls._get_ko_protein_name_map(ko_dict) #type:ignore
        
        
        return cls(
            pathway_dict, #type:ignore
            module_dict, #type:ignore
            ko_dict, #type:ignore
            brite_dict, #type:ignore
            pathway_module_dict, #type:ignore
            pathway_ko_dict, #type:ignore
            brite_ko_dict, #type:ignore
            module_ko_dict, #type:ignore
            ko_module_dict, #type:ignore
            ko_ec_dict, #type:ignore
            ko_symbol_dict) 


    @classmethod
    def _get_ko_protein_name_map(cls,
                                 ko_list: Dict[str, str]) -> Dict:
        def extract_symbol(line: str) -> str | None:
            symbol = re.findall(cls.GENE_SYMBOL_LIST, line)
            if len(symbol) == 0:
                return None
            else:
                return symbol[0]
        
        out_dict = dict()
        for key, value in ko_list.items():
            symbol_str = extract_symbol(value)
            if symbol_str is not None:
                out_dict[key] = symbol_str

        return out_dict


    # @classmethod
    # def _parse_brite_ko(cls, json_data: Dict[str, Any]) -> Tuple[pd.DataFrame,
    #                                                              pd.DataFrame,
    #                                                              pd.DataFrame,
    #                                                              Dict[str, List[str]],
    #                                                              Dict[str, str]]:
    #     """Parse brite ko dataset and store ko mapping in dataframe.
        
    #     Args:
    #         json_data (Dict[str, Any]): Brite hierarchy dataset as dict.

    #     Returns:
    #         pd.DataFrame: ko mapping data.
    #     """
    #     pathway_ko_list = cls._element_parser(json_data,
    #                                           [],
    #                                           cls._fetch_ko_data)

    #     pathway_ko_dict = defaultdict(list)
    #     pathway_ec_dict = defaultdict(set)
    #     brite_ko_dict = defaultdict(list)
    #     ko_ec_dict = defaultdict(set)
    #     ko_symbol_dict = dict()
        
    #     # map data from parsed dataset to dicionaries
    #     for (pathway, brite, ec_list, ko, symbol) in pathway_ko_list:
    #         # rename pathway from ko*** to map****
    #         if pathway is not None:
    #             pathway = pathway.replace("ko", "map")
            
    #         if pathway is not None:
    #             pathway_ko_dict[pathway].append(ko)
    #         if pathway is not None and len(ec_list) > 0:
    #             pathway_ec_dict[pathway].update(ec_list)
    #         if brite is not None:
    #             brite_ko_dict[brite].append(ko)
    #         if len(ec_list) > 0:
    #             ko_ec_dict[ko].update(ec_list)
    #         if symbol is not None:
    #             ko_symbol_dict[ko] = symbol
                    
        
    #     path_df_dict = {'Pathway': list(pathway_ko_dict.keys()),
    #                     'KO': list(pathway_ko_dict.values())}
    #     path_ec_df_dict = {'Pathway': list(pathway_ec_dict.keys()),
    #                        'EC': list(pathway_ec_dict.values())}
    #     brite_df_dict = {'Brite': list(brite_ko_dict.keys()),
    #                      'KO': list(brite_ko_dict.values())}
        
    #     # in ec dict, convert sets to lists
    #     ko_ec_dict = {i: list(j) for i, j in ko_ec_dict.items()}
    #     # ec_df_dict = {'KO': list(ko_ec_dict.keys()),
    #     #               'EC': list(ko_ec_dict.values())}
    #     # symbol_df_dict = {'KO': list(ko_symbol_dict.keys()),
    #     #                   'Symbol': list(ko_symbol_dict.values())}
        
    #     return (pd.DataFrame(path_df_dict).explode('KO'),
    #             pd.DataFrame(path_ec_df_dict).explode('EC'),
    #             pd.DataFrame(brite_df_dict).explode('KO'),
    #             ko_ec_dict,
    #             ko_symbol_dict)
    
    # @classmethod
    # def _parse_brite_module(cls, json_data: Dict[str, Any]) -> pd.DataFrame:
    #     """Parse brite module dataset and store module mapping in dataframe.

    #     Args:
    #         json_data (Dict[str, Any]): Brite hierarchy dataset as dict.

    #     Returns:
    #         pd.DataFrame: Module mapping data.
    #     """
    #     pathway_module_list = cls._element_parser(json_data,
    #                                               [],
    #                                               cls._fetch_module_data)

    #     pathway_module_dict = defaultdict(set)
    #     for (pathway, module) in pathway_module_list:
    #         pathway_module_dict[pathway].update(set([module]))
        
        
    #     df_dict = {'Pathway': list(pathway_module_dict.keys()),
    #                'Module': list([list(x) for x in pathway_module_dict.values()])} 
    #     return pd.DataFrame(df_dict).explode('Module')
            
    
    # @classmethod
    # def _element_parser(cls,
    #                     brite_json: Dict[str, Any],
    #                     parent_names: List[str],
    #                     data_fetcher: Callable[[str, List[Any], List[str]], List[Any]]) -> List[Tuple[str, ...]]:
    #     """Parse KEGG Brite json hierarchy and fetch element data.

    #     Args:
    #         brite_json (Dict[str, Any]): List of brite elements

    #     Returns:
    #         List[Tuple[str, str]]: _description_
    #     """
    #     output_data = list()

    #     # get list of all modules
    #     brite_module_list = brite_json.get("children")
    #     if brite_module_list is None:
    #         raise ValueError("Invalid format of brite module file...")

    #     # parse through hierarchy and fetch module and path link
    #     for element in brite_module_list:
    #         if "children" in element.keys():
    #             parent_names.append(element['name'])
    #             output_data += cls._element_parser(element,
    #                                                parent_names,
    #                                                data_fetcher)
    #         # at end of hierarchy, parse and collect data 
    #         else:
    #             element_row = element['name']
    #             output_data = data_fetcher(element_row, output_data, parent_names)
        
    #     return output_data
    
    # @classmethod
    # def _fetch_module_data(cls,
    #                        row: str,
    #                        current_output: List[Any],
    #                        *args) -> List[Any]:
    #     """Fetch module name and pathway data from row in KEGG Brite module
    #     (ko00002) dataset.

    #     Args:
    #         row (str): Module row.
    #         current_output (List[Any]): Output dataset of all previously collected rows

    #     Raises:
    #         ValueError: Invalid format of input row.

    #     Returns:
    #         List[Any]: Output dataset updated with row data if valid format.
    #     """
    #     name = re.match(cls.MODULE_NAME, row)

    #     # do not update output if no name found
    #     if name is None:
    #         return current_output
    #     else:
    #         name = name.group()
            
    #     pathway_list = re.findall(cls.BRITE_PATHWAY_LIST, row)

    #     # need to find path links, else continue
    #     if len(pathway_list) == 0:
    #         return current_output
    #     elif len(pathway_list) == 1:
    #         pathway_list = pathway_list[0].strip('PATH:').split(' ')
    #     else:
    #         raise ValueError("Failed to correctly parse pathways in brite module file...")
        
    #     # regroup data to pathway: module mapping
    #     for pathway in pathway_list:
    #         current_output.append((pathway, name))
        
    #     return current_output

    # @classmethod
    # def _fetch_ko_data(cls,
    #                    row: str,
    #                    current_output: List[Any],
    #                    parent_rows: List[str]) -> List[Any]:
    #     """Fetch ko name and pathway data from row in KEGG Brite ko
    #     (ko00001) dataset.

    #     Args:
    #         row (str): KO row.
    #         current_output (List[Any]): Output dataset of all previously collected rows
    #         parent_rows (List[str]): Names of parent elements higher up in the hierarchy.

    #     Raises:
    #         ValueError: Invalid format of input row.

    #     Returns:
    #         List[Any]: Output dataset updated with row data if valid format.
    #     """
    #     name = re.match(cls.KO_NAME, row)
        
    #     # do not update output if no name found
    #     if name is None:
    #         return current_output
    #     else:
    #         name = name.group()
        
    #     # fetch symbol name
    #     symbol = re.findall(cls.GENE_SYMBOL_LIST, row)
    #     if len(symbol) == 0:
    #         symbol = None
    #     else:
    #         symbol = symbol[0]
        
    #     # fetch pathway from parent name  
    #     pathway = re.findall(cls.BRITE_KO_PATHWAY, parent_rows[-1])
    #     brite = re.findall(cls.BRITE_KO_BRITE, parent_rows[-1])
        
    #     ec_list = re.findall(cls.BRITE_EC_LIST, row)

    #     # need to find path links, else continue
    #     if len(ec_list) == 0:
    #         return current_output
    #     elif len(ec_list) == 1:
    #         ec_list_str = ec_list[0].strip('EC:')
    #         ec_list: List[str] = re.findall(cls.BRITE_EC_ELEMS, ec_list_str)
    #     else:
    #         raise ValueError("Failed to correctly parse pathways in brite module file...")

    #     # need to find path links, else continue
    #     out_pathway, out_brite = None, None
    #     if len(pathway) == 0 and len(brite) == 0:
    #         return current_output
    #     elif len(pathway) > 1 or len(brite) > 1:
    #         raise ValueError("Failed to correctly parse pathways in brite ko file...")

    #     # strip prefix from id's
    #     if len(pathway) == 1:
    #         out_pathway = pathway[0].strip('PATH:')
    #     if len(brite) == 1:
    #         out_brite = brite[0].strip('BR:')
         
    #     current_output.append((out_pathway, out_brite, ec_list, name, symbol))
        
    #     return current_output


    @staticmethod
    def _validate_brite_module(file: Path) -> Tuple[bool, str | None]:
        """Valida brite module dataset. It checks the format of the brite data 
        to be json and checks that the root name is as expected.

        Args:
            file (Path): File location

        Returns:
            Tuple[bool, str | None]: True if valid format. Potential error
                message if validation failed.
        """
        return _validate_brite_module(file)

    @staticmethod
    def _validate_brite_ko(file: Path) -> Tuple[bool, str | None]:
        """Valida brite ko dataset. It checks the format of the brite data to be
        json and checks that the root name is as expected.

        Args:
            file (Path): File location

        Returns:
            Tuple[bool, str | None]: True if valid format. Potential error
                message if validation failed.
        """
        return _validate_brite_ko(file)
        
    @staticmethod
    def _validate_pathway_list(file: Path) -> Tuple[bool, str | None]:
        """Validata pathway file format. It simply checks that the file can be
        parsed as a tsv file with two columns.

        Args:
            file (Path): File location

        Returns:
            Tuple[bool, str | None]: True if valid format. Potential error
                message if validation failed.
        """
        return _validate_pathway_list(file)
