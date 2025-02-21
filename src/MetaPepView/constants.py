from enum import Enum
import re
from typing import List, Literal, Dict, Tuple
import plotly.express as px


# type literals 

RankType = Literal['Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
AbundanceMetricType = Literal['Match Count', 'Area']
KeggDbs = Literal['pathway', 'brite', 'module', 'ko', 'enzyme']


class PhysicalConstants:
    proton_mass = 1.00727

class GlobalConstants:
    # remote annotation databases
    standard_lineage_ranks: List[RankType] = ["Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    lineage_ranks_short: List[str] = ['K', 'P', 'C', 'O', 'F', 'G', 'S']
    ncbi_taxonomy_dir: str = r"./data/local/ncbi_taxonomy"
    ncbi_taxonomy_url: str = "ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    ncbi_taxonomy_files: Tuple[str, ...] = ("nodes.dmp", "names.dmp", "taxidlineage.dmp")
    
    gtdb_taxonomy_dir: str = r"./data/local/gtdb_taxonomy"
    gtdb_taxonomy_url: str = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest"
    gtdb_taxonomy_files: Tuple[str, ...] = ("bac120_taxonomy.tsv", "ar53_taxonomy.tsv")

    # lineage rank column suffix for global search
    global_annot_suffix: str = " (global search)"
    
    experiment_sample_table_cols: List[str] = ['Sample Name',
                                               'DB Search Imported',
                                               'De Novo Imported',
                                               'Taxonomy DB Name',
                                               'Functional Annotation DB Name']

    # kegg dataset locations
    kegg_map_dir: str = r"./data/local/kegg_datasets"
    kegg_map_url: str = "https://rest.kegg.jp/list/ko"
    kegg_map_file_name: str = "ko_mapping.tsv"
    kegg_brite_ko_url: str = "https://rest.kegg.jp/get/br:ko00001/json"
    kegg_brite_ko_file_name: str = "brite_ko.json"
    kegg_brite_mod_url: str = "https://rest.kegg.jp/get/br:ko00002/json"
    kegg_brite_mod_file_name: str = "brite_module.json"
    kegg_pathway_list_url: str = "https://rest.kegg.jp/list/pathway"
    kegg_pathway_list_file_name: str = "pathways.tsv"
    kegg_module_list_url: str = "https://rest.kegg.jp/list/module"
    kegg_module_list_file_name: str = "modules.tsv"
    kegg_ko_list_url: str = "https://rest.kegg.jp/list/ko"
    kegg_ko_list_file_name: str = "orthology_groups.tsv"
    kegg_brite_list_url: str = "https://rest.kegg.jp/list/brite"
    kegg_brite_list_file_name: str = "brite_groups.tsv"
    kegg_mod_ko_link_url: str = "https://rest.kegg.jp/link/ko/"
    kegg_mod_ko_link_name: str = "module_ko_map.tsv"
    kegg_path_ko_link_name: str = "pathway_ko_map.tsv"
    kegg_brite_ko_link_name: str = "brite_ko_map.tsv"
    kegg_path_mod_link_name: str = "pathway_module_map.tsv"
    kegg_ko_ec_link_name: str = "ko_ec_map.tsv"
    kegg_dropdown_desc_limit: int = 50
    kegg_db_prefix_map: Dict[str, List[str]] = {
        'pathway': ['path:'],
        'brite': ['br:'],
        'module': ['md:'],
        'ko': ['ko:'],
        'enzyme': ['ec:']
    }
    kegg_map_color_base_url: str = "https://www.kegg.jp/kegg-bin/show_pathway?"
    
    
    # Input datasets formats
    func_db_combine_delimiter: str = ";"                 # delimiter used during merging of rows within eggnog df
    func_db_combine_nan_fill: str = "-"                  # nan values need to be converted to str to be concatenated

    peptides_accession_delimiter: str = ";"
    
    # MetapepTable columns
    metapep_table_db_search_fields = ['PSM Count',
                                      'RT',
                                      'Scan',
                                      'm/z',
                                      'Charge',
                                      'ppm',
                                      'Length',
                                      'Feature Id',
                                      'Confidence',
                                      'Area',
                                      'Mass',
                                      'Accession',
                                      'PTM',
                                      'Source File']
    
    metapep_table_de_novo_fields = ['De Novo Confidence',
                                    'De Novo Area', 
                                    'De Novo Match Count', 
                                    'De Novo Scan', 
                                    'De Novo Source File']

    metapep_de_novo_fields = ['Confidence',
                              'Area', 
                              'Match Count', 
                              'Scan', 
                              'Source File']
    
    metapep_table_taxonomy_lineage = [i + ' Id' for i in standard_lineage_ranks] + \
                                     [i + ' Name' for i in standard_lineage_ranks]
    metapep_table_taxonomy_fields = ['Taxonomy Id',
                                     'Taxonomy Name'] + \
                                    metapep_table_taxonomy_lineage
    
    metapep_table_global_taxonomy_lineage = [i + f' Id (global search)' for i in standard_lineage_ranks] + \
                                            [i + f' Name (global search)' for i in standard_lineage_ranks]
    metapep_table_global_taxonomy_fields = [
        'Global LCA', 'Global LCA Rank'] + \
        metapep_table_global_taxonomy_lineage
    
    metapep_table_function_fields = ['KEGG_ko']


    db_search_dropdown_options = [
        {'label': 'Peaks 11', 'value': 'Peaks 11'},
        {'label': 'Peaks 10', 'value': 'Peaks 10'},
        {'label': 'MaxQuant', 'value': 'MaxQuant'},
        # {'label': 'Proteome Discoverer', 'value': 'ProteomeDiscoverer'},
        {'label': 'Sage', 'value': 'Sage'}
        ]
    de_novo_dropdown_options = [
        {'label': 'Peaks 11', 'value': 'Peaks 11'},
        {'label': 'Peaks 10', 'value': 'Peaks 10'},
        {'label': 'Novor (SearchGUI)', 'value': 'Novor'}
    ]
    
    # local datasets
    reference_dataset_loc: str = r"./data/share/qc_refs"
    crap_fasta_loc: str = "./data/share/crap.fasta"

    # annotation settings
    min_pept_len = 7
    trypsin_cleave_rule: re.Pattern[str] = re.compile(r"(?<=[RK])")
    sequence_regex: re.Pattern[str] = re.compile(r"[A-Z]+")
    global_annot_min_pept_tax = 3 # minimal number of unique peptides towards taxonomy to be included


class GraphConstants:
    # plotly color palette
    default_template = "plotly_white"
    color_palette = px.colors.qualitative.T10
    wide_color_palette = px.colors.qualitative.Dark24
    continuous_color_scale = px.colors.sequential.Blues
    #primary_color = "Black"
    primary_color = px.colors.qualitative.T10[0]
    gridcolor="slategray"
    secondary_grid_color="LightBlue"
    gridwidth=1
    default_layout = dict(
        margin=dict(l=20, r=20, t=10, b=10),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(
            family='Arial',
            size=14
        )
    )
    
    # Reference QA plots
    sample_trace_color = px.colors.qualitative.T10[2]



class StyleConstants:
    # Color styles
    plot_color = "#eeeeee"
    sidebar_color = "#dadada"
    header_color = "#483d45"
    tab_color = "#dddddd"

    import_failed_color = "#ffe5e5"
    import_success_color = "#ecffe6"

    # Common layout settings
    sidebar_width = "18rem"

    # define style arguments for the sidebar
    sidebar_style = {
        "position": "fixed",
        "top": 0,
        "left": 0,
        "bottom": 0,
        "width": sidebar_width,
        "padding": "2rem 1rem",
        "background-color": sidebar_color,
        "overflow-y": "auto",
        "overflow-x": "hidden",
        "display": "flex",
        "flex-direction": "column",
        "justify-content": "space-between"
    }

    # define style arguments for the header bar
    content_header_style = {
        "position": "fixed",
        "left": sidebar_width,
        "right": 0,
        "top": 0,
        "height": "4rem",
        "zIndex": 1, 
    }

    # define style arguments for the content box
    content_style = {
        "margin-top": "1rem",
        "margin-left": sidebar_width,
        "margin-right": "1rem",
        "padding": "0rem 0rem",
        "zIndex": 5
    }


    header_button_style = {
        "color": tab_color,
        "padding": "0 1rem",
        "cursor": "pointer",
        "display": "inline-block",
        "white-space": "nowrap"
    }


    # import container style
    # regular import container regular colour
    import_box_style = {"margin": ".5rem", "padding": "1rem 1rem",
        "background-color": plot_color,
        "border-radius": "1rem"}

    # regular import container success colour, if validation passed
    success_box_style = {"margin": ".5rem", "padding": "1rem 1rem",
        "background-color": import_success_color,
        "border-radius": "1rem"}

    # regular import container failed colour, if validation failed
    failed_box_style = {"margin": ".5rem", "padding": "1rem 1rem",
        "background-color": import_failed_color,
        "border-radius": "1rem"}
