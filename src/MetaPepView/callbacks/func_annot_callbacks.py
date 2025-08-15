from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from MetaPepView.server import app

# import layout elements
from MetaPepView.html_templates import hidden_graph_with_text, \
    sample_color_table_block

from backend import *
from backend.utils import truncate_end
from backend.plots import pathway_abundance_barplot
from constants import GlobalConstants

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict


def ko_to_symbol(kegg_orthology: str | float,
                 ko_map_df: pd.DataFrame,
                 ko_set: set) -> str | float:
    if isinstance(kegg_orthology, float):
        return kegg_orthology
    
    if kegg_orthology is None:
        return np.nan
    
    if kegg_orthology not in ko_set:
        return np.nan
    
    ko_info = ko_map_df.loc[kegg_orthology]
    if ko_info.empty:
        return np.nan
    
    else:
        return ko_info[0].split(";")[0]


@app.callback(
    Output("kegg_db_class_data", "data"),
    Input("database_present_status", "data"),
    Input("kegg_db_class_data", "data"),
)
def kegg_db_to_memory(db_status, loaded_data):
    if db_status["kegg_map"] == True and loaded_data is None:
        # load file that maps KO to gene symbol
        gc = GlobalConstants
        pathway_list = Path(gc.kegg_map_dir, gc.kegg_pathway_list_file_name)
        module_list = Path(gc.kegg_map_dir, gc.kegg_module_list_file_name)
        ko_list = Path(gc.kegg_map_dir, gc.kegg_ko_list_file_name)
        brite_list = Path(gc.kegg_map_dir, gc.kegg_brite_list_file_name)
        brite_ko = Path(gc.kegg_map_dir, gc.kegg_brite_ko_file_name)
        module_ko = Path(gc.kegg_map_dir, gc.kegg_mod_ko_link_name)
        brite_module = Path(gc.kegg_map_dir, gc.kegg_brite_mod_file_name)

        mod_ko_link = Path(gc.kegg_map_dir, gc.kegg_mod_ko_link_name)
        path_ko_link = Path(gc.kegg_map_dir, gc.kegg_path_ko_link_name)
        brite_ko_link  = Path(gc.kegg_map_dir, gc.kegg_brite_ko_link_name)
        path_module_link = Path(gc.kegg_map_dir, gc.kegg_path_mod_link_name)
        ko_ec_link = Path(gc.kegg_map_dir, gc.kegg_ko_ec_link_name)
        
        kegg_db = KeggDatabase.from_files(pathway_list,
                                          module_list,
                                          ko_list,
                                          brite_list,
                                          mod_ko_link,
                                          path_ko_link,
                                          brite_ko_link,
                                          path_module_link,
                                          ko_ec_link)
        return kegg_db.to_json()
    elif loaded_data is not None:
        raise PreventUpdate
    
    return None

@app.callback(
    Output("kegg_ko_map_data", "data"),
    Input("database_present_status", "data"),
    State('kegg_map_loc', 'value')
)
def store_kegg_ko_map(db_status, kegg_db_loc):
    if db_status["kegg_map"] == True:
        # load file that maps KO to gene symbol
        __ko_mapping = Path(kegg_db_loc, GlobalConstants.kegg_map_file_name)

        if Path(__ko_mapping).exists():
            ko_map_df = pd.read_csv(__ko_mapping,
                                    sep='\t',
                                    engine="python", 
                                    names=["KO", "Info"])
            ko_map_df = ko_map_df.set_index("KO", drop=True)
            ko_map_df["symbol"] = ko_map_df["Info"].apply(lambda x: x.split(";")[0])
            return ko_map_df.to_json()
    return None


@app.callback(
    Output("brite_group_dropdown", "options"),
    Output("brite_group_dropdown", "disabled"),
    Input("kegg_db_class_data", "data")
)
def display_brite_groups(kegg_db):
    if kegg_db is None:
        return [], True
    
    kegg_db = KeggDatabase.read_json(kegg_db)
    label_lim = GlobalConstants.kegg_dropdown_desc_limit
    
    return [{'label': truncate_end(path + " - " + desc, label_lim),
             'search': desc,
             'title': desc,
             'value': path} for path, desc in kegg_db.brite_dict.items()], False


@app.callback(
    Output("pathway_dropdown", "options"),
    Output("pathway_dropdown", "disabled"),
    Input("kegg_db_class_data", "data")
)
def display_pathways(kegg_db):
    if kegg_db is None:
        return [], True
    
    kegg_db = KeggDatabase.read_json(kegg_db)
    label_lim = GlobalConstants.kegg_dropdown_desc_limit
    
    return [{'label': truncate_end(path + " - " + desc, label_lim),
             'search': desc,
             'title': desc,
             'value': path} for path, desc in kegg_db.pathway_dict.items()], False
    
    
@app.callback(
    Output("kegg_module_dropdown", "options"),
    Output("kegg_module_dropdown", "value"),
    Input("kegg_db_class_data", "data"),
    Input("pathway_dropdown", "value"),
)
def display_modules(kegg_db, pathway):
    if kegg_db is None:
        return [], None
    
    kegg_db = KeggDatabase.read_json(kegg_db)
    label_lim = GlobalConstants.kegg_dropdown_desc_limit

    module_list = kegg_db.list_modules(pathway)
    module_desc_list = [kegg_db.module_dict.get(mod, "-") for mod in module_list]
    
    return [{'label': truncate_end(mod + " - " + desc, label_lim),
             'search': desc,
             'title': desc,
             'value': mod} for mod, desc in zip(module_list, module_desc_list)],\
            None
 

 
@app.callback(
    Output("custom_pathway_items", "options"),
    Output("custom_pathway_items", "placeholder"),
    Output("custom_pathway_items", "disabled"),
    Input("kegg_db_class_data", "data")
)
def custom_proteins_options(kegg_db):
    placeholder_text = 'Select...'
    if kegg_db is not None:
        kegg_db = KeggDatabase.read_json(kegg_db)
        options = kegg_db.list_ko()
        return options, placeholder_text, False
    else:
        return [], "Import kegg map dataset...", True


@app.callback(
    Output("brite_group_dropdown_container", "hidden"),
    Output("pathway_dropdown_container", "hidden"),
    Output("kegg_module_container", "hidden"),
    Output("custom_pathway_container", "hidden"),
    Input("kegg_group_type", "value"),
)
def enable_pathway_selector(kegg_group_type):
    if kegg_group_type == "Pathway":
        return True, False, False, True
    elif kegg_group_type == "BRITE":
        return False, True, True, True
    elif kegg_group_type == "Manual":
        return True, True, True, False
    else:
        raise ValueError("invalid group selector type")
    

@app.callback(
    Output("pathway_normalize_sample", "options"),
    Output("pathway_normalize_sample", "disabled"),
    Input("peptides", "data")
)
def update_pathway_sample_normalization(peptide_json):
    if peptide_json is None:
        return [], True
    
    # convert dataset to pandas df
    peptide_df = MetaPepTable.read_json(peptide_json).data
    
    # get sample names and update dropdown menu
    return peptide_df['Sample Name'].unique().tolist(), False


@app.callback(
    Output("kegg_module_dropdown", "disabled"),
    Input("kegg_db_class_data", "data"),
    Input("kegg_display_format_radio", "value"),
)
def disable_module_dropdown(kegg_db, kegg_display_format):
    if kegg_display_format == "Module" or kegg_db is None:
        return True
    else:
        return False


@app.callback(
    Output("pathway_barplot_graph", "children"),
    Output("pathway_barplot_graph", "style"),
    Output('func_annot_figure_title', 'children'),
    Input('peptides', 'data'),
    Input("kegg_group_type", "value"),
    Input("kegg_display_format_radio", "value"),
    Input("brite_group_dropdown", "value"),
    Input("pathway_dropdown", "value"),
    Input("kegg_module_dropdown", "value"),
    Input("custom_pathway_items", "value"),
    Input("barplot_pathway_fraction_checkbox", "value"),    
    Input("barplot_pathway_include_taxa_checkbox", "value"),
    Input('barplot_func_quantification_column', 'value'),
    Input("pathway_normalize_sample", "value"),
    Input("func_annot_combine_duplicates", "value"),
    Input('tax_barplot_clade_selection_taxa', 'value'),
    Input('tax_barplot_clade_selection_rank', 'value'),
    State('kegg_db_class_data', 'data')
)
def update_pathway_barplot(peptide_json, 
                           kegg_group_method,
                           kegg_format,
                           predifined_brite_group,
                           predifined_pathway,
                           predifined_module,
                           custom_prot,
                           fractional_abundance,
                           include_taxa,
                           quant_method,
                           normalize_to_sample,
                           combine_annotations,
                           filter_clade,
                           clade_rank,
                           kegg_db):
    # set y axis column based on quantification method
    if quant_method == "Match Count":
        ycol = "PSM Count"
    else:
        ycol = "Area"


    title = "functional abundance"
    if kegg_db is not None:
        kegg_db = KeggDatabase.read_json(kegg_db)
    else:
        block_element = hidden_graph_with_text("pathway_barplot_figure",
                                               "Import KEGG dataset (sidebar)...")
        return block_element, dict(), "Figure"
    
    if peptide_json is None:
        block_element = hidden_graph_with_text("pathway_barplot_figure",
                                               "Import DB Search and functional annotation datasets...")
        return block_element, dict(), "Figure"
    
    metapep_obj = MetaPepTable.read_json(peptide_json)
    if metapep_obj.functional_annotation_present is False:
        block_element = hidden_graph_with_text("pathway_barplot_figure",
                                               "No samples with functional annotation in dataset...")
        return block_element, dict(), "Figure"
    
    peptide_df = metapep_obj.data
    
    # divide for each sample the psm value by the sum of psm for that sample, if specified
    if fractional_abundance is True:
        for sample in peptide_df["Sample Name"].unique():
            sample_idx = peptide_df[peptide_df["Sample Name"] == sample].index.to_list()
            peptide_df.loc[sample_idx, ycol] /= peptide_df.loc[sample_idx, ycol].sum() # type:ignore
    
    # if taxa selection made on the taxonomy barplot, filter protein abundances by taxa
    # get taxonomy id from barplot selection point
    if filter_clade and clade_rank and clade_rank != 'Root':
        peptide_df = filter_taxonomy_clade(peptide_df, filter_clade, clade_rank, 'Name')
        # set title name for graph    
        title += f" from {filter_clade} clade"
       
    # define column lists required during processing when taxonomy is included and when not 
    if include_taxa is True:
        cols_filter_df = [ycol, "Sample Name", "Taxonomy Name", "KEGG_ko"]              # use to filter unrelevant columns
        cols_group_peptides = ["Peptide Index", ycol, "Taxonomy Name", "Sample Name"]   # use to merge protein names by peptide sequence
        cols_group_names = ["Protein Name", "Taxonomy Name", "Sample Name"]                                   # use to group and count psm's by function and taxonomy
        tax_col = "Taxonomy Name"
    else:
        cols_filter_df = [ycol, "Sample Name", "KEGG_ko"]                             # use to filter unrelevant columns
        cols_group_peptides = ["Peptide Index", ycol, "Sample Name"]                  # use to merge protein names by peptide sequence
        cols_group_names = ["Protein Name", "Sample Name"]                                                  # use to group and count psm's by function
        tax_col = None
        
    # filter unrelevant columns out of the dataset    
    peptide_df = peptide_df[cols_filter_df]
    
    # split cells with multiple kegg annotations into separate rows
    peptide_df["KEGG_ko"] = peptide_df["KEGG_ko"].str.split(",") 

    # ensure unique index prior to explode, so that duplicate id's correspond to same original row
    peptide_df = peptide_df.reset_index(drop=True) 
    # put peptides annotated towards multiple kegg ko's in separate rows
    peptide_df = peptide_df.explode("KEGG_ko")
    
    # filter peptide dataset outside of predifined pathways/modules/brite
    if kegg_group_method == "Pathway":
        valid_ko = kegg_db.list_ko(pathway=predifined_pathway,
                                   module=predifined_module)
    elif kegg_group_method == "BRITE":
        valid_ko = kegg_db.list_ko(brite=predifined_brite_group)
    elif custom_prot is None or custom_prot == []:
        block_element = hidden_graph_with_text("pathway_barplot_figure",
                                                "Select protein terms in filter settings...")
        return block_element, dict(), "Figure"
    else:
        valid_ko = custom_prot
    peptide_df = peptide_df[peptide_df["KEGG_ko"].isin(valid_ko)]
    
    # filter dataset with annotations belonging to nitrogen cycle
    if kegg_format == "Protein/Gene name":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_symbol)
        title = "Protein " + title
    elif kegg_format == "EC":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_ec)
        title = "Enzyme " + title
    elif kegg_format == "Module":
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"].apply(kegg_db.ko_to_module)
        title = "KEGG module " + title
    else:
        peptide_df["Protein Name"] = peptide_df["KEGG_ko"]
        title = "KEGG orthology " + title
        
    # one ko may correspond to multiple modules or ec, explode to separate rows
    peptide_df = peptide_df.explode("Protein Name")
    # ko under pathway may also match towards modules outside of pathway, filter these
    if kegg_format == "Module" and predifined_pathway is not None:
        peptide_df[peptide_df["Protein Name"].isin(kegg_db.list_modules(predifined_pathway))]
    
    peptide_df.dropna(subset="Protein Name", inplace=True)
    
    # drop kegg ko column
    peptide_df.drop("KEGG_ko", axis=1, inplace=True)
    
    # merge duplicate peptides whose protein name are the same
    peptide_df["Peptide Index"] = peptide_df.index
    peptide_df.drop_duplicates(subset=["Peptide Index", "Protein Name"], inplace=True)
    
    # for peptide duplicates with different protein names there are two options:
    #   Count psm towards both names: this does result in inflated PSM numbers
    #   Combine names (alphabetically): Most realistic data, but results in more categories that may not show in the plot.
    if combine_annotations is True:
        peptide_df = peptide_df.groupby(by=cols_group_peptides)["Protein Name"]\
            .agg(lambda x: ",".join(sorted(x)))
        peptide_df = peptide_df.to_frame().reset_index(names=cols_group_peptides)
        
    # with combined protein name categories, group by these categories to obtain psm counts
    peptide_df = peptide_df.groupby(by=cols_group_names)[ycol].agg('sum')
    peptide_df = peptide_df.to_frame().reset_index(names=cols_group_names)

    # keep only top n ko's based on largest contribution across samples
    top_n_func = peptide_df.groupby(by="Protein Name")[ycol]\
        .agg('sum')\
        .sort_values(ascending=False)\
    
    # check if elements exceeds limit, set in plot title
    if top_n_func.shape[0] > 20:
        title += " (top 20 abundant)"
        top_n_func = top_n_func.head(20)
    
    peptide_df = peptide_df[peptide_df["Protein Name"].isin(set(top_n_func.index))]

    # rescale each protein type to normalize towards one specified sample
    if normalize_to_sample is not None:
        peptide_df = quantity_to_fold_change(peptide_df,
                                             normalize_to_sample,
                                             'Sample Name',
                                             'Protein Name',
                                             ycol)
        # drop all zero values
        peptide_df = peptide_df[peptide_df[ycol] > 0.0]
    
    # configure y-axis title
    ytitle = "Peptide spectrum matches" if ycol == "PSM Count" else "Signal area"
    if fractional_abundance is True:
        ytitle = "Fraction " + ytitle.lower()
    if normalize_to_sample is not None:
        ytitle += f"/ {normalize_to_sample}"
    
    plot = pathway_abundance_barplot(peptide_df,
                                    "Sample Name",
                                    ycol,
                                    "Protein Name",
                                    tax_col=tax_col,
                                    custom_title=title,
                                    custom_xname=kegg_format,
                                    custom_yname=ytitle)
    
    plot.update_layout(title="")

    return dcc.Graph(figure=plot,
                     id="pathway_barplot_figure",
                     style={"height": "40rem"}), dict(), title


@app.callback(
    Output("kegg_export_samples", "options"),
    Input('peptides', 'data')
)
def update_kegg_export_sample_selector(peptide_json):
    """Update sample to visualize in the db search vs de novo
    taxonomy comparison.
    """
    if peptide_json is None:
        return []
    
    # TODO: This data should be in peptide metadata table
    metapep_obj = MetaPepTable.read_json(peptide_json)
    metapep_df = metapep_obj.data

    # get sample names and if they have de novo taxonomy annotation
    metapep_df = metapep_df[metapep_df["Functional Annotation"] == True]
    unique_samples = metapep_df['Sample Name'].drop_duplicates()
        
    return unique_samples.to_list()


@app.callback(
    Output("kegg_pathway_map_link", "href"),
    Output("kegg_pathway_map_link", "disabled"),
    Output("sample_color_table_func_export", "children"),
    Output("sample_color_table_func_export", "style"),
    Input('peptides', 'data'),
    Input("kegg_group_type", "value"),
    Input("kegg_export_samples", "value"),
    Input("pathway_dropdown", "value"),
    Input('tax_barplot_clade_selection_taxa', 'value'),
    Input('tax_barplot_clade_selection_rank', 'value'),
    State('kegg_db_class_data', 'data')
)
def construct_pathway_url(peptide_json, 
                          kegg_group_method,
                          selected_samples,
                          predifined_pathway,
                          filter_clade,
                          clade_rank,
                          kegg_db):
    if predifined_pathway is None or\
        kegg_group_method != "Pathway" or\
        kegg_db is None or\
        peptide_json is None or\
        selected_samples is None:
        return None, True, None, None
    
    # do not allow query construction if more than 4 samples are selected
    if len(selected_samples) > 4:
        return None, True, None, None

    kegg_db = KeggDatabase.read_json(kegg_db)
    
    peptide_df = MetaPepTable.read_json(peptide_json).data
    
    # Keep only taxa in peptide df that are part of selected clade
    if filter_clade and clade_rank and clade_rank != 'Root':
        peptide_df = filter_taxonomy_clade(peptide_df, filter_clade, clade_rank, 'Name')                                              # use to group and count psm's by function
        
    # only keep kegg ko and sample name
    peptide_df = peptide_df[["KEGG_ko", "Sample Name"]]
    peptide_df = peptide_df[peptide_df["Sample Name"].isin(selected_samples)]
    
    # split cells with multiple kegg annotations into separate rows
    peptide_df["KEGG_ko"] = peptide_df["KEGG_ko"].str.split(",") 
    peptide_df = peptide_df.explode("KEGG_ko")
    
    # get all KO's present in pathway, get all ko's in peptide dataset part of pathway
    valid_ko = kegg_db.list_ko(pathway=predifined_pathway)
    peptide_df = peptide_df[peptide_df["KEGG_ko"].isin(valid_ko)]
    peptide_df.drop_duplicates(subset=["KEGG_ko", "Sample Name"], inplace=True)

    # map sample names to color
    cmap = GraphConstants.color_palette
    sample_color_map = {selected_samples[i]: cmap[i] for i in range(len(selected_samples))}
    
    peptide_df["Color"] = peptide_df['Sample Name'].apply(lambda x: sample_color_map[x])
    peptide_df['Color'] = peptide_df['Color'].str.replace("#", f"%{hex(ord('#'))[2:]}")
    
    ko_color_df = peptide_df.pivot(index="KEGG_ko", columns="Sample Name", values="Color")
    
    new_line = '%0A'
    carriage_return = '%0D'
    space = "%20"
    
    kegg_url = GlobalConstants.kegg_map_color_base_url
    
    # specify pathway to display
    kegg_url += f"map={predifined_pathway}&multi_query="
    
    for idx, row in ko_color_df.iterrows():
        colors = row.dropna().to_list()
        colors = space.join(colors)
        kegg_url += str(idx) + space + colors + new_line
        
    # remove trailing new_line
    kegg_url = kegg_url[:-len(new_line)]
    
    # configure table that shows color for each sample
    table_block = sample_color_table_block(sample_color_map)
    return kegg_url, False, table_block, {"margin": "2rem 0rem"}
    