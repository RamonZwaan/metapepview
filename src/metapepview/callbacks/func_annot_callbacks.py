from dash import Dash, dash_table, html, dcc, callback, Output, Input, State, ctx
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

from metapepview.server import app

# import layout elements
from metapepview.html_templates import hidden_graph_with_text, \
    sample_color_table_block

from metapepview.backend import *
from metapepview.backend.utils import truncate_end
from metapepview.backend.plots import pathway_abundance_barplot
from metapepview.constants import GlobalConstants
from metapepview.backend.utils.functional_plot_utils import *

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict


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

    # check if everything is present to display plot
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
    if kegg_group_method == "Manual" and (custom_prot is None or custom_prot == []):
        block_element = hidden_graph_with_text("pathway_barplot_figure",
                                                "Select protein terms in filter settings...")
        return block_element, dict(), "Figure"
    

    peptide_df = metapep_obj.data

    # divide for each sample the psm value by the sum of psm for that sample, if specified
    if fractional_abundance is True:
        peptide_df = calculate_frac_abundance(peptide_df, ycol)
    
    # if taxa selection made on the taxonomy barplot, filter protein abundances by taxa
    # get taxonomy id from barplot selection point
    if filter_clade and clade_rank and clade_rank != 'Root':
        peptide_df = filter_taxonomy_clade(peptide_df, filter_clade, clade_rank, 'Name')

    peptide_df = process_kegg_groups(
        peptide_df=peptide_df,
        kegg_db=kegg_db,
        ycol=ycol,
        include_taxa=include_taxa,
        kegg_group_method=kegg_group_method,
        kegg_group_format=kegg_format,
        custom_prot=custom_prot,
        predifined_pathway=predifined_pathway,
        predifined_module=predifined_module,
        predifined_brite_group=predifined_brite_group,
        combine_annotations=combine_annotations
    )

    # keep only top n ko's based on largest contribution across samples
    func_abundances = peptide_df.groupby(by="Protein Name")[ycol].agg('sum')
    peptide_df = get_top_functions(peptide_df, func_abundances)

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
    
    # define column lists required during processing when taxonomy is included and when not 
    if include_taxa is True:
        tax_col = "Taxonomy Name"
    else:
        tax_col = None

    plot_title = configure_plot_title(kegg_group_format=kegg_format,
                                      filter_clade=filter_clade,
                                      clade_rank=clade_rank,
                                      func_abundances=func_abundances)    
    
    plot = pathway_abundance_barplot(peptide_df,
                                    "Sample Name",
                                    ycol,
                                    "Protein Name",
                                    tax_col=tax_col,
                                    custom_title=plot_title,
                                    custom_xname=kegg_format,
                                    custom_yname=ytitle)
    
    plot.update_layout(title="")

    return dcc.Graph(figure=plot,
                     id="pathway_barplot_figure",
                     style={"height": "40rem"}), dict(), plot_title


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
    if kegg_db is None:
        table_block = [html.P(
            "Import KEGG dataset (see sidebar)",
            className="fst-italic ms-5")]
        return None, True, table_block, None

    if peptide_json is None:
        table_block = [html.P(
            "Import or create project with functional annotation data",
            className="fst-italic ms-5")]
        return None, True, table_block, None

    if predifined_pathway is None or\
        kegg_group_method != "Pathway":
        table_block = [html.P(
            "Select KEGG pathway map...",
            className="fst-italic ms-5")]
        return None, True, table_block, None
    
    if selected_samples is None or len(selected_samples) == 0:
        table_block = [html.P(
            "Select sample in dropdown menu...",
            className="fst-italic ms-5")]
        return None, True, table_block, None
    
    # do not allow query construction if more than 4 samples are selected
    if len(selected_samples) > 4 :
        table_block = [html.P(
            "Too many samples selected...",
            className="fst-italic ms-5")]
        return None, True, table_block, None

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
    kegg_url = kegg_url.rstrip(new_line)

    # configure table that shows color for each sample
    table_block = sample_color_table_block(sample_color_map)
    # check length of URL, it should not cross the element limit
    if len(kegg_url.split(new_line)) > GlobalConstants.kegg_url_max_elems:
        table_block = [html.P(
            "Too many KO elements to export to KEGG",
            className="fst-italic ms-5")]
        return None, True, table_block, None

    return kegg_url, False, table_block, {"margin": "2rem 0rem"}


@app.callback(
    Output('export_functions_button', 'disabled'),
    Input('peptides', 'data')
)
def toggle_export_button(peptide_json):
    if peptide_json is None:
        return True
    else:
        return False


@app.callback(
    Output('download_community_functions_csv', 'data'),
    Input('export_functions_button', 'n_clicks'),
    State('peptides', 'data'),
    State("kegg_group_type", "value"),
    State("kegg_display_format_radio", "value"),
    State("brite_group_dropdown", "value"),
    State("pathway_dropdown", "value"),
    State("kegg_module_dropdown", "value"),
    State("custom_pathway_items", "value"),
    State("barplot_pathway_fraction_checkbox", "value"),    
    State("barplot_pathway_include_taxa_checkbox", "value"),
    State('barplot_func_quantification_column', 'value'),
    State("func_annot_combine_duplicates", "value"),
    State('tax_barplot_clade_selection_taxa', 'value'),
    State('tax_barplot_clade_selection_rank', 'value'),
    State('kegg_db_class_data', 'data'),
    prevent_initial_call=True
)
def export_func_composition(button_click, 
                            peptide_json, 
                            kegg_group_method,
                            kegg_format,
                            predifined_brite_group,
                            predifined_pathway,
                            predifined_module,
                            custom_prot,
                            fractional_abundance,
                            include_taxa,
                            quant_method,
                            combine_annotations,
                            filter_clade,
                            clade_rank,
                            kegg_db):
    # set y axis column based on quantification method
    if quant_method == "Match Count":
        ycol = "PSM Count"
    else:
        ycol = "Area"

    # check if everything is present to display plot
    if kegg_db is not None:
        kegg_db = KeggDatabase.read_json(kegg_db)
    else:
        raise PreventUpdate
    if peptide_json is None:
        raise PreventUpdate
    metapep_obj = MetaPepTable.read_json(peptide_json)
    if metapep_obj.functional_annotation_present is False:
        raise PreventUpdate
    if kegg_group_method == "Manual" and (custom_prot is None or custom_prot == []):
        raise PreventUpdate
    

    peptide_df = metapep_obj.data

    # divide for each sample the psm value by the sum of psm for that sample, if specified
    if fractional_abundance is True:
        peptide_df = calculate_frac_abundance(peptide_df, ycol)
    
    # if taxa selection made on the taxonomy barplot, filter protein abundances by taxa
    # get taxonomy id from barplot selection point
    if filter_clade and clade_rank and clade_rank != 'Root':
        peptide_df = filter_taxonomy_clade(peptide_df, filter_clade, clade_rank, 'Name')

    peptide_df = process_kegg_groups(
        peptide_df=peptide_df,
        kegg_db=kegg_db,
        ycol=ycol,
        include_taxa=include_taxa,
        kegg_group_method=kegg_group_method,
        kegg_group_format=kegg_format,
        custom_prot=custom_prot,
        predifined_pathway=predifined_pathway,
        predifined_module=predifined_module,
        predifined_brite_group=predifined_brite_group,
        combine_annotations=combine_annotations
    )

    # Order all functions as follows: 
    #   First order by sample name
    #   Then order by function name from most abundant to least abundant.
    #   For each function name, keep all rows together.
    ord_df = []
    for sample_name in peptide_df["Sample Name"].unique():
        sample_pept_df = peptide_df[peptide_df["Sample Name"] == sample_name]
        sample_func_abundances = (
            sample_pept_df.groupby(by=["Protein Name"])[ycol]
            .agg('sum')
            .sort_values(ascending=False)
        )
        sample_pept_df = (sample_pept_df
            .set_index("Protein Name")
            .sort_values("PSM Count", ascending=False)
            .loc[sample_func_abundances.index.to_list(), :]
            .reset_index()
        )
        ord_df.append(sample_pept_df)
    
    ord_df = pd.concat(ord_df).reset_index(drop=True)

    # Give "Protein Name" a more suitable name
    # filter dataset with annotations belonging to nitrogen cycle
    ord_df = ord_df.rename(columns={"Protein Name": kegg_format})

    return dcc.send_data_frame(ord_df.to_csv, "community_functions.tsv", sep="\t")
