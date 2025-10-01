# Community function visualization

MetaPepView provides the *community functions* module for visualization of protein functions. This module is available for all samples for which functional annotation data is supplemented from EggNOG or KEGG format. The module quantifies [KEGG](https://www.kegg.jp/) terms from multiple samples for comparative functional analysis. In addition, the KEGG Pathway and BRITE database is integrated in the dashboard for browsing and visualization of various types of functional groups.

![[Function_module.PNG]]

The quantification of KEGG terms between samples is visualized in a grouped barplot. All samples in the project for which a functional annotation source was provided during import are shown. When opening the module, initially, the top 20 most abundant (averaged over all samples) KEGG KO terms, quantified by PSM counts, are displayed. Left from the graph, a range of options is provided to visualize specific functional groups, display different grouping types, and to scale the quantifications:

- **KEGG ID display format**: Format of the functional groups to display per sample (categories on x-axis).
- **Select groups**: For displaying only a specific functional group, filter by:
	- *Pathway*: display KEGG terms belonging to certain pathway, selection from *Select pathway* and/or *Select module* menu's
	- *BRITE*: display KEGG terms belonging to BRITE hierarchy group, selection from *Select BRITE group* menu.
	- *Manual*: Manually select KEGG terms to display.
- **Quantification**: Quantify function groups as summed PSM counts or as summed signal intensities.
- **Normalize sample**: Select sample as normalization factor: For all samples, function expression is quantified as over-/underrepresentation factor for the selected sample.
- **Combine multiple annotations**: For peptides with multiple annotations, if checked, combine annotation names into a separate category. If unchecked, count peptide abundance towards each functional element separately. For example: `Pept_A` has KEGG KO annotation: `K00370,K07306,K17050`. When the option is checked, `K00370,K07306,K17050` becomes a separate function group (separate from `K00370` or any other separate element). If unchecked, the abundance of `Pept_A` is counted towards all of `K00370`, `K07306`, `K17050`.

!!! note
	When displaying KO terms, peptides may only have multiple annotations if multiple proteins were matched. When EC or Module terms are displayed, a single KO may match to multiple EC or Module terms. After translation of KO term to EC or Module term('s), MetaPepView will combine them in the same method when *Combine multiple annotations* is checked. However, in cases of peptides with multiple KO terms (due to multiple protein matches), MetaPepView will ensure that a peptide sequence will not be counted double within a single EC or Module term.
	Therefore, while redundancy is present when *Combine multiple annotations* is unchecked due to potential counting a single peptide towards multiple groups, peptides will only be counted once within a single term.

- **Include taxonomies**: Supplement functional annotation chart with taxonomy information. This shows which taxa contribute to the functional term expression. Taxonomy names can be displayed by hovering over the bars in the chart.
- **Fractional abundances**: Quantify functional terms as fraction of the total expressed functional terms.
- **Filter by clade**: Filter expressed functional terms to only those from a specified taxonomy group; Specify a *Root taxonomy* group belonging to the selected *Clade rank* level, and only terms from that taxonomy are displayed.
- **Visualize pathway map**: Export KEGG KO terms to a pathway map provided by KEGG. Select samples to visualize in the map in *Select Samples* and export the data to the map with the *Show pathway map* button. This directs you to the KEGG website. Thus, an internet connection is required.
