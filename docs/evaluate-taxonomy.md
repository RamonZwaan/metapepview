The number of peptides identified from a MS experiment and the accuracy of the taxonomy composition depends greatly on the quality of the protein sequence database used for DB search matching and taxonomy mapping; peptides not present in the database will not be matched and species not present in the database will not be represented in the community composition. A potential indicator for a bad protein database is a low number of DB search matches, especially when combined with many confident *de novo* peptides. However, the number of DB search matches considered acceptable depends strongly on the type of sample. In addition, other experimental factors may influence the number of DB search matches.

One method to evaluate the quality/coverage of the protein sequence database used for DB search matching and taxonomy annotation is by comparing taxonomy composition derived from the protein database, to the taxonomy composition obtained from [Unipept](https://unipept.ugent.be/) annotation of *de novo* peptide sequences. Using this method, we can evaluate absence of proteins/taxonomies from the protein sequence database. 

## The *Community composition evaluation* module

MetaPepView provides a separate module for this evaluation: *Community composition evaluation*. Its interface is similar to the *Community composition* module interface, except here, only one sample is visualized at a time. To visualize a sample in this module, both DB search and *de novo* metaproteomics data has to be supplied for the sample, and taxonomy annotation from the local accession-taxonomy map, as well as [Unipept](https://unipept.ugent.be/) (check *Annotate peptides to Unipept* in *Taxonomy annotation* options) has to be performed.

![[Taxonomy_db_search_de_novo_module.PNG]]

On the right side, the top chart shows the taxonomy compositions for a single sample: The left bar shows the composition from DB search matches, annotated with the user provided taxonomy map. The right bar shows the taxonomy composition as annotated by Unipept. This may be from DB search + *de novo* peptides or from *de novo* peptides only (see **Unipept composition de novo only** setting). The bottom chart (Abundance ratio) shows over- and underrepresentation of taxonomy groups between the taxonomy compositions.

The abundance ratio chart calculates for each taxonomy group the "hyperbolic tangent of the log ratio": 
$$
abundance\ ratio=tanh(log(x_{de\ novo}) - log(x_{db\ search}))
$$ 

Here, $x$ represents the fraction of the taxonomy group in either the DB search or *de novo* composition. The ratio is capped between a range of $[-1, 1]$, where -1 means a taxonomy group is only observed in DB search composition, +1 means the taxonomy group is only observed in *de novo* composition, and 0 means that the taxonomy group has identical relative abundance between the compositions. Thus, a higher number implies a larger fraction of a taxonomy group in the *de novo* composition relative to the DB search composition.

Similarly to the *Community composition* module, there are several options to modify the visualization:

- **Select sample**: Select the sample of the project to display (only valid samples that have all data necessary to visualize are shown).
- **Select taxa**: Select specific taxonomies to display (Only if *Taxa display type* is set to *Custom taxa*)
- **Taxa display type**: Specify what taxonomies to display: *Top 10/23* abundant taxonomies, or user specified *Custom taxa* (Taxa provided in *Select taxa*)
- **Display rank**: Select taxonomy level to display.
- **De novo only global annotation**: For the *de novo* (Unipept) composition, filter out any peptide sequences that were also matched with DB search from the abundance calculation.
- **Quantification**: Quantify taxonomy groups as summed PSM counts or as summed signal intensities.
- **Normalize abundances**: Check to display taxonomies as fraction of total composition per sample.
- **Include unannotated**: Add peptides without annotation (at specified taxonomy rank) as separate "Undefined" category in the composition.
- **Filter by clade**: Select a specific taxonomy group to visualize; Specify a *Root taxonomy* group belonging to the selected *Clade rank* level, and only "offspring" taxonomy groups belonging to the root taxonomy will be visualized. Note: *Display rank* should be at a lower level than the *Clade rank*.
