# Pseudobulk and differential expression

## Pseudobulk

A pseudobulk sample is formed by aggregating the expression values from
a group of cells from the same individual. The cells are typically
grouped by clustering or cell type assignment. Individual refers to the
experimental unit of replication (e.g., the individual mice or
patients).

Forming pseudobulk samples is important to perform accurate differential
expression analysis. Cells from the same individual are more similar to
each other than to cells from another individual. This means treating
each cell as an independent sample leads to underestimation of the
variance and misleadingly small p-values. Working on the level of
pseudobulks ensures reliable statistical tests because the samples
correspond to the units of replication.

We can use pseudobulks for example to find the expression changes
between two conditions for one cell type.

## Example

I load a `SingleCellExperiment` object containing gene expression counts
from eight Lupus patient before and after interferon beta stimulation.
The creator of the dataset has already annotated the cell types and if
cell is a singlet.

\
`sce`` ``<-`` ``muscData``::`[`Kang18_8vs8`](https://rdrr.io/pkg/muscData/man/Kang18_8vs8.html)`(``)`\
`# Keep only genes with more than 5 counts`\
`sce`` ``<-`` ``sce``[`[`rowSums`](https://rdrr.io/pkg/MatrixGenerics/man/rowSums.html)`(`[`counts`](https://rdrr.io/pkg/BiocGenerics/man/dge.html)`(``sce``)``)`` ``>`` ``5``, ``]`\
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``sce``)`

The `pseudobulk` functions emulates the `group_by` and `summarize`
pattern popularized by the `tidyverse`. You provide the columns from the
`colData` that you want to use for grouping the data (akin to
`group_by`) and named arguments specifiying how you summarize the
remaining columns (akin to `summarize`). Using the
`aggregation_functions` you can set how the `assay`’s and `reducedDim`’s
are summarized with a named list.

Here, I create a pseudobulk sample for each patient, condition, and cell
type. This means for example that the counts of the 119 B-cells from
patient 101 in the control condition are summed to one column in the
reduced dataset.

The first argument is a `SingleCellExperiment` object. The `group_by`
argument uses [`vars()`](../reference/vars.md) to quote the grouping
columns The `fraction_singlet` and `n_cells` arguments demonstrate how
additional columns from the `colData` are summarized. For
`fraction_singlet`, I use the fact that `mean` automatically coerces a
boolean vector to zeros and ones and `n_cells` demonstrates the
[`n()`](https://dplyr.tidyverse.org/reference/context.html) function
that returns the number of cells that are aggregated for each group.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`glmGamPoi2`](https://github.com/goeva-lab/glmGamPoi)`)`\
`reduced_sce`` ``<-`` `[`pseudobulk`](../reference/pseudobulk.md)`(``sce``,`\
`  group_by ``=`` `[`vars`](../reference/vars.md)`(``ind``, condition ``=`` ``stim``, ``cell``)``,`\
`  fraction_singlet ``=`` `[`mean`](https://rdrr.io/r/base/mean.html)`(``multiplets`` ``==`` ``"singlet"``)``, n_cells ``=`` ``n``(``)`\
`)`\
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``reduced_sce``)`

You can simulate the pseudobulk sample generation and check if you are
using the correct arguments by calling
[`dplyr::group_by`](https://dplyr.tidyverse.org/reference/group_by.html).
Note that the order of the output differs because `group_by`
automatically sorts the keys.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`dplyr`](https://dplyr.tidyverse.org)`, warn.conflicts ``=`` ``FALSE``)`\
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)`(``sce``)`` ``|>`\
`  `[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)`(``)`` ``|>`\
`  `[`group_by`](https://dplyr.tidyverse.org/reference/group_by.html)`(``ind``, condition ``=`` ``stim``, ``cell``)`` ``|>`\
`  `[`summarize`](https://dplyr.tidyverse.org/reference/summarise.html)`(``n_cells ``=`` `[`n`](https://dplyr.tidyverse.org/reference/context.html)`(``)``, .groups ``=`` ``"drop"``)`

With the reduced data, we can conduct differential expression analysis
the same way we would analyze bulk RNA-seq data (using tools like
`DESeq2` and `edgeR`). For example we can find the genes that change
most upon treatment in the B-cells

\
`# Remove NA's`\
`reduced_sce`` ``<-`` ``reduced_sce``[``, ``!`[`is.na`](https://rdrr.io/r/base/NA.html)`(``reduced_sce``$``cell``)``]`\
`# Use DESeq2's size factor calculation procedure`\
`fit`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``reduced_sce``, design ``=`` ``~`` ``condition`` ``*`` ``cell`` ``+`` ``ind``, size_factor ``=`` ``"ratio"``, verbose ``=`` ``TRUE``)`\
`res`` ``<-`` `[`test_de`](../reference/test_de.md)`(``fit``, contrast ``=`` ``cond``(``cell ``=`` ``"B cells"``, condition ``=`` ``"stim"``)`` ``-`` ``cond``(``cell ``=`` ``"B cells"``, condition ``=`` ``"ctrl"``)``)`

A volcano plot gives a quick impression of the overall distribution of
the expression changes.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ggplot2`](https://ggplot2.tidyverse.org)`, warn.conflicts ``=`` ``FALSE``)`\
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)`(``res``, `[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``x ``=`` ``lfc``, y ``=`` ``-`[`log10`](https://rdrr.io/r/base/Log.html)`(``pval``)``)``)`` ``+`\
`  `[`geom_point`](https://ggplot2.tidyverse.org/reference/geom_point.html)`(`[`aes`](https://ggplot2.tidyverse.org/reference/aes.html)`(``color ``=`` ``adj_pval`` ``<`` ``0.01``)``, size ``=`` ``0.5``)`

## Legacy

Originally, `glmGamPoi`’s API encouraged forming pseudobulks after
fitting the model (i.e., within [`test_de()`](../reference/test_de.md)).
The advantage was that this reduced the number of functions. Yet,
internally `glmGamPoi` basically threw away the original fit and re-ran
it on the aggregated data. This meant that computation time was wasted.
Thus the original approach forming the pseudobulk in `test_de` is now
deprecated in favor of first calling
[`pseudobulk()`](../reference/pseudobulk.md) and then proceed by calling
[`glm_gp()`](../reference/glm_gp.md) and
[`test_de()`](../reference/test_de.md) on the aggregated data.

## Session Info

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`
