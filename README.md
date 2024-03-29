# CS 1640: Bioinformatics Software Design

A prior way to explore the statistical dependencies between genes in single-cell gene expression
datasets is to use multivariate information (MVI). Single-cell datasets are large enough to allow
accurate estimations of probability distributions between more than two variables, allowing a
useful measure of MVI in combination with partial information decomposition (PID).
Algorithms incorporating PID allow for analysis of the statistical relationships between triplets
of variables, creating undirected networks. These networks include the functional interactions
between the genes. An inference algorithm based on the MVI measure of PID and in silico
analyses provide a strong basis for analyzing single-cell data. By utilizing in silico methods,
network inference approaches can be quantitatively assessed. This approach can be preferable to
in vitro methods of analyses as the “true” gene regulatory networks is known and can be
compared to. PIDC, a PID-based algorithm and context, identifies functional relationships
between genes more efficiently and effectively than other networking algorithms. However,
there are some limitations to this proposed method. One limitation includes relationships
between genes only being detected when there is sufficient variability in gene expression.
Additionally, co-regulatory relationships are likely to be detected where genes under the
influence the same regulator. Without any further assumptions, no causal relationships can be
determined. 

Assessing the accuracy
To show that the proposed method, the gene regulatory networks need to be identified with
reasonable degree of accuracy that is better than other methods. To do this, the ground truth of
the gene regulatory network must be known. However, the ground truth is often not known. One
solution to not knowing the ground truth is to utilize simulated data networks or curated
networks where the interactions are well known. These well-known interactions for curated
networks can act as the ground truths.
