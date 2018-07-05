### K-means Cluster from RPKMs values

#####  Normalize before applying k-means, because the differences in the variances across samples will bias the results
#####  Normalize the expression patterns (normalize w.r.t. the rows).
##### All the RPKM values were normalized (or scaled) such that the RPKMs in each sample had 0 mean and 1 variance. That prevents any sample from having a bigger or smaller weight on the clustering step. After this, clustering was applied using k-means.


###"Genes that were differentially expressed (q<XYZ ) between two or more time points were selected  (7353 genes). The    expression of these genes across the six time points was normalized (0 mean, 1 variance) to avoid that peaks of expression, caused by differences in sample variances, could mask the weight of other samples during clustering. The normalized expression values were then clustered using the Hartigan and Wong K-means clustering algorithm , with 20 clusters, a maximum number of allowed iterations of 25, and 100 random starts. 
