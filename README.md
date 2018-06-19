### K-means Cluster from RPKMs values

#####  Normalize before applying k-means, because the differences in the variances across samples will bias the results
#####  Normalize the expression patterns (normalize w.r.t. the rows).
##### All the RPKM values were normalized (or scaled) such that the RPKMs in each sample had 0 mean and 1 variance. That prevents any sample from having a bigger or smaller weight on the clustering step. After this, clustering was applied using k-means.
