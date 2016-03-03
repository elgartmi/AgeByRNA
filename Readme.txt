These files are matlab scripts that were used to compute the data in Figure 2 in the 
Nature Communications article "Impact of gut microbiota on the fly's germline" by Elgart et al.
If you use this code, please cite both this repository and the article.

The code estimates the developmental stage based on the transcriptional profile (while using 
data from previously published timeseries measurements)
The developmental stage of Drosophila melanogaster embryos at 2hr AED (with and without bacterial removal the in preceding generation) 
was estimated by comparing the RNA-seq mRNA data to a reference time-series transcriptional data from Lott et al. 
The estimation was based on the method developed by Efroni et al. with a few modifications. The estimation requires 
identification of a subset of genes with a single peak of expression which is preceded and followed by a monotonic 
change along the reference time-series. This subset is used as a ruler for estimating the developmental stage of a 
query sample by determining, for each gene, its location along the time-series based on its measured expression in 
the sample. The stage of the sample as a whole is then determined by averaging the time estimations for all the genes 
in the subset. We adjusted the procedure to allow work with RNA-seq data and used two different schemes for estimating 
the stage: The first employed a stringent reference set of monotonically decreasing maternal genes (attached to the project in MS-Excel file) 
used as in Efroni et al. This reference was defined by selecting from Lott et al. all the genes which satisfy the 
criteria in Efroni et al., and narrowing this list by intersecting it with curated experimental maternal data from 
Thomsen et al. The second scheme was based on representing the time course data of Lott et al.12  as a set of time 
varying vectors, each comprising averaged expression data for a particular time point. For every vector of average 
transcript levels in the sample, we computed the Euclidian distance from all the reference vectors in the time series 
data. The stage of each sample was then estimated by identifying the two reference vectors which correspond to the lowest 
distance from the sample and interpolating the stage between them based on the relative distance from them.
The estimations by the two schemes yielded very similar stages (no more than half a cycle difference). In this paper we used 
the estimation based on the first scheme and verified the main findings using the second scheme.  

The entry point to the code is the function "buildmodels.m" which accepts two parameters:
	skip_timepoint	-	this is for testing purposes: if you have a timeseries, then you first want to check if the code
						runs proper - i.e. you can "hide" one of the timepoins and use it as the data. It supposed to fall
						in between other points. If not, use other parameters (see below).
	matzygall		-	this signal what kind of genes to use as the ruler:
						'zyg' - zygotic (rising)
						'mat' - maternal (falling)
						'spk' - spiking (rising then at some point falling)
						'all' - uses all genes for simple euclidean distance estimation (this is actually quite robust 
								and can be used as a sanity check)