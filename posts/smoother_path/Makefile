
path.png: lhac-stdout post-processing.sh path-plotter.r
	bash post-processing.sh && Rscript path-plotter.r

lhac-stdout: sim.fasta designMat2 init.mod
	bash run-lhac.sh
