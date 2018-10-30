# Hyline-Genome


# giving permission to Plnik 
	chmod +x plink


# Reading Data
	head -n 110  d2.ped | cut -f1-10

# Saving the IDs in text file
	cut  -f1 d2.ped > newFile && mv newFile file


   
# Convert ped to bed
	./plink --no-fid --no-parents --no-sex --no-pheno --chr-set 95 --allow-extra-chr --file d1  --make-bed --out d1
	./plink --no-fid --no-parents --no-sex --no-pheno --chr-set 95 --allow-extra-chr --file d2  --make-bed --out d2



# Merging the two bed
	./plink --chr-set 95 --allow-extra-chr --bfile d1 --bmerge d2.bed d2.bim d2.fam --make-bed --out data

# Removing unknown location and sex chromosomes
    	./plink --chr-set 95 --bfile data --allow-extra-chr --chr 1-28 --make-bed --out d


# Sample Quality Control

# Data Missingness
	./plink --dog --bfile d --missing --out out/d

# Heterozygosity rate
	./plink --dog --bfile d --het --out out/d

# Related and duplicate individuals

	./plink --dog --bfile d --extract raw-GWA-data.prune.in --genome --min 0.2 --out pihat_min0.2

	awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome

	Rscript --no-save Relatedness.R

	cat fail-* | sort -k1 | uniq > fail-qc-inds.txt

	./plink --dog --bfile d --remove fail-qc-inds.txt --make-bed --out ci_d


# Markers Quality Control

	# Delete SNPs with missingness >0.2.
		./plink --dog --bfile ci_d --geno 0.1 --make-bed --out x_d

 
	#Generate a plot of the MAF distribution.
		./plink --dog --bfile x_d --freq --out MAF_check
		Rscript --no-save MAF_check.R


 	#Delete SNPs with MAF <0.21.
       ./plink --dog --bfile x_d --maf 0.021 --make-bed --out x2_d


	# Check the distribution of HWE p-values of all SNPs.

		./plink --dog --bfile x2_d --hardy
		awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
		Rscript --no-save hwe.R



	# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).

		./plink --dog --bfile x2_d --hwe 1e-6 --make-bed --out x3_d
	

###################################

	#Update FID
		./plink --dog --bfile cd --update-ids R/ufi.txt --make-bed --out cd_x
	#Update sex
		./plink --dog --bfile cd_x --update-sex R/usex.txt --make-bed --out cd_x2


###################################
############## GWAS ###############
###################################

		# Making a GRM
			./gcta64  --bfile HL --autosome --autosome-num 28 --make-grm  --out HL

		# Principal component analysis
			./gcta64  --grm HL --autosome --autosome-num 28 --pca 3 --out HL

		# GCTA-GREML: Estimate variance explained by all the SNPs

			# Model one y(1:3) = mean + fixed effects (Dam+ sex + Chamber) + GRM + e 
			for i in {1..3}
			do
			./gcta64 --reml  --grm HL --autosome --autosome-num 28 --reml-pred-rand --mpheno $i --pheno R/phe.txt --covar R/cov2.txt --out HL$i
			done

			# Model two y4 = mean + fixed effects (Dam+ sex + Chamber + 3 PCs) + GRM + e 
			./gcta64 --reml  --grm HL --autosome --autosome-num 28 --reml-pred-rand --mpheno 4 --pheno R/phe.txt --covar R/cov2.txt --qcovar HL.eigenvec --out HL4


		# Single SNP genome wide association SSGWA

			# Model one y(1:3) = mean + fixed effects (Dam+ sex + Chamber) + GRM + e 
			for i in {1..3}
			do
   			./gcta64 --bfile HL --grm HL --autosome --autosome-num 28 --mlma-loco --mpheno $i --pheno R/phe.txt --covar R/cov2.txt --out HL$i
			done

			# Model two y4 = mean + fixed effects (Dam+ sex + Chamber + 3 PCs) + GRM + e 
			./gcta64 --bfile HL --grm HL --autosome --autosome-num 28 --mlma-loco --mpheno 4 --pheno R/phe.txt --covar R/cov2.txt --qcovar HL.eigenvec --out HL4












