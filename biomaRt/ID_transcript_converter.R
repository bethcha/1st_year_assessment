# documentation from https://asia.ensembl.org/info/data/biomart/biomart_r_package.html

# load packages 
library(biomaRt)

# clear biomart cache
biomartCacheClear()

# initiate biomart ensembl - archived version for GRCz10
listEnsembl(version=80)
ensembl <- useEnsembl(biomart = 'genes')

# using Ensembl Biomart version 80 (which has GRCz10)
ensembl_80 <- useEnsembl(biomart='genes', dataset='drerio_gene_ensembl', 
                         version = 80)

# Ensembl - filter the search
# attributes = vector of attributes you want to retrieve (output of the query)
# filters = vector of filters used as an input to the query
# values = vector of values for the filters

# check available filters and attributes for biomart version 80
listAttributes(ensembl_80)
listFilters(ensembl_80)

# attributes needed: chromosome_name, start_position, end_position, strand


################################################################################


# load in the files needed. 
# These can be found in the '../Plotting rhythmic transcripts/input' folder 
ldwt_all_rhythmic <- read.csv("trueLDWT_all_genes", sep=' ')
sdwt_all_rhythmic <- read.csv("trueSDWT_all_genes", sep=' ')
ldclock_all_rhythmic <-read.csv("trueLDCLOCK_all_genes", sep=' ')
sdclock_all_rhythmic <-read.csv("trueSDCLOCK_all_genes", sep=' ')

# order columns from smallest to largest, just good practice
ldwt_rhythmic_sorted <- ldwt_all_rhythmic[order(ldwt_all_rhythmic$name),]
sdwt_rhythmic_sorted <- sdwt_all_rhythmic[order(sdwt_all_rhythmic$name),]
ldclock_rhythmic_sorted <- ldclock_all_rhythmic[order(ldclock_all_rhythmic$name),]
sdclock_rhythmic_sorted <- sdclock_all_rhythmic[order(sdclock_all_rhythmic$name),]


# load in the columns needed 
ldwt_rhythmic = ldwt_rhythmic_sorted$name
sdwt_rhythmic = sdwt_rhythmic_sorted$name
ldclock_rhythmic = ldclock_rhythmic_sorted$name
sdclock_rhythmic = sdclock_rhythmic_sorted$name


# remove the version number from the transcript ID
ldwt_cleaned = sub("*\\.[0-9]", "", ldwt_rhythmic)
sdwt_cleaned = sub("*\\.[0-9]", "", sdwt_rhythmic)
ldclock_cleaned = sub("*\\.[0-9]", "", ldclock_rhythmic)
sdclock_cleaned = sub("*\\.[0-9]", "", sdclock_rhythmic)


################################################################################


# run the conversion for GrCz10

# filter - input
bmart_filters = "ensembl_transcript_id"
# attributes - output
bmart_attributes = c("ensembl_transcript_id", "external_transcript_name", "ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand")

# run for LDWT
# run the Ensembl Biomart function (calls to server)
ldwt_output <- getBM(attributes = bmart_attributes, filters = bmart_filters, 
                     values = ldwt_cleaned, mart = ensembl_80)
ldwt_output

# run for SDWT
sdwt_output <- getBM(attributes = bmart_attributes, filters = bmart_filters, 
                     values = sdwt_cleaned, mart = ensembl_80)
sdwt_output

# run for LDCLOCK
ldclock_output <- getBM(attributes = bmart_attributes, filters = bmart_filters, 
                     values = ldclock_cleaned, mart = ensembl_80)
ldclock_output

# run for SDCLOCK
sdclock_output <- getBM(attributes = bmart_attributes, filters = bmart_filters, 
                     values = sdclock_cleaned, mart = ensembl_80)
sdclock_output


################################################################################


# save to file:
# these files are already available in the 'Plotting rhythmic transcripts/input' folder

write.csv(ldwt_output, "LDWT_all_rhythmic_chr_GRCz10.csv", row.names = FALSE,
          quote=FALSE)
write.csv(sdwt_output, "SDWT_all_rhythmic_chr_GRCz10.csv", row.names = FALSE,
          quote=FALSE)
write.csv(ldclock_output, "LDCLOCK_all_rhythmic_chr_GRCz10.csv", row.names = FALSE,
          quote=FALSE)
write.csv(sdclock_output, "SDCLOCK_all_rhythmic_chr_GRCz10.csv", row.names = FALSE,
          quote=FALSE)
