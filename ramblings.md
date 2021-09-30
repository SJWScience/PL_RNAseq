for QC checks

fastqc outputs a zip file that contains a 'fastqc_data.txt' and 'summary.txt'

fastqc_data.txt contains a line with total sequences and %GC - these can be extracted and averaged across all samples in the experiment and used to identify if there is an outlier (not same species or too few reads)

summary.txt contains pass, fail, warning for each of the standard fastqc metrics

add in an option to be verbose about the QC data and output a QC log file for the user to check. - from the summary.txt of each fastqc output zip files. Add in pre- and post- trimming QC??? not sure yet.
