# dedup
This tool can be used to predict duplicates in peaks as signal versus noise from single-end ChIP-seq data.
It requires two input files: alignment file in BAM and peak file. File name should not contain space. Put the two files in the same directory.

The BAM file needs to be coordinates sorted and duplicates should be marked rather than removed (via the "REMOVE_DUPLICATES=false" option in picard MarkDuplicates tool). The peak file is a tab delimited text file without a header line. The first 3 columns are chr, start position and end position. Below are some key information about how to use this tool:

1). Put dedup.sh, lowess.fitting.R and tool_info.txt in the same directory (e.g., within dedup_tool/)
2). Change Line 10 in dedup.sh by providing "/PATH to dedup.sh/", like "/Volumes/Users/James/tools/dedup_tool"
3). Modify the "tool_info.txt" file (provide the path for the tools required by dedup.sh)
4). To get the usage, simply type /PATH to dedup.sh/dedup.sh (replace PATH to dedup.sh with the actual path)
5). To start the analysis, open a terminal, then type something like: /PATH_to_dedup.sh/dedup.sh WORK_DIR=/PATH_to_DIRECTORY_holding_the_BAM_and_PEAK_file BAM=test.bam PEAK=test.encodePeak MAPQ=20 maxDUP=5


MAPQ: minimal mapping quality, reads with mapping quality score below MAPQ will be filtered out
maxDUP: maximum number of duplicates allowed per position/strand. If a position has more than maxDUP duplicates, only maxDUP duplicates will be used. This is based in the notion that certain position have high number of duplicates that more likely represent noise. Set between 3 and 5.

6). output files
*.summary.txt lists the number of non-duplicates, number of duplicates, etc.
*.duplicate_fitted.txt contains the number of non_duplicates (nonredundant reads), total duplicates, duplicates after applying maxDUP filtering (see above), number of duplicates predicted as signal, and number of duplicates predicted as noise in peaks. The sum of the last 2 columns equals to "total duplicates".
*.deduplicated.bam contains alignments for all nonredundant reads in IP (not just from peak regions) and alignments for the duplicates that are predicted as signal within peaks. Theis file is already coordinates sorted and duplicates marked.
