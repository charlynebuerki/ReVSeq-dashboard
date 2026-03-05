#read in hq data for dashboard data, fetch pileups and assemble them by virus
library(dplyr)
library(stringr)
library(readr) # readr is faster for large TSVs

#read in metadata needed for viruses
virus <- "flu-b"
metadata_path <- sprintf("/Users/cbuerk/Documents/Proj-ReVSeq/ReVSeq-dashboard/Nextstrain-pipelines/data/local_data/%s/metadata.tsv", virus)
metadata <- read_tsv(metadata_path)

# get plates
plates_data <-read.csv("Users/cbuerk/Documents/Proj-ReVSeq/swiss_co-circulating_viruses/Data/data/hq_data.csv") %>% 
  select(pseudonymized_id, barcode) %>% 
  filter(pseudonymized_id %in% (metadata %>% pull(SampleID))) %>% 
  mutate(plate = str_sub(barcode, 1, -4)) %>% 
  select(pseudonymized_id, plate) 

metadata_all <- left_join(metadata, plates_data, by=join_by(SampleID == pseudonymized_id)) %>% 
  filter(!(is.na(plate)))

#now should be able to gather the filepath for all samples to pileup data
metadata_all <-metadata_all %>% 
  mutate(
    id = str_sub(SampleID,4),
    file_path = sprintf("/Users/cbuerk/Documents/Proj-ReVSeq/swiss_co-circulating_viruses/Data/data/depth_files/%s/%s/depth/%s_depth.tsv",
                        plate, id, id )
    )

# from filepaths, direct to folder for all viruses to extract the depth tsv files 
dest_folder <- sprintf("/Users/cbuerk/Documents/Proj-ReVSeq/ReVSeq-dashboard/Extract-pileup/z_get_pileups_from_local_revseq/%s", virus)

#Ensure the directory exists
if (!dir.exists(dest_folder)) {
  dir.create(dest_folder, recursive = TRUE)
}

# file.copy returns a vector of TRUE/FALSE indicating success for each file
success <- file.copy(from = metadata_all$file_path, to = dest_folder, overwrite = TRUE)

# Optional: Check if any failed
if(any(!success)) {
  warning("Some files failed to copy.")
}

# now that it's copied, open files and filter for only matching reference fasta name on pileup

#ad hoc function
#given a virus reference, read from the lookup table and translate to substrain name
get_virus_reference_conversion_table<-function()
{
  virus_lookup<-read.table("/Users/cbuerk/Documents/Proj-ReVSeq/swiss_co-circulating_viruses/Data/resources/virus_lookup_table.bed", sep="\t", header = FALSE)
  names(virus_lookup) <- c("accession_number", "start", "end", "substrain_name")
  
  return(virus_lookup)
}

virus_to_lookup <- metadata_all$virus_identified %>% unique()

#get the accession number associated to virus identified
lookup.tbl<-get_virus_reference_conversion_table() %>% 
  filter(grepl(virus_to_lookup, substrain_name))

#if lookup.tbl is dim = 0, try for perfect matching:
if(dim(lookup.tbl)[1]==0) {
  lookup.tbl<-get_virus_reference_conversion_table() %>% 
    filter(substrain_name == virus_to_lookup)
}


#filter the pileup file one by one by accession number, keep only those from lookup tbl (contained in frst column of read dataframe)

#get accession to filter on:
target_accessions <- lookup.tbl %>% pull(accession_number)
#files to process:
files_to_process <- list.files(dest_folder, full.names = TRUE)

message(sprintf("Starting filtering for %s samples...", length(files_to_process)))

#rewrite to folder modified pileup 
# 3. Iterate through each file, filter, and overwrite


purrr::walk(files_to_process, function(current_path) {
  
  # Read the pileup/depth file
  # Note: col_names = FALSE if your depth files don't have headers
  # Typical samtools depth format: [Reference, Position, Depth]
  temp_depth <- read_tsv(current_path, col_names = FALSE, show_col_types = FALSE)
  
  # Filter: Keep rows where the first column matches our target accession(s)
  # X1 is the default name given by readr for the first column if no header exists
  filtered_depth <- temp_depth %>%
    filter(X1 %in% target_accessions)
  
  # Check if filtering resulted in data
  if(nrow(filtered_depth) > 0) {
    # Overwrite the file with the filtered version
    write_tsv(filtered_depth, current_path, col_names = FALSE)
  } else {
    warning(sprintf("No matching accessions found in: %s", basename(current_path)))
  }
})

message("Filtering complete. Files are ready for the dashboard.")
