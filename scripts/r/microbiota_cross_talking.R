# file_microbio <- "/Users/jaespinosa/Google Drive/microbiota_elena/FA microbiota/Results microbiota_16.04.20.csv"

# microbio_tbl <- read.csv(file_microbio,
#                          dec=",",
#                          sep=";",
#                          stringsAsFactors = F)
# head (microbio_tbl)

### microbiota data
rel_abundance_by_phylum <- "/Users/jaespinosa/git/food_addiction_analysis/data/microbiota/relative_abundances_by_phylum.csv"
microbiota_by_phylum_ori <- read.csv(rel_abundance_by_phylum,
                         dec=",",
                         # sep=";",
                         sep="\t",
                         check.names = F,
                         stringsAsFactors = F)
head (microbiota_by_phylum_ori)

transpose_df <- function(df) {
    # keep the first column 
    names <-  df[,1]

    # Transpose everything other than the first column
    df.T <- as.data.frame(as.matrix(t(df[,-1])))

    # Assign first column as the column names of the transposed dataframe
    colnames(df.T) <- names
    return(df.T)
}

microbiota_by_phylum_tmp <- transpose_df(microbiota_by_phylum_ori)
write.csv(microbiota_by_phylum_tmp, "/Users/jaespinosa/tmp.csv")
microbiota_by_phylum <- read.csv("/Users/jaespinosa/tmp.csv",
         dec=",",
         check.names = F,
         stringsAsFactors = F)

head(microbiota_by_phylum)


# microbiota_by_phylum$mouse_id <- row.names(microbiota_by_phylum)

## Behavioral data

behavioral_data_path <- "/Users/jaespinosa/git/food_addiction_analysis/data/microbiota/behavioral_data_from_results_microbiota_16_04_20_withoutCategorical.csv"

behavioral_data <- read.csv(behavioral_data_path,
                            dec=",",
                            sep=";",
                            check.names = F,
                            stringsAsFactors = F)
head(behavioral_data)

## merge behavior with microbiota
microbio_behavioral_merged <- merge (microbiota_by_phylum, behavioral_data, by.x= "mouse_id", by.y = "Mice")
head(microbio_behavioral_merged)

cor(microbio_behavioral_merged[,4], microbio_behavioral_merged[,14])

## Hacer todas las correlaciones incluso entre microbiota y behavioral itself
## y luego ya eliminar las que no quiera

m2 <- m %>%
  # convert data to long format
  gather(key="state",value="value",-YEAR,-WEEK) %>%
  # rename columns
  setNames(c("year","week","state","value")) %>%
  # convert year to factor
  mutate(year=factor(year)) %>%
  # convert week to factor
  mutate(week=factor(week)) %>%
  # convert value to numeric (also converts '-' to NA, gives a warning)
  mutate(value=as.numeric(value))













