
# global functions for the faster-report rmarkdown
# input is fastq file, output is stdout and varies
# if true, the saveraw argument saves raw data to text file {function_name}-{filename}.out

faster_gc <- function(x, saveraw = FALSE) {
  require(dplyr)
  require(data.table)

  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  density_obj <- system2("faster2", args = c("--gc", x), stdout = TRUE) %>%
    as.numeric() %>%
    # actually use density() here, not hist(). It returns a density list object with x and y, x is fixed from 1 to 100
    density(from = 0, to = 1, n = 100, na.rm = TRUE) # n is the number of equally spaced points at which the density is to be estimated.

  if(isTRUE(saveraw)) {
    #if(!dir.exists("rawdata")) {dir.create("rawdata")}
    data.table::fwrite(data.frame(x = density_obj$x, y = density_obj$y),
                       file = paste0("rawdata/", "faster_gc", "-", samplename, ".tsv"),
                       sep = "\t")
  }
  # return
  density_obj
}


faster_qscore <- function(x, saveraw = FALSE) {
  require(dplyr)
  require(data.table)
  
  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  density_obj <- system2("faster2", args = c("--qual", x), stdout = TRUE) %>%
    as.numeric() %>%
    # actually use density() here, not hist(). It returns a density list object with x and y, x is fixed from 1 to 50
    density(from = 1, to = 60, n = 60, na.rm = TRUE) # n is the number of equally spaced points at which the density is to be estimated.
  #
  if(isTRUE(saveraw)) {
    #if(!dir.exists("rawdata")) {dir.create("rawdata")}
    data.table::fwrite(data.frame(x = density_obj$x, y = density_obj$y), 
                       file = paste0("rawdata/", "faster_qscore", "-", samplename, ".tsv"), 
                       sep = "\t")
    }
  density_obj
}

faster_len <- function(x, saveraw = FALSE) {
  require(dplyr)
  require(data.table)
  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  
  lens <- system2("faster2", args = c("--len", x), stdout = TRUE) %>% as.numeric()
  lens[lens >= 50000] <- 50000 # all bigger than 50k are set to 50k
    
  density_obj <- lens[lens > 1 & lens <= 50000] %>%
    hist(breaks = seq(1, 50000, length.out = 61), plot = FALSE) # x = breaks or mids, y = counts
    #log10() %>%
    #density(from = 100, to = 60000, n = 60, na.rm = TRUE)
    #density(from = 200, to = 50000, n = 100, na.rm = TRUE)
  if(isTRUE(saveraw)) {
    #if(!dir.exists("rawdata")) {dir.create("rawdata")}
    data.table::fwrite(data.frame(x = density_obj$mids, y = density_obj$counts), 
                       file = paste0("rawdata/", "faster_len", "-", samplename, ".tsv"), 
                       sep = "\t")
  }
  density_obj
}


fastkmers <- function(x, saveraw = FALSE) {
  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  kmers <- system2("fastkmers", args = c("-k 3", "-v", x), stdout = TRUE)
  kmers_tbl <- read.table(text = kmers, sep = "\t", header = TRUE, col.names = c("kmer", "counts")) %>%
    dplyr::arrange(kmer) # in order to compare across files
  if(isTRUE(saveraw)) {
    data.table::fwrite(kmers_tbl, 
                       file = paste0("rawdata/", "fastkmers", "-", samplename, ".tsv"), 
                       sep = "\t"
                       )
  }
  kmers_tbl
}


# duplication level (Illumina only) is calculated with fastkmers -f, using the max read length as kmer 
# returns a data table with occ count percent
duplevel <- function(x, saveraw = FALSE) {
  require(dplyr)
  require(data.table)
  
  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  #faster -l x | sort -n | head -n 1
  faster_out <- system2("faster", args = c("-t", x), stdout = TRUE)
  max_len <- read.table(text = faster_out, header = T)$max_len
  if (max_len > 610) {
    stop(paste("max read length is", max_len, "and you selected Illumina, which can be maximum 610"))
  }
  fastkmers_out <- system2("fastkmers", args = c("-k", as.numeric(max_len), "-v", "-f", x), stdout = T)
  duplevel_tbl <- read.table(text = fastkmers_out, sep = "\t", header = T) %>%
    dplyr::arrange(occ) %>%
    dplyr::mutate(percent = round(count/sum(count), 4)*100) %>%
    head(10)
  
  if(isTRUE(saveraw)) {
    data.table::fwrite(duplevel_tbl, 
                      file = paste0("rawdata/", "duplevel", "-", samplename, ".tsv"), 
                      sep = "\t"
                      )
  }
  duplevel_tbl
}

content_percycle <- function(x, saveraw = FALSE) {
  require(dplyr)
  require(data.table)
  
  samplename <- basename(tools::file_path_sans_ext(x, compression = T))
  
  faster_out <- system2("faster", args = c("-t", x), stdout = TRUE)
  max_len <- read.table(text = faster_out, header = T)$max_len
  if (max_len > 610) {
    stop(paste("max read length is", max_len, "and you selected Illumina, which can be maximum 610"))
  }
  fastkmers_out <- system2("fastkmers", args = c("-k", max_len, "-c", x), stdout = TRUE)
  content_tbl <- read.table(text = fastkmers_out, sep = "\t", header = T) %>%
    dplyr::arrange(cycle)
  
  if(isTRUE(saveraw)) {
    data.table::fwrite(content_tbl, 
                       file = paste0("rawdata/", "content_percycle", "-", samplename, ".tsv"), 
                       sep = "\t"
    )
  }
  content_tbl
}

