
####  Merge_lines function  ####

# 1. load individual line profiles
# 2. normalize X-axis to Actin peaks and Y-axes to Baseline & Max intensity
#    > merge to long-format Table with Cell & ROI identifier
# 3. Define regions on line (Membrane-Peaks, Center, Minimum in Cell)
# 4. Calculate summary stats for each line profile (Mean Membrane, Mean Center, Mean Min, Mem/Center, Mem/Min)
#    > merge to Table with Cell & ROI identifier
# 5. Create quality control plot with Actin > Membrane detection, PH + summary lines


merge_lines <- function(folderdir) {
  
  #### Global definitions
  # Create Folders
  imagedir <- file.path(folderdir, "Quality_Control")
  dir.create(imagedir, showWarnings = FALSE)
  resultsdir <- file.path(folderdir, "Results")
  dir.create(resultsdir, showWarnings = FALSE)
  
  
  ## ADD FILESIZE ALREADY HERE (Memory efficiency)
  
  
  # Define global dataframes
  rawData <- tibble()
  
  Stats <- tibble(
    "name" = character(), 
    "ID" = character(), 
    "MeanMembrane_Actin" = numeric(), 
    "MeanCenter_Actin" = numeric(), 
    "MemCentRatio_Actin" = numeric(), 
    "MeanMembrane_PH" = numeric(), 
    "MeanCenter_PH" = numeric(), 
    "MemCentRatio_PH" = numeric(), 
    "MeanCellMin" = numeric(), 
    "MemMinRatio" = numeric(),
    "Responding" = factor()
  )
  
  
  # Define folder of line profiles
  folderlist <- list.files(path = folderdir, include.dirs = F, recursive = F, pattern = ".txt$") # pattern only includes txt-files
  
  
  
  ## find peaks
  ## developed by user "stats g": https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
  
  find_peaks <- function(x, m = 3) {
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i) {
      z <- i - m + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + m + 1
      w <- ifelse(w < length(x), w, length(x))
      if (all(x[c(z:i, (i + 2):w)] <= x[i + 1])) {
        return(i + 1)
      } else {
        return(numeric(0))
      }
    })
    pks <- unlist(pks)
    pks
  }
  
  
  
  
  ##### loop through each line profile individually
  for (i in 1:length(folderlist)) {
    
    
    ### 1. load profile
    temp <- read.table(file.path(folderdir, folderlist[i]), header = TRUE, row.names = 1, sep = "\t")
    
    # Add Image & cell identifier (ROI-name) to each profile and to Summary Stats table
    Stats[i, 1] <- folderlist[i]
    Stats[i, 2] <- folderlist[i]  # second time to keep individual ID (in case that two ROIs of differnt images have the same name)
    temp$name <- folderlist[i]
    temp$ID <- folderlist[i]
    
    # split by underscore "_"
    temp <- separate(data = temp, col = name, into = c("Image", "ROI"), sep = ".tif_")
    
    
    
    ### 2. Normalize Axes
    # normalize X-axes
    peak_ID <- find_peaks(temp$Actin, m = 4)                 # detect many local maxima with m=4
    Maxima <- temp[peak_ID, ]                                # subset line profile by local Maxima IDs
    Maxima <- cbind(Maxima, peak_ID)                         # required for subsetting temp by IDs for region definition
    
    
    Maxima <- Maxima %>% arrange(desc(Actin))                # sort for two global maxima
    Maxima <- Maxima[1:2, ]                                  # remove remaining maxima
    Maxima <- Maxima %>% arrange(X)                          # sort by X-axis again to start with first peak (otherwise aligning to -1 or +1)
    
    # rescale X-Axis
    X_dist <- abs(Maxima[1, 1] - Maxima[2, 1])               # distance between peaks in ?m
    temp$X_norm <- ((temp$X - Maxima[1, 1]) / X_dist) + 1    # (X-X(peak1))/X-dist > peak 1 = 1, peak 2 = 2
    
    # rescale Actin-Intensity    >> this should be a function
    Act_Min <- min(temp$Actin)                               # minimal Intensity
    Act_dist <- max(Maxima$Actin) - Act_Min                  # dynamic range of intensity: maximal intensity - min
    temp$Actin_norm <- (temp$Actin - Act_Min) / Act_dist     # normalized intensity = (Intensity - min) / dynamic range
    
    # rescale PH-Intensity
    PH_Min <- min(temp$PH)
    PH_dist <- max(temp$PH) - PH_Min                         # maxima$PH instead of temp&PH would normalize to membrane-intensity
    temp$PH_norm <- (temp$PH - PH_Min) / PH_dist
    
    # rescale Dapi-Intensity
    Dapi_Min <- min(temp$Dapi)
    Dapi_dist <- max(temp$Dapi) - Dapi_Min
    temp$Dapi_norm <- (temp$Dapi - Dapi_Min) / Dapi_dist
    
    
    # merge with existing Dataframe
    rawData <- rbind(rawData, temp)
    
    
    
    ### 3. Define regions on line profile
    # Membrane peak bins: Actin-peaks + two pixels towards cell center (Peak1+2, Peak2-2)
    Membrane <- rbind(temp[Maxima$peak_ID[1], ], 
                      temp[Maxima$peak_ID[1] + 1, ],
                      temp[Maxima$peak_ID[1] + 2, ],
                      temp[Maxima$peak_ID[2] - 2, ], 
                      temp[Maxima$peak_ID[2] - 1, ], 
                      temp[Maxima$peak_ID[2], ])
    
    # Center of cell bins: Middle between Peaks +/- 25% of cell diameter
    Cell_dia <- abs(Maxima$peak_ID[1] - Maxima$peak_ID[2])
    Cent_ID <- Maxima$peak_ID[1] + floor(Cell_dia/2)
    Cent_IDs <- c(floor(Cent_ID - Cell_dia/4):floor(Cent_ID + Cell_dia/4))
    Center <- temp[Cent_IDs, ]
    
    # Minima inside cell: sort between-peaks region by PH-intensity, keep only 25% dimmest
    # to avoid Nucleus/debris signals in cytosol
    Between_maxima <- temp[(Maxima$peak_ID[1] + 3):(Maxima$peak_ID[2] - 3), ]
    Between_maxima <- Between_maxima %>% arrange(PH)
    Between_maxima <- Between_maxima[1:floor(abs((Cell_dia - 4)/4)), ]   # -4 to account for peaks of 3px each
    
    
    ### 4. Summary Statistics Mean intensities in bins and ratio
    # Actin
    Stats[i, 3] <- median(Membrane$Actin)         # median Membrane
    Stats[i, 4] <- median(Center$Actin)           # median Center
    Stats[i, 5] <- Stats[i, 3] / Stats[i, 4]      # Membrane / Center ratio
    
    # PH domains
    Stats[i, 6] <- median(Membrane$PH)            # median Membrane
    Stats[i, 7] <- median(Center$PH)              # median Center
    Stats[i, 8] <- Stats[i, 6] / Stats[i, 7]      # Membrane / Center ratio
    Stats[i, 9] <- median(Between_maxima$PH)      # median MinCell
    Stats[i, 10] <- Stats[i, 6] / Stats[i, 9]     # Membrane / MinCell ratio
    
    
    ### 5. Quality control plots: Individual PH-plot with Dapi & F-Actin
    cols <- c("PH-Domain" = "black", "Dapi" = "grey", "F-Actin" = "#62c76b")
    
    plot <- ggplot(temp, aes(x = X, y = PH_norm)) +
      # line profiles
      geom_line(aes(y = PH_norm, color = "PH-Domain"), size = 1.5) +
      geom_line(aes(y = (Dapi_norm), color = "Dapi"), alpha = 0.5) +
      geom_line(aes(y = (Actin_norm), color = "F-Actin"), alpha = 0.5) +
      
      # Datapoints for Mean lines
      geom_point(data = Membrane, color = "#3591d1", size = 3) +
      geom_point(data = Between_maxima, color = "#f04546", size = 3) +
      geom_point(data = Center, shape = 4, color = "grey", size = 3, stroke = 1.5) +
      
      # Mean lines
      geom_line(y = mean(Membrane$PH_norm), size = 1, color = "#3591d1") +
      geom_line(y = mean(Center$PH_norm), size = 1, color = "grey") +
      geom_line(y = mean(Between_maxima$PH_norm), size = 1, color = "#f04546") +
      
      coord_cartesian(ylim = c(0, max(temp$PH_norm))) +
      
      
      # Annotate mean lines
      annotate("text", x = (max(temp$X)), y = (mean(Membrane$PH_norm) - 0.025), hjust = 1, alpha = 0.7, label = "Membrane") +
      annotate("text", x = (max(temp$X)), y = (mean(Center$PH_norm) - 0.025), hjust = 1, alpha = 0.7, label = "Center") +
      annotate("text", x = (max(temp$X)), y = (mean(Between_maxima$PH_norm) - 0.025), hjust = 1, alpha = 0.7, label = "Min in Cell") +
      
      # Legend (line profiles)
      scale_colour_manual(name = "", values = cols, guide = guide_legend(override.aes = aes(fill = NA))) +
      
      # Axes and Title
      xlab("\nPosition on line profile in µm") +
      ylab("Normalized Intensity\n") +
      ggtitle(folderlist[i]) +
      
      # Theme
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "top", axis.text = element_text(face = "bold"),
        axis.title = element_text(size = rel(1.15)), plot.title = element_text(hjust = 0.5)
      )
    
    
    plotname <- paste(folderlist[i], ".jpg", sep = "")
    ggsave(plotname,
           plot = plot, device = "jpeg", path = imagedir,
           scale = 1, width = 15, height = 15, units = "cm"
    )
  }
  
  ## separate Image_ROI-ID into two columns
  Stats <- separate(data = Stats, col = name, into = c("Image", "ROI"), sep = ".tif_")
  
  
  # save Results
  write.table(Stats, file.path(resultsdir, "Membrane_recruitment_results.txt"), sep = "\t", row.names = F)
  write.table(rawData, file.path(resultsdir, "All_peaks.txt"), sep = "\t", row.names = F)
}
