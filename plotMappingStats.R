projectDir = "/home/ATAC/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "ATACseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)


# plot_mapping_stats <- function(manifest, outfile = "Mapping_Stats_Spret.pdf") {
#   
#   # Set up data frame for plotting
#   stats = t(data.frame(unique = manifest$Spret.unique.reads,
#                        unmapped = manifest$Spret.unmapped.reads,
#                        dupe = manifest$Spret.mapped.reads - manifest$Spret.unique.reads))
#   row.names(stats) = c("Unique", "Unmapped", "Duplicate")
#   
#   # Plot options
#   mai = c(1.4, 0.8, 0.3, 0.3) # Margins
#   mgp = c(2, 0.5, 0.5) # The margin line (in mex units) for the axis title, axis labels and axis line. mgp[1] affects title, mgp[2:3] affects axis.
#   cols = c(rgb(77, 175, 74, maxColorValue = 255),  # Green
#            rgb(99, 99, 99, maxColorValue = 255),   # Gray
#            rgb(31, 120, 180, maxColorValue = 255)) # Blue
#   options(scipen = 9)
#   
#   # Make the barplot
#   cairo_pdf(file = outfile, height = 5, width = 5)
#   par(mai = mai, mgp = mgp)
#   barplot(stats/1e6, ylab = "Reads (millions)", col = cols, legend = rownames(stats),
#           cex.lab = 0.8, cex.axis = 0.8, cex.names = 0.8, las = 3,
#           names.arg = gsub("[[:punct:]]", " ", manifest$sample))
#   dev.off()
# }
# 
# plot_mapping_stats(files)
# 
# plot_mapping_stats <- function(manifest, outfile = "Mapping_Stats_B6.pdf") {
#   
#   # Set up data frame for plotting
#   stats = t(data.frame(unique = manifest$B6.unique.reads,
#                        unmapped = manifest$B6.unmapped.reads,
#                        dupe = manifest$B6.mapped.reads - manifest$B6.unique.reads))
#   row.names(stats) = c("Unique", "Unmapped", "Duplicate")
#   
#   # Plot options
#   mai = c(1.4, 0.8, 0.3, 0.3) # Margins
#   mgp = c(2, 0.5, 0.5) # The margin line (in mex units) for the axis title, axis labels and axis line. mgp[1] affects title, mgp[2:3] affects axis.
#   cols = c(rgb(77, 175, 74, maxColorValue = 255),  # Green
#            rgb(99, 99, 99, maxColorValue = 255),   # Gray
#            rgb(31, 120, 180, maxColorValue = 255)) # Blue
#   options(scipen = 9)
#   
#   # Make the barplot
#   cairo_pdf(file = outfile, height = 5, width = 5)
#   par(mai = mai, mgp = mgp)
#   barplot(stats/1e6, ylab = "Reads (millions)", col = cols, legend = rownames(stats),
#           cex.lab = 0.8, cex.axis = 0.8, cex.names = 0.8, las = 3,
#           names.arg = gsub("[[:punct:]]", " ", manifest$sample))
#   dev.off()
# }
# 
# plot_mapping_stats(files)

plot_mapping_stats <- function(manifest, outfile = "Mapping_Stats_Merged.pdf") {
  
  # Set up data frame for plotting
  stats = t(data.frame(unique = manifest$Merged.unique.reads,
                       unmapped = manifest$Merged.unmapped.reads,
                       dupe = manifest$Merged.mapped.reads - manifest$Merged.unique.reads))
  row.names(stats) = c("Unique", "Unmapped", "Duplicate")
  
  # Plot options
  mai = c(1.4, 0.8, 0.3, 0.3) # Margins
  mgp = c(2, 0.5, 0.5) # The margin line (in mex units) for the axis title, axis labels and axis line. mgp[1] affects title, mgp[2:3] affects axis.
  cols = c(rgb(77, 175, 74, maxColorValue = 255),  # Green
           rgb(99, 99, 99, maxColorValue = 255),   # Gray
           rgb(31, 120, 180, maxColorValue = 255)) # Blue
  options(scipen = 9)
  
  # Make the barplot
  cairo_pdf(file = outfile, height = 5, width = 5)
  par(mai = mai, mgp = mgp)
  barplot(stats/1e6, ylab = "Reads (millions)", col = cols, legend = rownames(stats),
          cex.lab = 0.8, cex.axis = 0.8, cex.names = 0.8, las = 3,
          names.arg = gsub("[[:punct:]]", " ", manifest$sample))
  dev.off()
}

plot_mapping_stats(files)

