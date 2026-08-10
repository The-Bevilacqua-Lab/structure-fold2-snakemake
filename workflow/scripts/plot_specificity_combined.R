suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(purrr)
  library(patchwork)
})

args          <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
output_path   <- args[2]

manifest <- read_tsv(manifest_path, col_types = cols())

all_data <- map_dfr(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]
  df  <- read_csv(row$file, col_names = c("base", "count", "specificity"),
                  col_types = "cdd", skip = 1)
  df %>%
    select(base, specificity) %>%
    mutate(group      = row$group,
           condition  = row$condition,
           replicate  = row$replicate)
})

base_colors <- c("A" = "#4C72B0", "U" = "#DD8452", "G" = "#55A868", "C" = "#C44E52")

all_data <- all_data %>%
  mutate(
    base            = if_else(base == "T", "U", base),
    base            = factor(base, levels = c("U", "G", "C", "A")),
    condition_label = if_else(condition == "plus", "+ DMS", "- DMS"),
    # replicate is used as-is (whatever ID the samplesheet assigned it) --
    # no assumption about a "repN" naming convention
    sample_label    = paste0(replicate, " (", condition_label, ")")
  )

replicates   <- sort(unique(all_data$replicate))
minus_levels <- paste0(replicates, " (- DMS)")
plus_levels  <- paste0(replicates, " (+ DMS)")
# Reverse so first level appears at top of y-axis in ggplot
all_data <- all_data %>%
  mutate(sample_label = factor(sample_label,
                               levels = rev(c(minus_levels, plus_levels))))

# One panel per comparison group (config['specificity_comparison_runs'])
panel_order <- manifest %>% distinct(group)

make_panel <- function(df, title, is_first, is_last) {
  p <- ggplot(df, aes(x = specificity, y = sample_label, fill = base)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = base_colors, breaks = c("U", "G", "C", "A"), name = NULL) +
    scale_x_continuous(limits = c(0, 1), name = "Proportion of stops",
                       breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(title = title, y = if (is_first) "Sample" else NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      panel.grid       = element_blank(),
      legend.position  = if (is_last) "right" else "none",
      legend.key.size  = unit(14, "pt"),
      legend.text      = element_text(size = 11),
    )
  if (!is_first) {
    p <- p + theme(axis.text.y  = element_blank(),
                   axis.ticks.y = element_blank())
  }
  p
}

panels <- lapply(seq_along(panel_order$group), function(i) {
  df <- filter(all_data, group == panel_order$group[i])
  make_panel(df, panel_order$group[i],
             is_first = (i == 1),
             is_last  = (i == nrow(panel_order)))
})

combined <- wrap_plots(panels, nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 13))

n_panels <- nrow(panel_order)
n_rows   <- length(replicates) * 2
fig_w    <- n_panels * 3.2 + 2.5
fig_h    <- max(3.5, n_rows * 0.45 + 1.5)

ggsave(output_path, combined, width = fig_w, height = fig_h, dpi = 300)
message("Saved specificity comparison plot to ", output_path)
