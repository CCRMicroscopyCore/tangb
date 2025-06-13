#YAP/TAZ analysis for Binwu Tang
#02/21/20
#Dataset from 02_18. Nikon CSUW1 confocal 60x objective.
#Author: Andy D. Tran, Optical Microscopy Core, LCBG, CCR, NCI

#Libraries and themes ---------------------------------------------------------

library(tidyverse)
library(ggsci)

theme <-  theme(panel.grid.major=element_blank(), 
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                axis.line=element_line(color="black", size=1), 
                axis.ticks=element_line(color="black", size=1),
                text=element_text(size=18), 
                plot.title=element_text(size=24), 
                axis.text=element_text(color="black"),
                plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

#Upload data-------------------------------------------------------------------

path_in <- "/Users/glucetylation/Documents/CCR Confocal Core/Lalage Wakefield lab/Binwu Tang/output/yap_taz/02_18"

path_out <- "/Users/glucetylation/Documents/CCR Confocal Core/Lalage Wakefield lab/Binwu Tang/results/yap_taz/02_18"

cell_path <- file.path(path_in, "cell.csv")
cell <- read_csv(cell_path) %>% 
  rename(image = cell_image)

nuc_path <- file.path(path_in, "nucleus.csv")
nuc <- read_csv(nuc_path) %>% 
  rename(image = nuc_image) %>% 
  rename(cell_id = nuc_id)

cyto_path <- file.path(path_in, "cytoplasm.csv")
cyto <- read_csv(cyto_path) %>% 
  rename(image = cyto_image) %>% 
  rename(cell_id = cyto_id)

tmp1 <- cell %>% 
  full_join(nuc, by = c("image", "cell_id")) %>% 
  full_join(cyto, by = c("image", "cell_id")) %>% 
  na.omit() %>% 
  mutate(cyto_marker_total = cyto_area * cyto_marker_int) %>% 
  group_by(image, cell_id) %>% 
  summarise(cell_area = mean(cell_area),
            cell_sore6_int = mean(cell_sore6_int), 
            cell_marker_int = mean(cell_marker_int),
            nuc_area = mean(nuc_area),
            nuc_dapi_int = mean(nuc_dapi_int),
            nuc_marker_int = mean(nuc_marker_int),
            cyto_area = sum(cyto_area), 
            cyto_marker_total = sum(cyto_marker_total)) %>% 
  mutate(cyto_marker_int = cyto_marker_total / cyto_area)

dat <- tmp1 %>% 
  dplyr::rowwise() %>% 
  mutate(density = strsplit(image, split = " ")[[1]][2]) %>% 
  mutate(condition = case_when(
         grepl("CK", image) ~ "CK",
         grepl("neg", image) ~ "neg",
         TRUE ~ "control"))

dat2 <- dat %>% 
  dplyr::rowwise() %>% 
  mutate(marker = case_when(
    grepl("Yap", image) ~ "YAP",
    grepl("Taz", image) ~ "TAZ",
    TRUE ~ "No Ab"))

dat2 <- dat2 %>% 
  group_by() %>% 
  mutate(ratio = nuc_marker_int / cyto_marker_int)

dat2$condition <- factor(dat2$condition, levels = c("neg", "control", "CK"))

dat2$marker <- factor(dat2$marker, levels = c("No Ab", "YAP", "TAZ"))

dat2$density <- factor(dat2$density, levels = c(50, 200, 500))

dat3 <- dat2 %>% 
  filter(marker != "No Ab")

#Plot SORE6 intensities--------------------------------------------------------

tmp_plot <- dat3 %>% 
  filter(cell_sore6_int < 400)

ggplot(tmp_plot, aes(x = density, y = cell_sore6_int, col = marker)) +
  geom_jitter(width = 0.35, size = 1, alpha = 0.25, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_hline(yintercept = 160, linetype = 2, col = "red") +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(rows = vars(marker)) +
  labs(title = "SORE6 intensity", y = "Intensity [AU]")

plotname_tmp <- "sore6_intensity.png"
ggsave(plotname_tmp, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Filter out controls, segregate CSCs-------------------------------------------

dat4 <- dat3 %>% 
  filter(condition != "control") %>% 
  mutate(csc = ifelse(cell_sore6_int > 160, "CSC+", "CSC-"))

sore6 <- dat4 %>% 
  mutate(count = ifelse(csc == "CSC+", 1, 0)) %>% 
  group_by(density, marker) %>% 
  summarise(count = sum(count), total = n()) %>% 
  mutate(percent = (count/total)*100)

#Plot nuclear marker intensity-------------------------------------------------

ggplot(dat4, aes(x = density, y = nuc_marker_int, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(rows = vars(marker)) +
  labs(title = "Nuclear intensity - CK treated", y = "Intensity [AU]")

plotname_01 <- "ck_nuclear.png"
ggsave(plotname_01, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Plot cytoplasmic intensity----------------------------------------------------

ggplot(dat4, aes(x = density, y = cyto_marker_int, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(rows = vars(marker)) +
  labs(title = "Cytoplasmic intensity - CK treated", y = "Intensity [AU]")

plotname_02 <- "ck_cytoplasmic.png"
ggsave(plotname_02, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Plot nuclear/cytoplasmic ratio------------------------------------------------

ggplot(dat4, aes(x = density, y = ratio, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(rows = vars(marker)) +
  labs(title = "Nuclear/cytoplasmic ratio - CK treated", 
       y = "Nuclear/cytoplasmic ratio")

plotname_03 <- "ck_ratio.png"
ggsave(plotname_03, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Plot nuclear marker intensity by CSC state------------------------------------

ggplot(dat4, aes(x = density, y = nuc_marker_int, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(cols = vars(csc), rows = vars(marker)) +
  labs(title = "Nuclear intensity - CK treated", y = "Intensity [AU]")

plotname_04 <- "ck_nuclear_csc.png"
ggsave(plotname_04, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Plot cytoplasmic intensity by CSC state---------------------------------------

ggplot(dat4, aes(x = density, y = cyto_marker_int, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(cols = vars(csc), rows = vars(marker)) +
  labs(title = "Cytoplasmic intensity - CK treated", y = "Intensity [AU]")

plotname_05 <- "ck_cytoplasmic_csc.png"
ggsave(plotname_05, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Plot nuclear/cytoplasmic ratio by CSC state-----------------------------------

ggplot(dat4, aes(x = density, y = ratio, col = density)) +
  geom_jitter(width = 0.35, alpha = 0.5, size = 1, show.legend = F) +
  geom_violin(width = 0.5, alpha = 0.1, col = "black", show.legend = F) +
  geom_boxplot(width = 0.1, alpha = 0.1, col = "black", show.legend = F) +
  scale_color_npg() +
  scale_y_log10() +
  annotation_logticks(base = 10, sides = "l") +
  theme +
  facet_grid(cols = vars(csc), rows = vars(marker)) +
  labs(title = "Nuclear/cytoplasmic ratio - CK treated", 
       y = "Nuclear/cytoplasmic ratio")

plotname_06 <- "ck_ratio_csc.png"
ggsave(plotname_06, plot = last_plot(), path = path_out, 
       width = 11, height = 8, dpi = 600)

#Export spreadsheet------------------------------------------------------------

csvname_01 <- "yap_taz.csv"
path_out_01 <- file.path(path_out, csvname_01)

write_csv(dat4, path_out_01)

csvname_02 <- "csc_percent.csv"
path_out_02 <- file.path(path_out, csvname_02)

write_csv(sore6, path_out_02)
