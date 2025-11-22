setwd("D:/Salinity Stress Research/Systematic data-to publish/Profil hidrofobisitas")

# load libraries
library ("seqinr")
library("dplyr")
library("ggplot2")
library("reshape2")
library("zoo")
library("viridis")

# load file
fasta_file <- "Analyze_Collection data-protein.fasta"  # ubah sesuai nama file kamu
seqs <- read.fasta(file = fasta_file, seqtype = "AA", as.string = TRUE)

# --- STEP 2: Bersihkan nama & sekuens ---
# Hapus karakter aneh di nama (header)
names(seqs) <- gsub("[^A-Za-z0-9_]", "_", names(seqs))

# Hapus karakter non-amino acid dan ubah ke huruf besar
seqs <- lapply(seqs, function(x) gsub("[^A-Z]", "", toupper(getSequence(x, as.string = TRUE))))

# --- STEP 4: Hydrophobicity function (Kyte–Doolittle scale) ---
get_hydro_profile <- function(sequence, name, window = 19) {
  kd_scale <- c(
    A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5, 
    Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I = 4.5, 
    L = 3.8, K = -3.9, M = 1.9, F = 2.8, P = -1.6, 
    S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2
  )
  seq_vec <- unlist(strsplit(sequence, split = ""))
  seq_vec[!seq_vec %in% names(kd_scale)] <- "X"
  kd_scale[["X"]] <- 0
  hydro_values <- sapply(seq_vec, function(res) kd_scale[[res]])
  smooth_hydro <- zoo::rollmean(hydro_values, k = window, fill = NA, align = "center")
  df <- data.frame(Position = 1:length(smooth_hydro),
                   Hydrophobicity = smooth_hydro,
                   Species = name)
  return(df)
}

# --- STEP 5: Compute hydrophobicity for all sequences ---
profiles <- lapply(names(seqs), function(nm) get_hydro_profile(seqs[[nm]], nm))
profiles_df <- do.call(rbind, profiles)

# --- Save results to CSV ---
write.csv(profiles_df, "Hydrophobicity_profiles.csv", row.names = FALSE)

# --- STEP 7: Compute mean GRAVY per species ---
gravy_summary <- aggregate(Hydrophobicity ~ Species, data = profiles_df, mean, na.rm = TRUE)
write.csv(gravy_summary, "GRAVY_summary.csv", row.names = FALSE)
print(gravy_summary)

# --- STEP 8: Label cantik (untuk grafik publikasi) ---
species_labels <- c(
  "Oryza_sativa" = "Oryza sativa",
  "Indosasa_sinica" = "Indosasa sinica",
  "Aegilops_speltoides" = "Aegilops speltoides",
  "Triticum_turgidum_subsp._durum" = "Triticum turgidum subsp. durum",
  "Triticum_aestivum" = "Triticum aestivum",
  "Triticum_monococcum" = "Triticum monococcum",
  "Aeluropus_littoralis" = "Aeluropus littoralis",
  "Suaeda_salsa" = "Suaeda salsa",
  "Suaeda_japonica" = "Suaeda japonica",
  "Salicornia_dolichostachya" = "Salicornia dolichostachya",
  "Mesembryanthemum_crystallinum" = "Mesembryanthemum crystallinum",
  "Populus_ilicifolia" = "Populus ilicifolia",
  "Populus_alba" = "Populus alba",
  "Populus_deltoides" = "Populus deltoides",
  "Populus_tremuloides" = "Populus tremuloides",
  "Populus_euphratica" = "Populus euphratica",
  "Eutrema_salsugineum" = "Eutrema salsugineum"
)

# --- STEP 9: Plot profil hidrofobisitas ---
ggplot(profiles_df, aes(x = Position, y = Hydrophobicity, color = Species)) +
  geom_line(size = 0.8, alpha = 0.9) +
  theme_minimal() +
  scale_color_manual(
    values = viridis::viridis(
      length(unique(profiles_df$Species)),
      option = "E",   # Pilihan terbaik untuk >15 warna
      direction = -1
    ),
    labels = species_labels[names(species_labels) %in% unique(profiles_df$Species)]
  ) +
  labs(
    title = "Hydrophobicity Profiles (Kyte–Doolittle, window = 19)",
    x = "Residue Position",
    y = "Hydrophobicity",
    color = "Species"
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(face = "italic", size = 9),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85")
  ) +
  guides(color = guide_legend(ncol = 3))