# adipocyte differentiation transcription factors
# https://www.researchgate.net/figure/Transcriptional-regulation-of-adipogenesis-Transcription-factors-Zfp423-Zfp467-EBF1_fig1_342211336
# https://www.researchgate.net/figure/A-cascade-of-transcription-factors-that-regulate-adipogenesis-PPARg-is-one-of-the-key_fig2_343634595
# https://www.nature.com/articles/nrm3198/figures/1
# https://www.sciencedirect.com/science/article/pii/S1873506116300228
# https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2015.00124/full

adipocyte_differentiation_tfs <- c(
    # positive regulation
    "Cebpa", "Cebpb", "Cebpd", "Ebf1", "Klf4", "Klf5", "Klf6", "Klf9", 
    "Pparg", "Sox6", "Stat3", "Stat5a", "Stat5b", "Zbtb16",
    
    # insignificant or not expressed
    # Foxo1, Klf15, Lmo3, Mlxipl, Mmp3, Rorb, Znf117, Nr3c1, Jak2, Krox20, Asxl2, Noct
    
    # negative regulation
    "Asxl1", "Cebpg", "Gata3", "Smad2", "Smad3", "Foxc2"
    
    # insignificant or not expressed
    # Gata2
)

adipocyte_differentiation_tfs.df <- data.frame(Sign = c(rep("Positive", 14), rep("Negative", 6)))
adipocyte_differentiation_tfs.df$Sign <- factor(adipocyte_differentiation_tfs.df$Sign, levels = c("Positive", "Negative"))
rownames(adipocyte_differentiation_tfs.df) <- adipocyte_differentiation_tfs

adipocyte_differentiation_tfs.colors <- c("Positive" = "#FF0000", "Negative" = "#436EEE")



# adipocyte markers
# https://www.frontiersin.org/journals/endocrinology/articles/10.3389/fendo.2021.599134/full

adipocyte_markers <- c(
    # AT common
    "Slc27a1", "Nrg4", "P2rx5", "Slc36a2", "Tmem120b", "Mpzl2",
    
    # WAT
    "Slc7a10", "Fabp4", "Fbxo31", "Lep", "Lpl", "Nr1h3", "Nrip1",
    "Rb1", "Rbl1", "Retn", "Serpina3k", "Tcf21", "Wfdc21",
    
    # WAT & Beige
    "Ebf3", "Hoxc8", "Hoxc9", "Pdgfra", "Asc1", 
    
    # Beige
    "Aqp7", "Tnfrsf9", "Car4", "Cd40", "Cited1", "Ear2",
    "Nr2f6", "Shox2", "Sp100", "Tbx1", "Tmem26",
    
    # BAT & Beige
    "Adrb3", "Cidea", "Clstn3", "Cox7a1", "Epsti1", "Fgf21", "Hspb7", "Kcnk3",
    "Lhx8", "Mtus1", "Ppargc1a", "Plin5", "Prdm16", "Sirt1", "Ucp1",
    
    # BAT
    "Bmp7", "Ebf2", "Ednrb", "Mir133b", "Mir206", "Myf5", "Pdk4", "Prex1", "Zic1"
)

adipocyte_markers.df <- data.frame(Gene = adipocyte_markers,
                                   Class = c(rep("AT common", 6),
                                             rep("WAT", 13),
                                             rep("WAT & Beige", 5),
                                             rep("Beige", 11),
                                             rep("BAT & Beige", 15),
                                             rep("BAT", 9)))

adipocyte_markers.df$Class <- factor(adipocyte_markers.df$Class, 
                                     levels = c("AT common", "WAT", "WAT & Beige", "Beige", "BAT & Beige", "BAT"))

adipocyte_markers.colors <- c("AT common" = "#cded56", 
                              "WAT" = "#ffec95", 
                              "WAT & Beige" = "#ffa232",
                              "Beige" = "#f0c7a9",
                              "BAT & Beige" = "#ff6132",
                              "BAT" = "#90401e")



