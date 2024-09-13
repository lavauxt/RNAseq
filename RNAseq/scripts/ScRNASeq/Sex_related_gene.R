# from PBMC_example$bulk_DE_cors # 59 sex-specific DEGs in bulk PBMC (up-regulated = female-biased)
# https://cran.r-project.org/web/packages/scMappR/vignettes/scMappR_Vignette.html

# Small list of sex related gene
#male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1")
#female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Utx", "Rps4x", "Jarid1c", "Prdx4")

# Extended list of sex related gene
#female.features <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Rps4x", "Jarid1c", "Prdx4", "Med14","Zfx", "Ocrl", "Ddx39b", "Gpm6b", "Fthl17", "Esr1", "Foxl2", "Cyp19a1", "Wnt4", "Amhr2", "Eda2r")
#male.features <- c("Sry", "Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Usp9y", "Zfy1", "Zfy2", "Rps4y1", "Smcy", "Eif1ay", "Ssty1", "Ssty2", "Rbmy1a1", "Sox30", "Sox9", "Amh", "Insl3", "Hsd17b3", "Cyp17a1")

# male.features <- c("Abca13", "Adam22", "Arhgap31-as1", "Bcorp1", "Bpi", "Camp", "Cdy4p", "Ceacam19","Ceacam6", "Ceacam8", "Celsr2", "Crisp3", "Ctd-2308n23.2", "Ddx3y", "Defa3", "Defa4","Eif1ay", "Fbn1", "Gstt1", "Inhba", "Kalp", "Kdm5d", "Kif9-as1", "Lcn2", "Linc00278","Ltf", "Mmp8", "Mpo", "Prky", "Rn7skp282", "Rp11-424g14.1", "Rp6-99m1.2", "Rps4y1", "Slc1a2", "Tmsb4y", "Ttty14", "Ttty15", "Txlng2p", "Uty", "Usp9y", "Xist", "Ythdf3-as1", "Zfy", "Zfy-as1", "Znf839p1")
# female.features <- c("Ac084018.1", "Ighg3", "Ighv1-46", "Igkv1-27", "Igkv3-15", "Igkv3d-11", "Iglv1-40", "Iglv1-47", "Iglv3-1", "Iglv6-57", "Ift81", "Kdm5c", "Prkx", "Tsix", "Xist")

# curated
female.features_mouse <- c("Xist", "Tsix", "Kdm6a", "Eif2s3x", "Ddx3x", "Rps4x", "Prdx4", "Med14","Zfx", "Ocrl", "Ddx39b", "Gpm6b", "Esr1", "Foxl2", "Wnt4", "Amhr2", "Eda2r")
male.features_mouse <- c("Eif2s3y", "Ddx3y", "Uty", "Kdm5d", "Sox30", "Sox9", "Amh", "Insl3")

female_features_human <- c("IGKV1-27", "IGHG3", "IGHV1-46", "IGKV3-15", "IGKV3D-11", "IGLV3-1", "IGLV1-47", "IFT81", "IGLV1-40")
male_features_human <- c("XIST", "TSIX", "DDX3Y", "PRKY", "KDM5D", "RPS4Y1", "UTY", "USP9Y", "EIF1AY", "TXLNG2P", "ZFY", "TTTY15", "BCORP1", 
"LINC00278", "RP11-424G14.1", "TMSB4Y", "LTF", "RN7SKP282", "CEACAM8", "ZNF839P1", "ABCA13", "ZFY-AS1", "TTTY14", "DEFA4", "CEACAM6", "DEFA3", "BPI", "IGLV6-57", "RP6-99M1.2", 
 "KALP", "CAMP", "CRISP3", "KDM5C", "MPO", "PRKX", "KIF9-AS1", "ADAM22", "ARHGAP31-AS1", "CELSR2", "CDY4P", "AC084018.1", 
"LCN2", "CTD-2308N23.2", "INHBA", "FBN1", "SLC1A2", "GSTT1", "CEACAM19", "MMP8", "YTHDF3-AS1")