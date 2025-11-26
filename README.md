# MANTEL-BRAY-CURTIS-FST
````
setwd("D:/Marine_Iguanas_Project/MARINE_IGUANAS/SECOND PAPER 2025/RADSEQ/mantels")
library(vegan)
list.files()

### ==== 1. Load the matrices ====

# Bray–Curtis dissimilarity matrix (OTU-based, between populations)
bray <- as.matrix(read.table("BrayCurtis_Genus_by_Population2022_26112025.txt",
                             header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
head(bray)

# Allele-sharing / genetic distance matrix
fst <- as.matrix(read.table("Allele_sharing_distance_matrix.txt",
                            header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
head(fst)

### ==== 2. Match populations between matrices ====
cat("\nBEFORE cleaning — in FST not in Bray:\n")
print(setdiff(rownames(fst), colnames(bray)))
cat("\nBEFORE cleaning — in Bray not in FST:\n")
print(setdiff(colnames(bray), rownames(fst)))

# Shared populations
common <- intersect(rownames(fst), rownames(bray))

# Reorder both matrices to the same order
fst  <- fst[common, common]
bray <- bray[common, common]

# Check
rownames(fst)
rownames(bray)

### ==== 3. Symmetrize & set diagonals ====
# Bray–Curtis is already a dissimilarity, so diag = 0
bray[upper.tri(bray)] <- t(bray)[upper.tri(bray)]
diag(bray) <- 0

fst[upper.tri(fst)] <- t(fst)[upper.tri(fst)]
diag(fst) <- 0

### ==== 4. Convert to distance objects ====
# Bray–Curtis already is a distance:
dist_bray <- as.dist(bray)

# Genetic distance: choose raw or linearized form
# Option A: linearized (recommended for distance–distance correlations)
dist_fst <- as.dist(fst / (1 - fst))

# Option B: raw allele-sharing distance (alternative)
# dist_fst <- as.dist(fst)

### ==== 5. Mantel test ====
mantel_result <- mantel(dist_bray,
                        dist_fst,
                        method = "spearman",
                        permutations = 9999)

mantel_result

````
