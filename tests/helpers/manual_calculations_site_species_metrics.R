# Manual Calculation Script for site_species_metrics() Validation
# 
# This script contains simple calculations for validation 
# of site_species_metrics().
# 
# Test Case 1: Perfect Partitioning ----------------------------------------
# 
# This creates a 6-site, 4-species presence-absence matrix with 2 bioregions
# where species are perfectly partitioned (no overlap between bioregions).
# 
# Expected: All species should have affinity=1.0, fidelity=1.0, indval=1.0

calculate_test_case_1 <- function() {
  
  # Create matrix: Sites 1-3 in Bioregion 1 with Species A-B only
  #                Sites 4-6 in Bioregion 2 with Species C-D only
  comat <- matrix(c(
    1, 1, 0, 0,  # Site1: Species A, B present
    1, 1, 0, 0,  # Site2: Species A, B present
    1, 1, 0, 0,  # Site3: Species A, B present
    0, 0, 1, 1,  # Site4: Species C, D present
    0, 0, 1, 1,  # Site5: Species C, D present
    0, 0, 1, 1   # Site6: Species C, D present
  ), nrow = 6, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:6)
  colnames(comat) <- paste0("Sp_", LETTERS[1:4])
  
  # Define clusters: Sites 1-3 in Bioregion 1, Sites 4-6 in Bioregion 2
  clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c("1", "1", "1", "2", "2", "2"),
    stringsAsFactors = FALSE
  )
  
  # --- Manual Calculations ---
  
  # For Affinity: A_i = R_i / Z
  #   R_i = occurrence/range size of species i in its bioregion
  #   Z = total number of sites in the bioregion
  #
  # For Fidelity: F_i = R_i / D_i
  #   R_i = occurrence/range size of species i in its bioregion  
  #   D_i = total range size of species i (across all bioregions)
  #
  # For IndVal: IndVal_i = F_i * A_i
  
  # Species A: Present in Sites 1-3 (all in Bioregion 1)
  #   R_i = 3 (occurs in 3 sites in Bioregion 1)
  #   Z = 3 (Bioregion 1 has 3 sites)
  #   D_i = 3 (total occurrences across all sites)
  #   Affinity = 3/3 = 1.0
  #   Fidelity = 3/3 = 1.0
  #   IndVal = 1.0 * 1.0 = 1.0
  
  # Species B: Present in Sites 1-3 (all in Bioregion 1)
  #   Same calculations as Species A
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  # Species C: Present in Sites 4-6 (all in Bioregion 2)
  #   R_i = 3, Z = 3, D_i = 3
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  # Species D: Present in Sites 4-6 (all in Bioregion 2)
  #   Same calculations as Species C
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  expected <- data.frame(
    Bioregion = c("1", "1", "2", "2"),
    Species = c("Sp_A", "Sp_B", "Sp_C", "Sp_D"),
    affinity = c(1.0, 1.0, 1.0, 1.0),
    fidelity = c(1.0, 1.0, 1.0, 1.0),
    indval = c(1.0, 1.0, 1.0, 1.0),
    stringsAsFactors = FALSE
  )
  
  # --- Calculate Rho for each species-bioregion combination ---
  #
  # Rho formula: ρ_ij = (n_ij - (n_i * n_j / n)) / 
  #                     sqrt(((n - n_j)/(n - 1)) * (1 - n_j/n) * (n_i * n_j / n))
  # Where:
  #   n = total number of sites
  #   n_i = number of sites where species i is present
  #   n_j = number of sites in bioregion j
  #   n_ij = number of occurrences of species i in sites of bioregion j
  
  n <- nrow(comat)  # 6 sites total
  
  # Species A in Bioregion 1:
  n_i_A <- sum(comat[, "Sp_A"] > 0)  # 3 (Species A present in 3 sites)
  n_j_1 <- sum(clusters$K_2 == 1)    # 3 (Bioregion 1 has 3 sites)
  n_ij_A1 <- sum(comat[clusters$K_2 == 1, "Sp_A"] > 0)  # 3 (all in Br1)
  
  numerator_A1 <- n_ij_A1 - (n_i_A * n_j_1 / n)
  # = 3 - (3 * 3 / 6) = 3 - 1.5 = 1.5
  
  denominator_A1 <- sqrt(((n - n_j_1) / (n - 1)) * (1 - n_j_1 / n) * 
                          (n_i_A * n_j_1 / n))
  # = sqrt(((6-3)/(6-1)) * (1 - 3/6) * (3*3/6))
  # = sqrt((3/5) * (0.5) * (1.5))
  # = sqrt(0.6 * 0.5 * 1.5)
  # = sqrt(0.45)
  # ≈ 0.6708
  
  rho_A1 <- numerator_A1 / denominator_A1  # 1.5 / 0.6708 ≈ 2.236
  
  # Species A in Bioregion 2:
  n_j_2 <- sum(clusters$K_2 == 2)    # 3
  n_ij_A2 <- sum(comat[clusters$K_2 == 2, "Sp_A"] > 0)  # 0
  
  numerator_A2 <- n_ij_A2 - (n_i_A * n_j_2 / n)  # 0 - 1.5 = -1.5
  denominator_A2 <- sqrt(((n - n_j_2) / (n - 1)) * (1 - n_j_2 / n) * 
                          (n_i_A * n_j_2 / n))  # Same as above: ≈ 0.6708
  rho_A2 <- numerator_A2 / denominator_A2  # -1.5 / 0.6708 ≈ -2.236
  
  # By symmetry:
  # Species B: Same as A (rho_B1 ≈ 2.236, rho_B2 ≈ -2.236)
  # Species C: rho_C1 ≈ -2.236, rho_C2 ≈ 2.236
  # Species D: rho_D1 ≈ -2.236, rho_D2 ≈ 2.236
  
  # Create full rho table (all species-bioregion combinations)
  expected_rho <- data.frame(
    Bioregion = rep(c("1", "2"), each = 4),
    Species = rep(c("Sp_A", "Sp_B", "Sp_C", "Sp_D"), 2),
    rho = c(2.236068, 2.236068, -2.236068, -2.236068,  # Bioregion 1
            -2.236068, -2.236068, 2.236068, 2.236068), # Bioregion 2
    stringsAsFactors = FALSE
  )
  # Note: will need tolerance (1e-6) in tests due to rounding
  
  return(list(
    comat = comat,
    clusters = clusters,
    expected = expected,
    expected_rho = expected_rho
  ))
}


# Test Case 2: Partial Overlap ------------------------------------------
#
# Creates a 4-site, 3-species matrix with 2 bioregions
# Species have different fidelity/affinity patterns with overlap

calculate_test_case_2 <- function() {
  
  # Create matrix with partial overlap between bioregions
  comat <- matrix(c(
    1, 1, 0,  # Site1: Species A, B (Bioregion 1)
    1, 0, 0,  # Site2: Species A only (Bioregion 1)
    0, 1, 1,  # Site3: Species B, C (Bioregion 2)
    0, 0, 1   # Site4: Species C only (Bioregion 2)
  ), nrow = 4, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:4)
  colnames(comat) <- c("Sp_A", "Sp_B", "Sp_C")
  
  clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c("1", "1", "2", "2"),
    stringsAsFactors = FALSE
  )
  
  # --- Manual Calculations ---
  
  # Species A: Present in Sites 1-2 (both in Bioregion 1)
  #   Assigned to Bioregion 1 (2 occurrences in Br1, 0 in Br2)
  #   R_i = 2 (occurs in 2 sites in Bioregion 1)
  #   Z = 2 (Bioregion 1 has 2 sites)
  #   D_i = 2 (total occurrences: 2 in Br1 + 0 in Br2)
  #   Affinity = 2/2 = 1.0
  #   Fidelity = 2/2 = 1.0
  #   IndVal = 1.0 * 1.0 = 1.0
  
  # Species B: Present in Sites 1, 3 (one in each bioregion)
  #   Tie-breaking: assigned to bioregion with most occurrences (1 each)
  #   The function typically assigns to the first/smallest bioregion number
  #   Assuming assigned to Bioregion 1:
  #   R_i = 1 (occurs in 1 site in Bioregion 1)
  #   Z = 2 (Bioregion 1 has 2 sites)
  #   D_i = 2 (total occurrences: 1 in Br1 + 1 in Br2)
  #   Affinity = 1/2 = 0.5
  #   Fidelity = 1/2 = 0.5
  #   IndVal = 0.5 * 0.5 = 0.25
  
  # Species C: Present in Sites 3-4 (both in Bioregion 2)
  #   Assigned to Bioregion 2 (0 occurrences in Br1, 2 in Br2)
  #   R_i = 2, Z = 2, D_i = 2
  #   Affinity = 2/2 = 1.0
  #   Fidelity = 2/2 = 1.0
  #   IndVal = 1.0 * 1.0 = 1.0
  
  expected <- data.frame(
    Bioregion = c("1", "1", "2"),
    Species = c("Sp_A", "Sp_B", "Sp_C"),
    affinity = c(1.0, 0.5, 1.0),
    fidelity = c(1.0, 0.5, 1.0),
    indval = c(1.0, 0.25, 1.0),
    stringsAsFactors = FALSE
  )
  
  # --- Calculate Rho ---
  n <- 4
  
  # Species A in Bioregion 1:
  n_i_A <- 2
  n_j_1 <- 2
  n_ij_A1 <- 2
  numerator <- n_ij_A1 - (n_i_A * n_j_1 / n)  # 2 - (2*2/4) = 2 - 1 = 1
  denominator <- sqrt(((n - n_j_1)/(n - 1)) * (1 - n_j_1/n) * (n_i_A * n_j_1/n))
  # = sqrt(((4-2)/(4-1)) * (1 - 2/4) * (2*2/4))
  # = sqrt((2/3) * (0.5) * (1))
  # = sqrt(0.3333)
  # ≈ 0.5774
  rho_A1 <- numerator / denominator  # 1 / 0.5774 ≈ 1.732
  
  # Species A in Bioregion 2:
  n_j_2 <- 2
  n_ij_A2 <- 0
  numerator <- n_ij_A2 - (n_i_A * n_j_2 / n)  # 0 - 1 = -1
  rho_A2 <- numerator / denominator  # -1 / 0.5774 ≈ -1.732
  
  # Species B (present in both bioregions):
  n_i_B <- 2
  n_ij_B1 <- 1
  numerator_B1 <- n_ij_B1 - (n_i_B * n_j_1 / n)  # 1 - 1 = 0
  denominator_B <- sqrt(((n - n_j_1)/(n - 1)) * (1 - n_j_1/n) * (n_i_B * n_j_1/n))
  # = sqrt((2/3) * 0.5 * 1) = sqrt(0.3333) ≈ 0.5774
  rho_B1 <- numerator_B1 / denominator_B  # 0 / 0.5774 = 0
  
  n_ij_B2 <- 1
  numerator_B2 <- n_ij_B2 - (n_i_B * n_j_2 / n)  # 1 - 1 = 0
  rho_B2 <- numerator_B2 / denominator_B  # 0
  
  # Species C in Bioregion 2:
  n_i_C <- 2
  n_ij_C2 <- 2
  numerator_C2 <- n_ij_C2 - (n_i_C * n_j_2 / n)  # 2 - 1 = 1
  denominator_C <- sqrt(((n - n_j_2)/(n - 1)) * (1 - n_j_2/n) * (n_i_C * n_j_2/n))
  # = sqrt(0.3333) ≈ 0.5774
  rho_C2 <- numerator_C2 / denominator_C  # 1.732
  
  # Species C in Bioregion 1:
  n_ij_C1 <- 0
  rho_C1 <- -1.732  # By symmetry
  
  expected_rho <- data.frame(
    Bioregion = rep(c("1", "2"), each = 3),
    Species = rep(c("Sp_A", "Sp_B", "Sp_C"), 2),
    rho = c(1.732051, 0, -1.732051,   # Bioregion 1
            -1.732051, 0, 1.732051),  # Bioregion 2
    stringsAsFactors = FALSE
  )
  
  return(list(
    comat = comat,
    clusters = clusters,
    expected = expected,
    expected_rho = expected_rho
  ))
}


# Test Case 3: Three Bioregions with Generalists --------------------------
#
# Creates a 9-site, 5-species matrix with 3 bioregions
# Includes specialist species (high fidelity) and generalist species 
# (low fidelity, present in multiple bioregions)

calculate_test_case_3 <- function() {
  
  # Create matrix: 3 bioregions with 3 sites each
  # Species A: Specialist for Bioregion 1
  # Species B: Generalist (present in all 3 bioregions)
  # Species C: Specialist for Bioregion 2
  # Species D: Specialist for Bioregion 3
  # Species E: Specialist for Bioregion 3
  comat <- matrix(c(
    1, 1, 0, 0, 0,  # Site1 (Br1): A, B
    1, 1, 0, 0, 0,  # Site2 (Br1): A, B
    1, 1, 0, 0, 0,  # Site3 (Br1): A, B
    0, 1, 1, 0, 0,  # Site4 (Br2): B, C
    0, 1, 1, 0, 0,  # Site5 (Br2): B, C
    0, 1, 1, 0, 0,  # Site6 (Br2): B, C
    0, 1, 0, 1, 1,  # Site7 (Br3): B, D, E
    0, 1, 0, 1, 1,  # Site8 (Br3): B, D, E
    0, 1, 0, 1, 1   # Site9 (Br3): B, D, E
  ), nrow = 9, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:9)
  colnames(comat) <- paste0("Sp_", LETTERS[1:5])
  
  clusters <- data.frame(
    ID = rownames(comat),
    K_3 = c("1", "1", "1", "2", "2", "2", "3", "3", "3"),
    stringsAsFactors = FALSE
  )
  
  # --- Manual Calculations ---
  
  # Species A: Specialist for Bioregion 1
  #   R_i = 3 (all 3 sites in Br1), Z = 3, D_i = 3
  #   Affinity = 3/3 = 1.0, Fidelity = 3/3 = 1.0, IndVal = 1.0
  
  # Species B: Generalist (present in all 3 bioregions)
  #   Assigned to bioregion with most occurrences (all have 3, so Br1 by tie)
  #   R_i = 3 (in Br1), Z = 3, D_i = 9 (3+3+3)
  #   Affinity = 3/3 = 1.0, Fidelity = 3/9 = 0.333, IndVal = 0.333
  
  # Species C: Specialist for Bioregion 2
  #   R_i = 3, Z = 3, D_i = 3
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  # Species D: Specialist for Bioregion 3
  #   R_i = 3, Z = 3, D_i = 3
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  # Species E: Specialist for Bioregion 3
  #   R_i = 3, Z = 3, D_i = 3
  #   Affinity = 1.0, Fidelity = 1.0, IndVal = 1.0
  
  expected <- data.frame(
    Bioregion = c("1", "1", "2", "3", "3"),
    Species = c("Sp_A", "Sp_B", "Sp_C", "Sp_D", "Sp_E"),
    affinity = c(1.0, 1.0, 1.0, 1.0, 1.0),
    fidelity = c(1.0, 1/3, 1.0, 1.0, 1.0),
    indval = c(1.0, 1/3, 1.0, 1.0, 1.0),
    stringsAsFactors = FALSE
  )
  
  # Note: Full rho calculations omitted for brevity 
  # (similar pattern as Test Cases 1-2)
  
  return(list(
    comat = comat,
    clusters = clusters,
    expected = expected
  ))
}


# Test Case 4: Bipartite Network (Sites and Species Clustered) ------------
#
# Creates a bipartite network where both sites AND species are assigned to 
# bioregions. This tests the Cz indices (participation coefficient and z-score)

calculate_test_case_4 <- function() {
  
  # Create bipartite interaction matrix (weighted)
  # Sites 1-2 in Bioregion 1, Sites 3-4 in Bioregion 2
  # Species A-B primarily in Bioregion 1, Species C-D primarily in Bioregion 2
  # Species E-F are generalists that occur MORE in Br1 (their assigned bioregion)
  comat <- matrix(c(
    # Sp_A, Sp_B, Sp_C, Sp_D, Sp_E, Sp_F
    3,    2,    0,    0,    3,    2,   # Site1 (Br1): strong with A,B,E,F
    2,    4,    0,    0,    2,    3,   # Site2 (Br1): strong with A,B,E,F (varying weights)
    0,    0,    2,    3,    1,    0,   # Site3 (Br2): strong with C,D; weak with E
    0,    0,    3,    2,    0,    1    # Site4 (Br2): strong with C,D; weak with F
  ), nrow = 4, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:4)
  colnames(comat) <- paste0("Sp_", LETTERS[1:6])
  
  # Clusters for sites
  site_clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c("1", "1", "2", "2"),
    stringsAsFactors = FALSE
  )
  
  # Clusters for species
  # A,B are specialists in Br1; C,D are specialists in Br2
  # E,F are assigned to Br1 (by majority rule, equal occurrence so arbitrary)
  species_clusters <- data.frame(
    ID = colnames(comat),
    K_2 = c("1", "1", "2", "2", "1", "1"),
    stringsAsFactors = FALSE
  )
  
  # --- Manual Cz Calculations ---
  #
  # Participation Coefficient: C_i = 1 - sum((k_is / k_i)^2)
  #   k_is = number of links of node i to nodes in bioregion s
  #   k_i = total degree of node i
  #
  # Within-bioregion Degree Z-score: z_i = (k_i - mean(k_si)) / sd(k_si)
  #   k_i = number of links of node i to nodes in its bioregion
  #   mean(k_si) = average degree within bioregion
  #   sd(k_si) = standard deviation of degrees within bioregion
  
  # --- SITES ---
  
  # Site1 (in Br1):
  #   Total degree (sum of weights): 3 + 2 + 0 + 0 + 3 + 2 = 10
  #   Links to Br1 species (Sp_A, Sp_B, Sp_E, Sp_F): 3 + 2 + 3 + 2 = 10
  #   Links to Br2 species (Sp_C, Sp_D): 0 + 0 = 0
  #   C_Site1 = 1 - ((10/10)^2 + (0/10)^2) = 1 - (1 + 0) = 0
  #   (All links to species in own bioregion)
  
  # Site2 (in Br1):
  #   Total degree: 2 + 4 + 0 + 0 + 2 + 3 = 11
  #   Links to Br1 species: 2 + 4 + 2 + 3 = 11
  #   Links to Br2 species: 0 + 0 = 0
  #   C_Site2 = 1 - ((11/11)^2 + (0/11)^2) = 0
  
  # Within Br1 sites (Site1, Site2):
  #   Site1 degree within Br1 species: 10
  #   Site2 degree within Br1 species: 11
  #   mean = (10 + 11) / 2 = 10.5
  #   sd = sqrt(((10-10.5)^2 + (11-10.5)^2) / (2-1))
  #      = sqrt((0.25 + 0.25) / 1) = sqrt(0.5) ≈ 0.7071
  #   z_Site1 = (10 - 10.5) / 0.7071 ≈ -0.7071
  #   z_Site2 = (11 - 10.5) / 0.7071 ≈ 0.7071
  
  # Site3 (in Br2):
  #   Total degree: 0 + 0 + 2 + 3 + 1 + 0 = 6
  #   Links to Br2 species (Sp_C, Sp_D): 2 + 3 = 5
  #   Links to Br1 species (Sp_A, Sp_B, Sp_E, Sp_F): 0 + 0 + 1 + 0 = 1
  #   C_Site3 = 1 - ((5/6)^2 + (1/6)^2) 
  #          = 1 - (0.6944 + 0.0278) 
  #          = 1 - 0.7222 = 0.2778
  
  # Site4 (in Br2):
  #   Total degree: 0 + 0 + 3 + 2 + 0 + 1 = 6
  #   Links to Br2 species: 3 + 2 = 5
  #   Links to Br1 species: 0 + 0 + 0 + 1 = 1
  #   C_Site4 = 1 - ((5/6)^2 + (1/6)^2) = 0.2778
  
  # Within Br2 sites (Site3, Site4):
  #   Site3 degree within Br2 species: 5
  #   Site4 degree within Br2 species: 5
  #   mean = (5 + 5) / 2 = 5
  #   sd = 0
  #   z_Site3 = (5 - 5) / 0 = undefined → set to 0
  #   z_Site4 = (5 - 5) / 0 = undefined → set to 0
  
  # --- SPECIES ---
  
  # Sp_A (specialist in Br1):
  #   Total degree: 3 + 2 + 0 + 0 = 5
  #   Links to Br1 sites (Site1, Site2): 3 + 2 = 5
  #   Links to Br2 sites (Site3, Site4): 0 + 0 = 0
  #   C_SpA = 1 - ((5/5)^2 + (0/5)^2) = 0
  
  # Sp_B (specialist in Br1):
  #   Total degree: 2 + 4 + 0 + 0 = 6
  #   Links to Br1 sites: 2 + 4 = 6
  #   Links to Br2 sites: 0 + 0 = 0
  #   C_SpB = 1 - ((6/6)^2 + (0/6)^2) = 0
  
  # Sp_C (specialist in Br2):
  #   Total degree: 0 + 0 + 2 + 3 = 5
  #   Links to Br2 sites: 2 + 3 = 5
  #   Links to Br1 sites: 0 + 0 = 0
  #   C_SpC = 0
  
  # Sp_D (specialist in Br2):
  #   Total degree: 0 + 0 + 3 + 2 = 5
  #   Links to Br2 sites: 3 + 2 = 5
  #   Links to Br1 sites: 0 + 0 = 0
  #   C_SpD = 0
  
  # Within Br2 species (Sp_C, Sp_D):
  #   Sp_C degree within Br2 sites: 5
  #   Sp_D degree within Br2 sites: 5
  #   mean = (5 + 5) / 2 = 5
  #   sd = 0
  #   z_SpC = (5 - 5) / 0 = undefined → set to 0
  #   z_SpD = (5 - 5) / 0 = undefined → set to 0
  
  # Sp_E (generalist, assigned to Br1 because more occurrences there):
  #   Total degree: 3 + 2 + 1 + 0 = 6
  #   Links to Br1 sites (Site1, Site2): 3 + 2 = 5
  #   Links to Br2 sites (Site3, Site4): 1 + 0 = 1
  #   C_SpE = 1 - ((5/6)^2 + (1/6)^2) 
  #         = 1 - (0.6944 + 0.0278) 
  #         = 1 - 0.7222 = 0.2778
  
  # Sp_F (generalist, assigned to Br1 because more occurrences there):
  #   Total degree: 2 + 3 + 0 + 1 = 6
  #   Links to Br1 sites (Site1, Site2): 2 + 3 = 5
  #   Links to Br2 sites (Site3, Site4): 0 + 1 = 1
  #   C_SpF = 1 - ((5/6)^2 + (1/6)^2) = 0.2778
  
  # Within Br1 species (Sp_A, Sp_B, Sp_E, Sp_F):
  #   Sp_A degree within Br1 sites: 5
  #   Sp_B degree within Br1 sites: 6
  #   Sp_E degree within Br1 sites: 5
  #   Sp_F degree within Br1 sites: 5
  #   mean = (5 + 6 + 5 + 5) / 4 = 5.25
  #   sd = sqrt(((5-5.25)^2 + (6-5.25)^2 + (5-5.25)^2 + (5-5.25)^2) / (4-1))
  #      = sqrt((0.0625 + 0.5625 + 0.0625 + 0.0625) / 3)
  #      = sqrt(0.75 / 3) = sqrt(0.25) = 0.5
  #   z_SpA = (5 - 5.25) / 0.5 = -0.5
  #   z_SpB = (6 - 5.25) / 0.5 = 1.5
  #   z_SpE = (5 - 5.25) / 0.5 = -0.5
  #   z_SpF = (5 - 5.25) / 0.5 = -0.5
  
  expected_cz <- data.frame(
    Node = c(paste0("Site", 1:4), paste0("Sp_", LETTERS[1:6])),
    Bioregion = c("1", "1", "2", "2", "1", "1", "2", "2", "1", "1"),
    C = c(0, 0, 5/18, 5/18,                      # Sites: 1 - ((5/6)^2 + (1/6)^2) = 1 - 26/36 = 10/36 = 5/18
          0, 0, 0, 0, 5/18, 5/18),               # Species
    z = c(-1/sqrt(2), 1/sqrt(2), 0, 0,           # Sites: ±0.7071068
          -0.5, 1.5, 0, 0, -0.5, -0.5),          # Species
    stringsAsFactors = FALSE
  )
  
  return(list(
    comat = comat,
    site_clusters = site_clusters,
    species_clusters = species_clusters,
    expected_cz = expected_cz
  ))
}


# Test Case 5: Small Network for Rho Calculation --------------------------
#
# Minimal 3-site, 2-species, 2-bioregion case for detailed rho validation

calculate_test_case_5 <- function() {
  
  # Create minimal matrix for manual rho verification
  comat <- matrix(c(
    1, 0,  # Site1 (Br1): Species A only
    1, 1,  # Site2 (Br2): Species A and B
    0, 1   # Site3 (Br2): Species B only
  ), nrow = 3, byrow = TRUE)
  rownames(comat) <- paste0("Site", 1:3)
  colnames(comat) <- c("Sp_A", "Sp_B")
  
  clusters <- data.frame(
    ID = rownames(comat),
    K_2 = c("1", "2", "2"),
    stringsAsFactors = FALSE
  )
  
  # --- Manual Rho Calculations ---
  #
  # Rho formula: ρ_ij = (n_ij - (n_i * n_j / n)) / 
  #                     sqrt(((n - n_j)/(n - 1)) * (1 - n_j/n) * (n_i * n_j / n))
  
  n <- 3  # total sites
  
  # Species A in Bioregion 1:
  n_i_A <- 2      # Sp_A present in Sites 1, 2
  n_j_1 <- 1      # Bioregion 1 has 1 site (Site1)
  n_ij_A1 <- 1    # Sp_A occurs once in Br1 (Site1)
  
  numerator <- n_ij_A1 - (n_i_A * n_j_1 / n)
  # = 1 - (2 * 1 / 3) = 1 - 0.6667 = 0.3333
  
  denominator <- sqrt(((n - n_j_1) / (n - 1)) * (1 - n_j_1 / n) * 
                       (n_i_A * n_j_1 / n))
  # = sqrt(((3 - 1) / (3 - 1)) * (1 - 1/3) * (2 * 1 / 3))
  # = sqrt((2/2) * (2/3) * (2/3))
  # = sqrt(1 * 0.6667 * 0.6667)
  # = sqrt(0.4444)
  # ≈ 0.6667
  
  rho_A1 <- numerator / denominator  # 0.3333 / 0.6667 ≈ 0.5
  
  # Species A in Bioregion 2:
  n_j_2 <- 2      # Bioregion 2 has 2 sites
  n_ij_A2 <- 1    # Sp_A occurs once in Br2 (Site2)
  
  numerator_A2 <- n_ij_A2 - (n_i_A * n_j_2 / n)
  # = 1 - (2 * 2 / 3) = 1 - 1.3333 = -0.3333
  
  denominator_A2 <- sqrt(((n - n_j_2) / (n - 1)) * (1 - n_j_2 / n) * 
                          (n_i_A * n_j_2 / n))
  # = sqrt(((3 - 2) / 2) * (1 - 2/3) * (2 * 2 / 3))
  # = sqrt((1/2) * (1/3) * (4/3))
  # = sqrt(0.5 * 0.3333 * 1.3333)
  # = sqrt(0.2222)
  # ≈ 0.4714
  
  rho_A2 <- numerator_A2 / denominator_A2  # -0.3333 / 0.4714 ≈ -0.707
  
  # Species B in Bioregion 1:
  n_i_B <- 2      # Sp_B present in Sites 2, 3
  n_ij_B1 <- 0    # Sp_B does not occur in Br1
  
  numerator_B1 <- n_ij_B1 - (n_i_B * n_j_1 / n)
  # = 0 - (2 * 1 / 3) = -0.6667
  
  denominator_B1 <- sqrt(((n - n_j_1) / (n - 1)) * (1 - n_j_1 / n) * 
                          (n_i_B * n_j_1 / n))
  # = sqrt(1 * 0.6667 * 0.6667) ≈ 0.6667
  
  rho_B1 <- numerator_B1 / denominator_B1  # -0.6667 / 0.6667 ≈ -1.0
  
  # Species B in Bioregion 2:
  n_ij_B2 <- 2    # Sp_B occurs twice in Br2 (Sites 2, 3)
  
  numerator_B2 <- n_ij_B2 - (n_i_B * n_j_2 / n)
  # = 2 - (2 * 2 / 3) = 2 - 1.3333 = 0.6667
  
  denominator_B2 <- sqrt(((n - n_j_2) / (n - 1)) * (1 - n_j_2 / n) * 
                          (n_i_B * n_j_2 / n))
  # ≈ 0.4714 (same as for Species A in Br2)
  
  rho_B2 <- numerator_B2 / denominator_B2  # 0.6667 / 0.4714 ≈ 1.414
  
  expected_rho <- data.frame(
    Bioregion = c("1", "1", "2", "2"),
    Species = c("Sp_A", "Sp_B", "Sp_A", "Sp_B"),
    rho = c(0.5, -1.0, -0.7071068, 1.414214),
    stringsAsFactors = FALSE
  )
  
  return(list(
    comat = comat,
    clusters = clusters,
    expected_rho = expected_rho
  ))
}
