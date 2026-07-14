### Code for building up the metapopulation model.
# Written by Chrissy Hernandez
# July 13 2026
# Last updated:

# We're gonna build up from the absolute simplest model. Two patches, no dispersal:
twopatch_Umat<- matrix(data=0, nrow=8, ncol=8)
twopatch_Umat[1:4,1:4]<- Umat
twopatch_Umat[5:8,5:8]<- Umat

# fertility parameters are set in script 03_parameter_rate_setup.R (copied here for reference)
# params<- data.frame(Sc_max= rogers_cub_survival,   # max cub survival (equal for sexes), load from script rogers 1997
#                     b= beta,       # calculated from mcdonald 2016
#                     rep_K= rogers_k,          # max litter size (K), 
#                     h= 10)   # harem size per male - assume single male sufficient for groups

# I want the two patches to have different starting populations, so I'm just
# pasting together two of the initial pop vectors from the non-spatial model
n0_2patch<- c(initials[1,], initials[2,]) 

output<- proj_meta_nodisp_norem(twopatch_Umat, n0_2patch, params, stages, 
                         npatches=2, DDapply = "fertility", time=20, 
                         return.vec = TRUE)

# Okay, now two patches with dispersal: The problem with adding dispersal is
# that we want them to disperse before they mate, and then mating is
# density-dependent. So this is why we keep running into issues with the
# megamatrix approach, and switching to loops. I think we also wanted to work
# with discrete individuals, which is also a bit clunky in the matrix setting.
# So I've written a sort of hybrid approach, where we do a patch-based model
# where each patch follows an MPM.

output<- proj_meta_fixeddisp_norem(twopatch_Umat, n0_2patch, params, stages, 
                                   npatches=2, DDapply = 'fertility', time = 20,
                                   dispersal_prob=0.1, 
                                   dispersal_stages = c("Adult_f", "Adult_m"),
                                   return.vec = TRUE)



