# script for progress update figures 
# 

## old params----
# fixed dispersal proj without removals, old parameters 
ggsave(filename = "meta_fixeddisp_norem_oldparams.png",
       plot = meta_fixeddisp_norem_oldparams,
       device = "png",
       path = here("Figs"), 
       bg = "white")

# density dependent dispersal proj without removals, old parameters 
ggsave(filename = "meta_DDdisp_norem_oldparams.png",
       plot = meta_DDdisp_norem_oldparams,
       device = "png",
       path = here("Figs"), 
       bg = "white")



# new params ----
# fixed dispersal proj without removals, new parameters (tighter dd)
ggsave(filename = "meta_fixeddisp_norem_newparams.png",
       plot = meta_fixeddisp_norem_newparams,
       device = "png",
       path = here("Figs"), 
       bg = "white")

# density dependent dispersal proj without removals, old parameters (tighter dd)
ggsave(filename = "meta_DDdisp_norem_newparams.png",
       plot = meta_DDdisp_norem_newparams,
       device = "png",
       path = here("Figs"), 
       bg = "white")


## Removal projections ----
# no dispersal, removals at year 10, old params
# Summary fig - shows that when initial vecs equal, removals across multiple years introduce stochasticity
ggsave(filename = "removal_ex_oldparams.png",
       plot = removal_summary_plot,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# no dispersal, removals, new parameters (tighter dd)
ggsave(filename = "meta_nodisp_rem_newparams.png",
       plot = meta_nodisp_rem_newplot,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# density dependent dispersal, removals, old parameters 
ggsave(filename = "meta_DDdisp_rem_oldparams.png",
       plot = meta_DDdisp_rem_oldplot,
       device = "png",
       path = here("Figs"), 
       bg = "white")


# DD dispersal, removals, new parameters (tighter dd)
ggsave(filename = "meta_DDdisp_rem_newparams.png",
       plot = meta_DDdisp_rem_newplot,
       device = "png",
       path = here("Figs"), 
       bg = "white")


