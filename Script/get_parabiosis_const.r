# get_parabiosis_const.r

get_parabiosis_const <- function( )
{ 
  CELLTYPES <- c( "Tconv", "CD8", "Treg" )
  
  TISSUES <- c( "Adrenals", "BoneMarrow", "Brain", "IEL", "Kidney",
                "Liver", "LN", "LPL", "Lung", "MLN", "Muscle",
                "Pancreas", "PP", "Skin", "Spleen", "WAT" )                     
  
  TISSUES.NONLYMPH <- c( "Adrenals", "Brain", "Kidney", "Liver", 
                         "Lung", "Muscle", "Pancreas", "Skin", "WAT" )   
  
  TISSUES.GUT <- c( "PP", "IEL", "LPL" )
  
  TISSUES.NODES <- c( "LN", "MLN" )
  
  TISSUES.SHORT.CIRCUIT <- c( "Spleen", "BoneMarrow" )
  
  tissue.all.ordered <- c( 
    "Blood", 
    "BoneMarrow", "Spleen", "LN", "MLN", "PP", 
    "Adrenals", "Brain", "Kidney", "Lung", "Liver", 
    "Muscle", "Pancreas", "Skin", "WAT", 
    "IEL", "LPL" )
  
  tissue.all.ordered.full <- c( 
    "Blood", 
    "Bone Marrow", "Spleen", "Lymph Nodes", 
    "Mesenteric Lymph Nodes", "Peyer's Patches", 
    "Adrenals", "Brain", "Kidney", "Lungs", "Liver", 
    "Muscle", "Pancreas", "Skin", "White Adipose Tissue", 
    "Intraepithelial Lymphocytes", "Lamina Propria Lymphocytes" )    
  
  tibble( tissue.all.ordered = tissue.all.ordered ) %>% 
    mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
    filter( tissue.all.ordered != "Blood" ) %>% 
    mutate( tissue = factor( tissue.all.ordered, levels = tissue.all.ordered ) ) %>%
    mutate( group = "Non-lymphoid") %>% 
    mutate( group = ifelse( tissue %in% c( 
      "BoneMarrow", "Spleen", "LN", "MLN", "PP"), "Lymphoid", group ) ) %>%  
    mutate( group = ifelse( tissue %in% c( "IEL", "LPL" ), "GALT", group ) ) %>% 
    mutate( group = factor( 
      group, levels = c( "Lymphoid", "Non-lymphoid", "GALT" ) ) ) ->
    dTissueAllOrderedGroup
  
  tibble( tissue.all = tissue.all.ordered ) %>% 
    mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
    filter( tissue.all != "Blood" ) %>% 
    mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
    mutate( tissue.group = "Non-lymphoid" ) %>% 
    mutate( tissue.group = ifelse( 
      f.tissue %in% c( "BoneMarrow", "Spleen", "LN", "MLN", "PP"), 
      "Lymphoid", tissue.group )) %>%
    mutate( tissue.group = ifelse( 
      f.tissue %in% c( "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
    mutate( f.tissue.group = factor( 
      tissue.group, levels = c( "Lymphoid", "Non-lymphoid", "GALT" ) ) ) %>% 
    mutate( f.tissue.group.longer = 
              dplyr::recode( f.tissue.group, 
                             Lymphoid = "Lymphoid\ntissues", 
                             `Non-lymphoid` = "Non-lymphoid\ntissues",
                             GALT = "Gut-assoc. lymphoid\ntissues" ) ) ->
    dTissueAllOrderedGroup.f
  
  
  tibble( tissue.all = tissue.all.ordered ) %>% 
    mutate( tissue.all.ordered.full = tissue.all.ordered.full ) %>% 
    # filter( tissue.all != "Blood" ) %>% 
    mutate( f.tissue = factor( tissue.all, levels = tissue.all ) ) %>%
    mutate( tissue.group = "Non-lymphoid" ) %>% 
    mutate( tissue.group = ifelse( 
      f.tissue %in% c( "BoneMarrow", "Spleen", "LN", "MLN", "PP"), 
      "Lymphoid", tissue.group )) %>%
    mutate( tissue.group = ifelse( 
      f.tissue %in% c( "IEL", "LPL" ), "GALT", tissue.group ) ) %>% 
    mutate( tissue.group = ifelse( 
      f.tissue %in% c( "Blood" ), "Blood", tissue.group ) ) %>% 
    mutate( f.tissue.group = factor( 
      tissue.group, levels = c( "Lymphoid", "Non-lymphoid", "GALT", "Blood" ) ) ) %>% 
    mutate( f.tissue.group.longer = 
              dplyr::recode( f.tissue.group, 
                             Lymphoid = "Lymphoid\ntissues", 
                             `Non-lymphoid` = "Non-lymphoid\ntissues",
                             GALT = "Gut-assoc. lymphoid\ntissues",
                             Blood = "Blood" ), ) ->
    dTissueAllOrderedGroupBlood.f
  
  pmc <- list( CELLTYPES,
               TISSUES, TISSUES.NONLYMPH, TISSUES.GUT, 
               TISSUES.NODES, TISSUES.SHORT.CIRCUIT,
               tissue.all.ordered, tissue.all.ordered.full, 
               dTissueAllOrderedGroup, dTissueAllOrderedGroup.f, 
               dTissueAllOrderedGroupBlood.f )
  names( pmc ) <- c( "CELLTYPES",
                     "TISSUES", "TISSUES.NONLYMPH", "TISSUES.GUT", 
                     "TISSUES.NODES", "TISSUES.SHORT.CIRCUIT",
                     "tissue.all.ordered", "tissue.all.ordered.full", 
                     "dTissueAllOrderedGroup", "dTissueAllOrderedGroup.f",
                     "dTissueAllOrderedGroupBlood.f" )
  return( pmc )
}
