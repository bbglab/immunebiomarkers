color_map <- list(
    'Clinical' = '#FFFF99',
    'HLA' = '#FDB462',
    'RNA: Remaining' = '#D9D9D9',
    'RNA: TGFB' = '#BEBADA',
    'RNA: Proliferation' = '#8DD3C7',
    'RNA: T-cell' = '#FB8072',
    'Somatic: Mutation' = '#80B1D3',
    'Somatic: CNV' = '#B3DE69', 
    'Somatic: SV' = '#FCCDE5'
)

fill_map <- color_map

make_outline <- function(color_map){
    color_map[['Clinical']] <- 'black'
    color_map[['Main']] <- 'black'
    color_map
}
outline_map <- make_outline(color_map); 

alpha_map <- list(
    "Clinical" = 1,
    "HLA" = 1,
    "RNA: Remaining" = 1,
    "RNA: TGFB" = 1,
    "RNA: Proliferation" = 1,
    "RNA: T-cell" = 1,
    "Somatic: CNV" = 1,
    "Somatic: Mutation" = 1,
    "Somatic: SV" = 1,
    "Main" = 1
)
size_map <- list(
    "Clinical" = 2,
    "HLA" = 2,
    "RNA: Remaining" = 2,
    "RNA: TGFB" = 2,
    "RNA: Proliferation" = 2,
    "RNA: T-cell" = 2,
    "Somatic: CNV" = 2,
    "Somatic: Mutation" = 2,
    "Somatic: SV" = 2,
    "Main" = 8
)

x_labs <- list( 
    'bor' = 'Response Odds Ratio Estimate',
    'pfs' = '1 / PFS Hazard Estimate',
    'os' = '1 / OS Hazard Estimate',
    'lr' = 'Response Odds Ratio Estimate',
    'surv' = '1 / OS Hazard Estimate'
)
x_labeller <- function (i, x_axis) {
    if (x_axis == "cor_tmb") {
        "Correlation TMB"
    } else if (x_axis == "cor_tcell") {
        "Correlation T-cell Effector Gene Set"
    } else if (x_axis == "cor_prolif") {
        "Correlation Proliferation Gene Set"
    } else if (x_axis == "cor_pretreat") {
        "Correlation Prior Systemic Therapy"
    } else if (x_axis == "cor_tgfb") {
        "Correlation TGFB Gene Set"
    } else {
        x_labs[[i]]
    }
}
inherit_map <- list( "yes" = TRUE, "no" = FALSE)
x_translate <- function(i) ifelse( grepl("cor", i), inherit_map$yes, inherit_map$no)
