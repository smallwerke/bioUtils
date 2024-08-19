#' CT to RE
#'
#' @description
#' A simple script that takes in a tibble of CT values as produced by the RQdeltaCT package and
#' calculates the deltaCT, 2^-(delta delta CT), and relative expression values and returns them as
#' separate tibbles in a tibble that also contains lists of the mean for the control group CT, and
#' the mean of the control group for 2^-(delta delta CT) values. In addition to the CT values the
#' name of the control group and the desired housekeeping gene name must be supplied.
#'
#' Returned Tibble:
#'
#' CT_ctrl_mean: a list of the mean CT values for the control group of each gene
#'
#' dCT: a tibble of the deltaCT values
#'
#' ddCT: a tibble of the 2^-(delta delta CT) values
#'
#' ddCT_ctrl_mean: a list of the mean ddCT values for the control group of each gene
#'
#' RE: relative expression for each sample as the fold change compared to the control group
#'
#' Currently NO error checking of the inputs is being carried out!
#'
#' @param CTvals the ct values under examination from the RQdeltaCT package
#' @param ctrl string defining the control group in the CTvals tibble (Group col)
#' @param hskp string defining the housekeeping gene in the CTvals tibble (which Gene col)
#' @param CTund optional the ct values initially loaded by RQdletaCT before any filtering was carried out, this is used to set individual RE values to NA if Undetermined in original import
#'
#' @return a tibble containing 5 rows with the calculated data
#' @export
#'
#' @importFrom magrittr "%>%"
#' @importFrom data.table ":="
#' @importFrom rlang .data
CTtoRE <- function(CTvals, ctrl, hskp, CTund = NULL) {

    # for testing in scripts..
    #detach("package:bioUtils", unload=TRUE)
    #library(bioUtils, verbose=TRUE)

    CTvals.ctrl.mean = ""
    CTvals.dCT = ""
    CTvals.ddCT = ""
    CTvals.ddCT.ctrl.mean = ""
    CTvals.RE = ""
    CTvals.RE.Und = ""

    #print("RUNNING CTtoRE")
    # determine the mean CT for the control group
    CTvals.ctrl.mean = CTvals %>% dplyr::filter(.data$Group == ctrl) %>% dplyr::summarize(dplyr::across(dplyr::where(is.numeric), mean)) %>% dplyr::select(-dplyr::any_of(c("Group", "Sample")))

    # take the submitted CT values
    CTvals.dCT <- RQdeltaCT::delta_Ct(data = CTvals,
                           normalise = TRUE,
                           ref = hskp,
                           transform = FALSE)

    # calculate the 2^(-ddCT) values
    CTvals.ddCT = CTvals.dCT
    for (c in colnames(CTvals.ctrl.mean)) {
        #print(paste0("working on col ", c))
        if (c[1] != hskp) {
            CTvals.ddCT = CTvals.ddCT %>% dplyr::mutate(dplyr::across(c(c[1]), ~ 2^(-(.x - ((CTvals.ctrl.mean[c[1]] - CTvals.ctrl.mean[hskp])[1,1]) )) ))
        }
    }

    # calculate the average of all control groups for each gene
    CTvals.ddCT.ctrl.mean = CTvals.ddCT %>% dplyr::filter(.data$Group == ctrl) %>% dplyr::summarize(dplyr::across(dplyr::where(is.numeric), mean)) %>% dplyr::select(-dplyr::any_of(c("Group", "Sample")))
    # make a copy of the ddCt values and then replace them column by column with RE
    CTvals.RE = CTvals.ddCT
    for (c in colnames(CTvals.ddCT.ctrl.mean)) {
        if (c[1] != hskp) {
            CTvals.RE = CTvals.RE %>% dplyr::mutate(dplyr::across(c(c[1]), ~ (.x / pull(CTvals.ddCT.ctrl.mean[c[1]]) ) ))
        }
    }

    if (!is.null(CTund)) {

        CTvals.RE.Und = CTvals.RE

        # pull out the Undetermined samples from the CT input and then
        # set those genes & samples to NA
        # since the CTund is likely the initial 'raw' input table make sure that the gene & samples
        # are still being used in the RE table otherwise there will be an error...
        CT.UND = CTund %>% dplyr::filter(.data$Ct == "Undetermined", .data$Gene != hskp, .data$Sample %in% CTvals.RE$Sample)

        # use a basic loop to go through the above filtered Undetermined samples
        # check that this Gene is in the RE table and tehn set it to NA
        for (c in 1:nrow(CT.UND)) {
            #gene = CT.UND[c,2]
            #print(gene)
            #sample = CT.UND[c,1][1]
            #print(sample)

            # if the gene isn't in the CTRE tibble move on!
            #if (has_name(CTvals.RE, CT.UND[c,2])) { print("FOUND GENE") }
            if (!rlang::has_name(CTvals.RE, CT.UND[c,2])) { next }

            CTvals.RE.Und = dplyr::rows_update(CTvals.RE.Und, tibble::tibble(Sample=CT.UND[c,1], !!CT.UND[c,2] := NA), by="Sample")
        }
    }

    # return the results
    return(dplyr::tribble(~Name, ~Data,
                   "CT_ctrl_mean", CTvals.ctrl.mean,
                   "dCT", CTvals.dCT,
                   "ddCT", CTvals.ddCT,
                   "ddCT_ctrl_mean", CTvals.ddCT.ctrl.mean,
                   "RE", CTvals.RE,
                   "RE_drop_UND", CTvals.RE.Und
            )
    )
}
