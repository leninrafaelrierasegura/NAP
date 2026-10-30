# Matched-physics comparisons: only the coefficient family/order differs.
rough_compare_coefficients <- function(h_results,tau_results,
    baseline_h=readRDS("data_files/rough_comparison/brasil_h_results.rds"),
    baseline_tau=readRDS("data_files/rough_comparison/brasil_tau_results.rds")) {
  rows_h <- rows_tau <- list()
  for (i in seq_along(h_results)) {
    r <- h_results[[i]]
    for (field in c("J","alpha","time_steps","T_final","diffusion","kappa","fit_J"))
      stopifnot(identical(r$parameters[[field]],baseline_h$parameters[[field]]))
    stopifnot(identical(dim(r$errors),dim(baseline_h$errors)))
    n <- nrow(r$errors)
    rows_h[[i]] <- data.frame(alpha=r$parameters$alpha,order=r$parameters$rational_order,
      brasil_raw=baseline_h$rates$raw,operator_raw=r$rates$raw,
      brasil_corrected=baseline_h$rates$corrected,operator_corrected=r$rates$corrected,
      brasil_finest_error=baseline_h$errors[n,],operator_finest_error=r$errors[n,],
      finest_error_ratio=r$errors[n,]/baseline_h$errors[n,],
      max_error_change_percent=apply(abs(r$errors/baseline_h$errors-1)*100,2L,max))
  }
  for (i in seq_along(tau_results)) {
    r <- tau_results[[i]]
    for (field in c("J","alpha","N","T_final","diffusion","kappa","fit_indices"))
      stopifnot(identical(r$parameters[[field]],baseline_tau$parameters[[field]]))
    stopifnot(identical(dim(r$errors),dim(baseline_tau$errors)))
    n <- nrow(r$errors)
    rows_tau[[i]] <- data.frame(alpha=r$parameters$alpha,order=r$parameters$rational_orders,
      brasil_slope=baseline_tau$rates$observed,operator_slope=r$rates$observed,
      brasil_finest_error=baseline_tau$errors[n,],operator_finest_error=r$errors[n,],
      finest_error_ratio=r$errors[n,]/baseline_tau$errors[n,],
      max_error_change_percent=apply(abs(r$errors/baseline_tau$errors-1)*100,2L,max))
  }
  list(h=do.call(rbind,rows_h),tau=do.call(rbind,rows_tau),
    baseline_h=baseline_h,baseline_tau=baseline_tau,h_results=h_results,tau_results=tau_results)
}

rough_coefficient_comparison_plot <- function(comparison,study=c("h","tau")) {
  study <- match.arg(study)
  baseline <- if(study=="h")comparison$baseline_h else comparison$baseline_tau
  candidates <- if(study=="h")comparison$h_results else comparison$tau_results
  frame <- function(r,label) {
    x <- if(study=="h")2^-r$parameters$J else r$tau
    d <- expand.grid(index=seq_along(x),column=seq_along(r$parameters$alpha))
    d$x <- x[d$index]
    d$error <- as.vector(r$errors)
    d$alpha <- factor(r$parameters$alpha[d$column],levels=r$parameters$alpha)
    d$method <- label
    d
  }
  data <- frame(baseline,"BRASIL baseline")
  for(r in candidates) {
    m <- if(study=="h")r$parameters$rational_order else r$parameters$rational_orders[1L]
    data <- rbind(data,frame(r,paste0("Operator m=",m)))
  }
  data$method <- factor(data$method,levels=c("BRASIL baseline","Operator m=4","Operator m=8"))
  ggplot2::ggplot(data,ggplot2::aes(x,error,color=method,linetype=method,group=method))+
    ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(size=1.7)+
    ggplot2::scale_x_log10()+ggplot2::scale_y_log10()+
    ggplot2::scale_color_manual(values=c("#0072B2","#D55E00","#009E73"))+
    ggplot2::scale_linetype_manual(values=c("solid","dashed","dotdash"))+
    ggplot2::facet_wrap(~alpha,ncol=2L,scales="free_y",labeller=ggplot2::label_both)+
    ggplot2::labs(x=if(study=="h")expression(h) else expression(tau),
      y=expression(E),color=NULL,linetype=NULL,
      title=if(study=="h")"Spatial error: matched coefficient comparison" else
        "Temporal error: matched coefficient comparison")+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(legend.position="bottom",panel.grid.minor=ggplot2::element_blank())
}

rough_save_comparison <- function(comparison,output_dir="data_files/rough_comparison") {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(comparison,file.path(output_dir,"comparison.rds"))
  utils::write.csv(comparison$h,file.path(output_dir,"spatial_comparison.csv"),row.names=FALSE)
  utils::write.csv(comparison$tau,file.path(output_dir,"temporal_comparison.csv"),row.names=FALSE)
  for(study in c("h","tau"))for(extension in c("png","pdf"))
    ggplot2::ggsave(file.path(output_dir,paste0(study,"_coefficient_comparison.",extension)),
      rough_coefficient_comparison_plot(comparison,study),width=8,height=10,dpi=300,bg="white")
  invisible(comparison)
}
