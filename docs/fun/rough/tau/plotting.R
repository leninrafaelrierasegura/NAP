rough_tau_convergence_plot <- function(result) {
  p <- result$parameters
  errors <- result$errors
  fit <- p$fit_indices
  constants <- apply(errors[fit,,drop=FALSE]/result$tau[fit],2L,function(x)exp(mean(log(x))))
  frame <- expand.grid(index=seq_along(result$tau),column=seq_along(p$alpha))
  frame$tau <- result$tau[frame$index]
  frame$error <- as.vector(errors)
  frame$alpha <- factor(p$alpha[frame$column],levels=p$alpha)
  guide <- frame
  guide$error <- guide$tau*constants[guide$column]
  ggplot2::ggplot(frame,ggplot2::aes(tau,error,color=alpha,group=alpha))+
    ggplot2::geom_line(data=guide,linetype="dashed",linewidth=.65,alpha=.75)+
    ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(size=2.1)+
    ggplot2::scale_x_log10(breaks=result$tau,labels=function(x)paste0("1/",round(1/x)))+
    ggplot2::scale_y_log10()+
    ggplot2::scale_color_viridis_d(option="D",begin=.05,end=.85)+
    ggplot2::facet_wrap(~alpha,ncol=2L,labeller=ggplot2::label_both)+
    ggplot2::labs(x=expression(tau~"(log scale)"),
      y=expression(E[tau]),color=expression(alpha),
      title="Temporal convergence with rough forcing")+
    ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(legend.position="none",panel.grid.minor=ggplot2::element_blank(),
      plot.title=ggplot2::element_text(size=13,face="bold"))
}

rough_tau_save_figures <- function(result,output_dir="data_files/rough_tau") {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    plot <- rough_tau_convergence_plot(result)
    name <- "conv_tau_rough"
    for(extension in c("png","pdf"))ggplot2::ggsave(
      file.path(output_dir,paste0(name,".",extension)),plot,
      width=8,height=9,dpi=300,bg="white")
  invisible(output_dir)
}
