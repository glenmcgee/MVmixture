## function to summarize pairwise clustering
pairwise_clusters <- function(obj){


  dims <- dim(obj$Zbeta)
  matbetay <- matthetay <- matbetax <- matthetax <- matrix(NA,nrow=dims[1]*dims[2],ncol=dims[1]*dims[2])
  for(kk in 1:dims[1]){
    for(jj in 1:dims[2]){
      for(kk2 in 1:dims[1]){
        for(jj2 in 1:dims[2]){
          ## ordering by y
          matbetay[(kk-1)*dims[2]+jj,(kk2-1)*dims[2]+jj2] <- round(100*mean(obj$Zbeta[kk,jj,]==obj$Zbeta[kk2,jj2,]))
          matthetay[(kk-1)*dims[2]+jj,(kk2-1)*dims[2]+jj2] <- round(100*mean(obj$Ztheta[kk,jj,]==obj$Ztheta[kk2,jj2,]))
          ## ordering by x
          matbetax[(jj-1)*dims[1]+kk,(jj2-1)*dims[1]+kk2] <- round(100*mean(obj$Zbeta[kk,jj,]==obj$Zbeta[kk2,jj2,]))
          matthetax[(jj-1)*dims[1]+kk,(jj2-1)*dims[1]+kk2] <- round(100*mean(obj$Ztheta[kk,jj,]==obj$Ztheta[kk2,jj2,]))
        }
      }

    }
  }

  rownames(matbetay) <- colnames(matbetay) <- rownames(matthetay) <- colnames(matthetay) <-
    paste0("y",rep(1:dims[1],each=dims[2]),"x",rep(1:dims[2],dims[1]))
  rownames(matbetax) <- colnames(matbetax) <- rownames(matthetax) <- colnames(matthetax) <-
    paste0("x",rep(1:dims[2],each=dims[1]),"y",rep(1:dims[1],dims[2]))

  return(list(beta_y=matbetay,theta_y=matthetay, ## order by y
              beta_x=matbetax,theta_x=matthetax)) ## order by x
}





library(reshape)
make_heatplot <- function(heat){
  heatplot <- ggplot(data = melt(heat), aes(x=X1, y=X2,fill=value)) +
    geom_tile()+
    scale_fill_gradientn(colours = c("white", "coral4"),values = c(0,1))+
    geom_text(aes(X1, X2, label = value),color = "white", size = 4)+
    labs(x="", y="")+
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = 'top')+
    # guides(x =  guide_axis(angle = 45)) +
    # guides(fill=guide_legend(title="Prob"))+
    guides(fill="none")+
    theme(plot.margin = unit(c(-0.5,0, 0.5, 0), "cm"),
          panel.background = element_blank())+
    ggtitle("")

  return(heatplot)
}









