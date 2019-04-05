#' Plotting DVDtest
#' 
#' Plot a list of the DVDtest-related figures via \code{ggplot2}
#' 
#' Fig \code{pfig} illustrates a collection of the p value curves among all \code{k}.
#' Fig \code{kfig} illustrates a collection list of the figures with varying distributions among all 
#' \code{k}, highlighted the zones of small p values in dark, where \code{ydata1} and 
#' \code{ydata2} mark as red and blue, repectively. The dashed and dotted lines denote the smooth mean function and +/- 2 
#' sigma, respectively.
#' 
#' 
#' 
#' @param tobj a return test object of \code{\link{DVDtest}}
#' @param kxlab a title for the x axis, \code{.index}
#' @param kylab a title for the y axis, \code{.value}
#' @param kname a name for \code{k}, e.g. ROI in the references
#' @param thhd a number to mark out the zone of the p-value of interest. Defaults to 0.05
#' 
#' @return a list of ggplot objects on p value curves and varying distributions
#' @author Meng Xu
#' @note Please contact the maintainer if need more details.
#' @references reiss-EMR18.pdf
#' @keywords ggplot
#' @export
#' @import reshape2 ggplot2 mgcv gamlss gamlss.dist
#' @seealso Examples in \code{\link{DVDtest}}



DVDplot <-
  function(tobj, kxlab = NULL, kylab = NULL, kname= NULL,thhd=0.05){
    kfig <- list()
    p.dft<-as.data.frame(tobj$pval)
    p.dft$.index<-tobj$.index
    p.df<-melt(p.dft, id.var=".index")
    pfig <- ggplot(p.df, aes_(x=~.index,y=~value,group=~variable,colour=~variable))+
      geom_line(show.legend = FALSE)+
      geom_hline(yintercept = 0.05, linetype = 2, col = "red")+
      labs(y = "p value", x = kxlab)
    
    
    
    for (i in 1:NCOL(tobj$pval)) {
      for (k in 1:length(tobj$pval[,i])) if (tobj$pval[k,i]>thhd) tobj$pval[k,i]=1
      pl1<-data.frame(.index = tobj$.index, mn=tobj$vdparam[[i]]$pred1$mu,
                      std=tobj$vdparam[[i]]$pred1$sigma,
                      pval=tobj$pval[,i])
      pl2<-data.frame(.index = tobj$.index, mn=tobj$vdparam[[i]]$pred2$mu,
                      std=tobj$vdparam[[i]]$pred2$sigma,
                      pval=tobj$pval[,i])
      kfig[[i]]=ggplot()+geom_line(data=pl1,aes_(x=~.index,y = ~mn,alpha=~pval),
                                   size=2,color="red")+
        geom_line(data=pl1,aes_(x=~.index,y = ~mn),size=0.2,color="red",linetype=2)+
        geom_line(data=pl1,aes_(x=~.index,y =~(mn-1.96*std),alpha=~pval),size=0.5,color="red")+
        geom_line(data=pl1,aes_(x=~.index,y =~(mn-1.96*std)),size=0.2,color="red",linetype=3)+
        geom_line(data=pl1,aes_(x=~.index,y =~(mn+1.96*std),alpha=~pval),size=0.5,color="red")+
        geom_line(data=pl1,aes_(x=~.index,y =~(mn+1.96*std)),size=0.2,color="red",linetype=3)+
        geom_ribbon(data=pl1,aes_(x=~.index,ymin=~(mn-1.96*std),ymax=~(mn+1.96*std)),alpha=0.1,
                    fill="red")+
        geom_line(data=pl2,aes_(x=~.index,y = ~mn,alpha=~pval),size=2,color="blue")+
        geom_line(data=pl2,aes_(x=~.index,y = ~mn),size=0.2,color="blue",linetype=2)+
        geom_line(data=pl2,aes_(x=~.index,y = ~(mn-1.96*std),alpha=~pval),size=0.5,color="blue")+
        geom_line(data=pl2,aes_(x=~.index,y =~(mn-1.96*std)),size=0.2,color="blue",linetype=3)+
        geom_line(data=pl2,aes_(x=~.index,y = ~(mn+1.96*std),alpha=~pval),size=0.5,color="blue")+
        geom_line(data=pl2,aes_(x=~.index,y =~(mn+1.96*std)),size=0.2,color="blue",linetype=3)+
        geom_ribbon(data=pl2,aes_(x=~.index,ymin=~(mn-1.96*std),ymax=~(mn+1.96*std)),alpha=0.1,
                    fill="blue")+
        labs(title = paste0(kname,i),y=kylab,x=kxlab)+
        scale_alpha(range = c(.8, 0))+
        theme(legend.position="non")
    }
    return(list(pfig = pfig, kfig=kfig))
  }

