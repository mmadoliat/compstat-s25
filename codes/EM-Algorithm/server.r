library(shiny)

options(shiny.sanitize.errors = FALSE)

shinyServer(function(input, output, clientData, session) {
  
  dmnorm <- function(x, pars) {
    dmnorm <- 0;
    for (i in 1:length(pars$pi)) dmnorm <- dmnorm + pars$pi[i]*dnorm(x,pars$mu[i],sqrt(pars$sigma2[i]))
    return(dmnorm)
  }
  
  EM <- function(data=data, B=2, eps = 1e-05, iter.max = 100) {
    n <- length(data); k <- ans <- 1;
    par <- list(); par$pi <- par$mu <- par$sigma2 <- list();
    par$pi[[1]] <- par$mu[[1]] <- par$sigma2[[k]] <- rep(1/B,B); 
    if (input$fixms!=TRUE) {
      q <- quantile(data,seq(0,1,1/B));
      for (j in 1:B) {dat <- data[data>=q[j] & data<=q[(j+1)]]; par$mu[[k]][j] <- mean(dat); par$sigma2[[k]][j] <- (length(data)-1)/length(data)*var(dat)}
    } else {
      par$mu[[k]] <- min(data) + (1:B)*(max(data)-min(data))/(B+1); par$sigma2[[k]] <- rep((par$mu[[k]][1]-min(data))/2,B);
    }
    
    p.hat <- matrix(0, nr=n, nc=B);
    while (ans==1 & k < iter.max) {
      for (j in 1:B) {p.hat[,j] <- dnorm(data,mean = par$mu[[k]][j], sd=sqrt(par$sigma2[[k]][j]))}
      p.hat <- p.hat/apply(p.hat,1,sum); k <- k + 1; if (sum(is.nan(p.hat))>0) p.hat[is.nan(p.hat)] <- 0;
      par$pi[[k]] <- apply(p.hat,2,mean);
      if (input$fixms!=TRUE) {
        par$mu[[k]] <- apply(p.hat*data,2,sum)/(n*par$pi[[k]]);
        par$sigma2[[k]] <- apply(p.hat*outer(data,par$mu[[k]],"-")^2,2,sum)/(n*par$pi[[k]])
      } else {
        par$mu[[k]] <- par$mu[[k-1]]; par$sigma2[[k]] <- par$sigma2[[k-1]]
      }
      if (sum(abs(par$pi[[k]]-par$pi[[(k-1)]])) < eps) ans <- 0
    }
    updateNumericInput(session, "modes", max=ifelse(input$fixms==TRUE,1000,10))
    return(par);
  }
  
  output$data <- renderTable({
    inFile <- input$file;    if (is.null(inFile)) return(NULL)
    inData <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
    updateNumericInput(session, "column", max=ncol(inData));
    return(inData)
  })
  
  output$plot = renderPlot({
    inFile <- input$file;    if (is.null(inFile)) return(NULL)
    inData <- read.csv(inFile$datapath, header=input$header, sep=input$sep);
    inData <- as.vector(na.exclude(inData[,input$column]));
    bins <- seq(min(inData), max(inData), length.out = input$bins + 1)
    result <- EM(data=inData, B=input$modes); K <- length(result$pi);
    rData <- range(hist(inData, breaks = bins, probability = TRUE, plot = FALSE)$density)*1.25
    hist(inData, breaks = bins, col = 'darkgray', border = 'white', probability = TRUE, xlab="Data", main = "Histogram and GMM fits", ylim=c(0,max(rData)))
    if ("init" %in% input$show) {
      pars <- list(pi=result$pi[[1]],mu=result$mu[[1]],sigma2=result$sigma2[[1]]);
      curve(dmnorm(x,pars),add=TRUE,col=4,lwd=2,lty=2)
    }
    if ("final" %in% input$show) {
      pars <- list(pi=result$pi[[K]],mu=result$mu[[K]],sigma2=result$sigma2[[K]])
      curve(dmnorm(x,pars),add=TRUE,col=2,lwd=2)
    }
      k <- input$step; if (k==0 || k>K) {k <- K}; updateNumericInput(session, "step", value = k, max=K)
      pars <- list(pi=result$pi[[k]],mu=result$mu[[k]],sigma2=result$sigma2[[k]])
      curve(dmnorm(x,pars),add=TRUE,lwd=1.5)
    if ("legend" %in% input$show) {
      legend("topleft",c(paste("step",k),"initial","final"),lty=c(1,2,1),lwd=c(1.5,2,2),col=c(1,4,2))
    }
  })
  
  output$ic = renderPlot({ 
    inFile <- input$file;    if (is.null(inFile)) return(NULL); m <- ifelse(input$fixms==TRUE,100,10)
    inData <- read.csv(inFile$datapath, header=input$header, sep=input$sep);
    inData <- as.vector(na.exclude(inData[,input$column]));
    result <- pars <- list(); K <- logL <- AIC <- BIC <- NULL; n <- length(inData)
    for (B in 1:m) {
      result[[B]] <- EM(data=inData, B=B); K[B] <- length(result[[B]]$pi);
        pars[[B]] <- list(pi=result[[B]]$pi[[K[B]]],mu=result[[B]]$mu[[K[B]]],sigma2=result[[B]]$sigma2[[K[B]]])
        logL[B] <- sum(log(dmnorm(inData,pars[[B]]))); AIC[B] <- 2*3*B-2*logL[B]; BIC[B] <- 3*B*log(n)-2*logL[B]
    }
    plot(AIC, type="b", log="y", axes=FALSE, xlab="# of Modes", ylab="AIC/BIC", main = "Information Criterion", ylim=range(c(AIC,BIC),na.rm=TRUE))
    points(BIC, col=2, pch=20, type="b"); vAIC <- which(AIC==min(AIC,na.rm=TRUE)); vBIC <- which(BIC==min(BIC,na.rm=TRUE))
    abline(v=vAIC,lty=2); abline(v=vBIC,lty=3,col=2); 
    axis(1,at=c(1,vAIC,vBIC,m),labels=c(as.character(c(1,vAIC,vBIC,m)))); axis(2)
    legend("bottomright",c("AIC","BIC"),lty=c(2,3),col=c(1,2),pch=c(1,20))
    updateNumericInput(session, "modes", value=vAIC)
  })

  output$summary <- renderTable({
    inFile <- input$file;    if (is.null(inFile)) return(NULL)
    inData <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
    inData <- as.vector(na.exclude(inData[,input$column]));
    bins <- seq(min(inData), max(inData), length.out = input$bins + 1)
    result <- EM(data=inData, B=input$modes);
    tab <- cbind(t(matrix(unlist(result$pi),nr=input$modes)),t(matrix(unlist(result$mu),nr=input$modes)),t(matrix(unlist(result$sigma2),nr=input$modes)))
    colnames(tab) <- paste(rep(c("pi","mu","s2"),each=input$modes),rep(1:input$modes,3),sep=".")
    ifelse(input$fixms==TRUE,return(data.frame(tab[2,1:(ncol(tab)/3)])),return(data.frame(tab[,1:min(ncol(tab),100)])))
  })
})