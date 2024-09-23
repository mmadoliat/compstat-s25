library(fda); library(mgcv); library(dendextend); library(svd); library(shiny); 
source("fpca.r"); 

ui<-fluidPage(tags$head(tags$style(HTML("body { max-width: 1250px !important; }"))),
  titlePanel("FDA Workshop"),
  sidebarLayout(
    sidebarPanel(width=3, tags$head(tags$style(type="text/css", ".well { max-width: 300px; }")),
      sliderInput("xdeg", HTML("Degree of B-spline Basis:"), min = 0, max = 5, value = 3, step = 1, width="250px"),
      uiOutput("xdf", width="250px"), 
      fluidRow(column(6, checkboxInput('twoD', 'Two-Dim.', FALSE, width="125px")), column(6, checkboxInput('diffP', 'Diff. Pen.', FALSE, width="125px"))),
      uiOutput("ydeg", width="250px"), uiOutput("ydf", width="250px"), tags$hr(style="border-color: red;", width="150px"),
      sliderInput("dimn", HTML("# of Common Basis <br/>(# of Clusters for Clustering):"), min = 1, max = 10, value = 3, step = 1, width="210px"),
      uiOutput("plot", width="240px"), uiOutput("p.ang1", width="240px"), uiOutput("p.ang2", width="240px")),
    mainPanel(width=9, tags$style(type="text/css", ".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; },"),#".nav-tabs {font-size: 10px}"),
      tabsetPanel(id="Panel", type = "tabs", 
          tabPanel(title="Data", value="Data",
             column(12,uiOutput("ts.selected", align = "center"), style="color:red;"),
             fluidRow(column(4,radioButtons('f.choice', 'Choose from:', c("Server" = "server", "Upload" = "upload", "Simulate" = "sim"), selected= "sim", inline = TRUE, width="250px")), 
                      column(4,uiOutput("s.choice", width="250px")), column(4,uiOutput('sts.twoD', width="250px"))),
             column(4,uiOutput("file")), column(8,uiOutput("sep"),uiOutput("header")),
             column(4,uiOutput("model")), column(4,uiOutput("m.rep")), column(2,uiOutput("n.sd")), column(2,uiOutput("c.sd")),
             column(8,plotOutput("data.plot", height = 600, width = 600)), column(4,tableOutput('data'))),
          tabPanel("Basis Functions",
                   column(8,plotOutput("basis.desc", height = 600, width = 600)), column(4,uiOutput("basis.n", width="300px"))),
          tabPanel("Data Description (FDA Summary)", 
             column(4,uiOutput("desc", width="250px")), column(4,uiOutput("ams.choice", width="400px"), uiOutput("run.fda.gcv", width="200px")), column(4,uiOutput("sts.choice")),
             fluidRow(column(8,plotOutput("data.desc", height = 600, width = 600)), 
                      column(4,uiOutput("s.plot"), fluidRow(column(8,uiOutput("b.indx")), column(4,uiOutput("s.CI"))), column(4,uiOutput("diff")), column(8,uiOutput("tuning.fda"))))),
          tabPanel("FPCA",
             column(4,uiOutput("runit", width="200px"),tags$hr(style="border-color: red;")), column(4,uiOutput("ams.FPCA", width="400px")), column(4,uiOutput("sts.FPCA")), 
             fluidRow(column(8, plotOutput("FPCA.out", height = 600, width = 600)), 
                      column(4, uiOutput("results"), uiOutput("sn.plot"), fluidRow(column(6,uiOutput("bn.indx")), column(6,uiOutput("tuning.fpca"))),
                             uiOutput("so.plot"), uiOutput("diffo"), column(6,uiOutput("bo.indx")), column(6,uiOutput("tuningo.fda"))))),
          tabPanel("Manual", includeMarkdown("report.Rmd"))))
  )
)

server<-function(input, output, clientData, session) {

  load("data/servshiny.Rda"); ind.lam <- list(); ind.lam$fsa <- ind.lam$fpca <- ind.lam$D <- temp <- B <- P <- Ts <- Trs <- NULL; 
  df1 <- 100; df2 <- 15; vf1 <- 20; vf2 <- 10; lambdas <- lambdas1 <- c(-75,25); lambdas2 <- c(-25,25); T1 <- 100; T2 <- 35;

  simulate <- function() {if (memory.size()>700) gc();
    Trs <- NULL; T <- ifelse(input$twoD, T2, T1); set.seed(T*input$m.rep*input$c.sd*input$n.sd)
    t <- seq(0,1,length.out=T); f1 <- t+sin(2*pi*t); f2 <- cos(6*pi*t); 
    a <- rnorm(input$m.rep,1,sd=input$c.sd); b <- rnorm(input$m.rep,1,sd=input$c.sd); 
    c <- rnorm(input$m.rep,1,sd=input$c.sd); d <- rnorm(input$m.rep,1,sd=input$c.sd)
    if (!input$twoD) {
      if ("f1" %in% input$model) {Trs <- matrix(a%x%f1,nr=T)}
      if ("f2" %in% input$model) {Trs <- cbind(Trs,matrix(b%x%f2,nr=T))}
      if ("f12" %in% input$model) {Trs <- cbind(Trs,matrix(c%x%f1+d%x%f2,nr=T))}
    } else {
      if ("f1" %in% input$model) {Trs <- matrix(a%x%f1%x%(f1),nr=T^2)}
      if ("f2" %in% input$model) {Trs <- cbind(Trs,matrix(b%x%f2%x%(f2),nr=T^2))}
      if ("f12" %in% input$model) {Trs <- cbind(Trs,matrix(c%x%f1%x%f2,nr=T^2))}
      if ("f21" %in% input$model) {Trs <- cbind(Trs,matrix(d%x%f2%x%f1,nr=T^2))}
    }
    Trs <- scale(Trs,scale=F)
    return(list(Trs=Trs, noise=rnorm(input$m.rep*ifelse(input$twoD,T^2,T)*length(input$model),sd=input$n.sd)))
  }

  observeEvent(input$twoD, {
    updateTabsetPanel(session, "Panel", selected = "Data")
  })
    
  output$xdf <- renderUI({
    sliderInput("xdf", HTML("Deg. of freedom of B-spline Basis:"), min = input$xdeg+1, max = ifelse(input$twoD,df2,df1), value = ifelse(input$twoD,vf2,vf1), step = 1, width="250px")
  })

  output$ydeg <- renderUI({
    if (!input$twoD) return();
    sliderInput("ydeg", HTML("Degree of B-spline Basis:"), min = 0, max = 5, value = 3, step = 1, width="250px")
  })
  
  output$ydf <- renderUI({
    if (!input$twoD) return()
    sliderInput("ydf", HTML("Deg. of freedom of B-spline Basis:"), min = input$ydeg+1, max = ifelse(input$twoD,df2,df1), value = ifelse(input$twoD,vf2,vf1), step = 1, width="250px")
  })
  
  output$plot <- renderUI({
    if (!input$twoD) return();
    selectInput("plot","3D Visualition: ", choices = c("Contour Plot"="contour", "Heat Map"="hmap", "Perspective"="persp"), selected="contour", width="240px")
  })
  
  output$p.ang1 <- renderUI({
    if (input$plot!="persp" || !input$twoD) return();
    sliderInput("p.ang1", "Azimuthal Angle:", min = -180, max = 180, value = 40, step = 1, width="240px")
  })
  
  output$p.ang2 <- renderUI({
    if (input$plot!="persp" || !input$twoD) return();
    sliderInput("p.ang2", "Colatitude Angle:", min = -180, max = 180, value = 20, step = 1, width="240px")
  })
    
  output$s.choice <- renderUI({
    if (input$f.choice!="server") return(); s.choices <- 1:length(Xs); name <- NULL;
    for (i in s.choices) name[i] <- paste0("TS-m=",length(Xs[[i]]),"-n=",length(Xs[[i]][[1]])); names(s.choices) <- names(Xs) <<- name; 
    selectInput("s.choice","Select a file from server: ", choices = s.choices, width="250px"); 
  })
  
  output$sts.twoD <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || !input$twoD) return();
    if (input$f.choice=="upload") {Ts <- as.matrix(read.table(input$file$datapath, header=input$header, sep=input$sep))
    } else if (input$f.choice=="server") {i <- as.numeric(input$s.choice); Ts <- matrix(unlist(Xs[[i]]),nr=length(Xs[[i]][[1]]))}
    if (input$f.choice=="sim") {Max <- length(input$model)*input$m.rep} else {Max <- ncol(Ts)}
    sliderInput("sts.twoD","2D Function #:", min = 1, max = Max, value = 1, step=1, width="250px")
  })

  output$model <- renderUI({
    if (input$f.choice!="sim") return();
    if (!input$twoD) {choices <- c("a.f(t)" = "f1", "b.g(t)" = "f2", "c.f(t) + d.g(t)" = "f12")} else {choices <- c("a.f1" = "f1", "b.f2" = "f2", "c.f12" = "f12", "d.f21" = "f21")}
    checkboxGroupInput('model', 'Model:', choices=choices, selected="f1", inline = TRUE, width="250px")
  })

  output$m.rep <- renderUI({
    if (input$f.choice!="sim") return();
    sliderInput("m.rep","Replicates of each", min = 1, max = 20, value = 1, width="250px")
  })
  
  output$n.sd <- renderUI({
    if (input$f.choice!="sim") return(); 
    sliderInput("n.sd","Noise SD:", min=0, max = 1, value = 0.2, width="125px");
  })
  
  output$c.sd <- renderUI({
    if (input$f.choice!="sim") return(); 
    sliderInput("c.sd","Coeff. SD:", min=0, max = 1, value = 0.2, width="125px");
  })
  
  output$file <- renderUI({
    if (input$f.choice!="upload") return();
    fileInput('file', 'Choose CSV File', accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
  })
  output$sep <- renderUI({if (input$f.choice!="upload") return(); radioButtons('sep', 'Separator', c(","=',', ":"=':', ";"=';', Tab='\t'), ',', inline = TRUE)})
  output$header <- renderUI({if (input$f.choice!="upload") return(); checkboxInput('header', 'Header', TRUE)})
  
  output$ts.selected = renderText({
    if (input$f.choice=="upload" && is.null(input$file)) return("<b>Select a 'csv' file that contain the time series in its columns</b>")
    if (input$f.choice=="upload") {Ts <- as.matrix(read.table(input$file$datapath, header=input$header, sep=input$sep))
    } else if (input$f.choice=="server") {i <- as.numeric(input$s.choice); Ts <- matrix(unlist(Xs[[i]]),nr=length(Xs[[i]][[1]]))
    } else if (input$f.choice=="sim") {Ts <- matrix(0,nrow=ifelse(input$twoD,T2^2,T1),ncol=length(input$model)*input$m.rep)}
    text <- paste("<b>",ncol(Ts),"Time series of length",nrow(Ts),"</b>"); return(text)
  })

  output$data <- renderTable({if (memory.size()>700) gc();
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="upload") {Ts <<- as.matrix(read.table(input$file$datapath, header=input$header, sep=input$sep));
    } else if (input$f.choice=="server") {i <- as.numeric(input$s.choice); Ts <<- matrix(unlist(Xs[[i]]),nr=length(Xs[[i]][[1]]))
    } else {simul <- simulate(); Trs <<- simul$Trs; Ts <<- Trs + simul$noise}; if (is.null(colnames(Ts))) {colnames(Ts) <<- paste("TS:",1:ncol(Ts))}; Ts <<- scale(Ts,scale=F);
    updateSelectInput(session, "desc", selected = "gcv"); updateSliderInput(session, "dimn", max=min(10,ncol(Ts)), value=min(2,ncol(Ts)));
    return(head(as.matrix(Ts[,1:min(9,ncol(Ts))]),15))
  })

  output$data.plot = renderPlot({if (memory.size()>700) gc();
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="server") {fname <- names(Xs)[as.numeric(input$s.choice)]; i <- as.numeric(input$s.choice); Ts <- matrix(unlist(Xs[[i]]),nr=length(Xs[[i]][[1]]))} 
    else if (input$f.choice=="upload") {fname <- input$file$name; Ts <- as.matrix(read.table(input$file$datapath, header=input$header, sep=input$sep))} 
    else {fname <- "Simulation"; simul <- simulate(); Ts <- simul$Trs+simul$noise}; Ts <- scale(Ts,scale=F);
    if (!input$twoD) {
      ts.plot(Ts, main=paste("Time Series -", fname), ylab="", ylim=range(Ts), gpars=list(xaxt="n"))
      if (input$f.choice=="sim") for (i in 1:ncol(simul$Trs)) points(simul$Trs[,i],type="l",col=2)
    } else {
      if (input$plot=="contour") contour(matrix(Ts[,input$sts.twoD],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
      else if (input$plot=="hmap") filled.contour(matrix(Ts[,input$sts.twoD],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
      else if (input$plot=="persp") persp(matrix(Ts[,input$sts.twoD],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
	  title(main=paste("2D Function -",fname,"- #",input$sts.twoD))
    }
  })
  
  output$basis.n <- renderUI({
    sliderInput("basis.n", "Basis #:", min = 1, max = ifelse(!input$twoD,input$xdf,input$xdf*input$ydf), value = 1, step=1, width="400px")
  })
  
  output$basis.desc = renderPlot({if (memory.size()>700) gc();
    if (!input$twoD) {
      xs <- seq(0,1,length.out=1000); Bx <- fda::bsplineS(xs,breaks=seq(0,1,length.out=input$xdf-input$xdeg+1),norder=input$xdeg+1) 
      ts.plot(Bx, col=8, main="Bspline Basis", xlab="Grid Points", gpars=list(xaxt="n"))
      points(Bx[,input$basis.n], type="l", lwd=2, col=2)
	    axis(1,trunc(summary(1:nrow(Bx))[-4]))
	  } else {
	    a <- 5; b <- 100; xs <- seq(0,a,length.out = b); ys <- seq(0,a,length.out = b);
	    Bx <- fda::bsplineS(xs,breaks=seq(0,a,length.out=input$xdf-input$xdeg+1),norder=input$xdeg+1)
	    By <- fda::bsplineS(ys,breaks=seq(0,a,length.out=input$ydf-input$ydeg+1),norder=input$ydeg+1) 
	    grid.out <- expand.grid(xs,ys); B <- By%x%Bx; x.range <- range(xs); y.range <- range(ys);
      if (input$plot=="contour") {
        plot(0,col=0,ylim=c(ys[1]-max(By),ys[length(ys)]),xlim=c(xs[1]-max(Bx),xs[length(xs)]), main="2D Basis Function", xaxt="n", yaxt="n",
             xlab="x", ylab="y", cex.lab=2, frame.plot=F); polygon(c(0,a,a,0),c(0,0,a,a), border=4, lwd=2,col=7);
          axis(1,summary(xs)[-4],round(summary(xs)[-4]/a,2)); axis(2,summary(ys)[-4],round(summary(ys)[-4]/a,2))
        for (i in 1:input$xdf) points(xs, Bx[,i]+y.range[1]-max(Bx), col=3, type="l")
        for (j in 1:input$ydf) points(By[,j]+x.range[1]-max(By), ys, type="l", col=3)
        contour(xs,ys,matrix(B[,input$basis.n],nr=b),add=T,col=2,lwd=2);
          j <- (input$basis.n %/% input$xdf)+1; i <- (input$basis.n %% input$xdf); if (i==0) {i <- input$xdf; j <- j-1};
          points(xs, Bx[,i]+y.range[1]-max(Bx), col=2, type="l", lwd=2)
          points(By[,j]+x.range[1]-max(By), ys, col=2, type="l", lwd=2)
      } else if (input$plot=="hmap") {
        filled.contour(xs,ys,matrix(B[,input$basis.n],nr=b), main="2D Basis Function", 
                       xlab="x", ylab="y", xlim=x.range, ylim=y.range, cex.lab=2)
      } else if (input$plot=="persp") {
        persp(xs, ys, matrix(B[,input$basis.n],nr=b), main="2D Basis Function", cex.lab=2,
              xlab="x", ylab="y", zlab="Basis", xlim=x.range, ylim=y.range, theta=input$p.ang1, phi=input$p.ang2)
      }
	  }
  })
  
  output$desc <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    choices <- list(Summary=c("Time Series" = "ts"), More=c("Tuning Parameter (GCV)"="gcv", "Scree Plot" = "scree")); 
    if (input$twoD) {choices[[1]] <- c("2D Functions" = "2f")}
    selectInput("desc","Select",choices=choices, width="250px")
  })

  output$ams.choice <- renderUI({ 
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$desc!="ts") return();
    radioButtons('ams.choice', 'Plot Choices:', c("All" = "all", "Multiple" = "multiple", "Single" = "single"), selected= "all", inline = TRUE, width="400px")
  })
  
  output$sts.choice <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (((input$desc=="ts" && input$ams.choice!="all") || input$desc=="2f") && !(input$s.plot=="bf" && length(input$s.plot)==1) && length(input$s.plot))
      if (input$ams.choice=="single" || input$desc=="2f") {sliderInput("sts.choice", "Choose function:", min = 1, max = ncol(Ts), value = ifelse(is.null(input$sts.choice),1,input$sts.choice), step = 1, width="400px")} else {
        sliderInput("sts.choice", "Choose a range( '< 10' ):", min = 1, max = ncol(Ts), value = c(1,min(10,ncol(Ts))), dragRange=TRUE, step = 1, width="400px")
      }
  })
  
  output$s.plot <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || !input$desc%in%c("ts","2f")) return();
    choices <- c("Time Series" = "ts", "True Functions" = "tf", "Basis Func." = "bf", "B-Spline Smoothing" = "bss", "Pen. B-Sp. Smooth." = "pbss")
    if (input$f.choice!="sim") {choices <- choices[-2]}; if (input$desc=="2f") choices <- c("2D Functions" = "2f", choices[-1]); 
    if (input$ams.choice=="multiple" || input$desc=="2f") {radioButtons('s.plot', 'Plot:', choices=choices, selected=choices[1], width="250px")} 
    else {checkboxGroupInput('s.plot', 'Plot:', choices=choices, selected=choices[1], width="250px")}
  })

  output$b.indx <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || !input$desc%in%c("ts","2f") || !length(intersect(input$s.plot,c("bf","bss","pbss")))) return();
    val <- c(1,ifelse(input$desc!="2f",input$xdf,input$xdf*input$ydf)); if(length(input$b.indx)==2) val <- input$b.indx
    sliderInput("b.indx", "Basis Contr.", min = 1, max = ifelse(input$desc!="2f",input$xdf,input$xdf*input$ydf), value = val, dragRange=TRUE, step = 1, width="200px")
  })
  
  output$s.CI <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$desc=="ts" && length(intersect(input$s.plot,c("bss","pbss")))) 
      if (input$ams.choice=="single") checkboxInput('s.CI', 'Show CI', ifelse(is.null(input$s.CI), FALSE, input$s.CI), width="200px")
  })
  
  output$diff <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if ((input$desc%in%c("ts","2f") && !"pbss"%in%input$s.plot) || input$desc=="scree" || !input$diffP) return()
    sliderInput('diff','Pen.Dif.Ord.', value = 2, min = 0, max=min(4,ifelse(input$twoD,min(input$xdf,input$ydf),input$xdf)-1), step=1, width="200px")
  })
  
  output$tuning.fda <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if ((input$desc=="gcv" && !(input$run.fda.gcv)) || (input$desc%in%c("ts","2f") && !"pbss"%in%input$s.plot) || input$desc=="scree") return()
    if (input$twoD) lambdas <<- lambdas2 else lambdas <<- lambdas1;
    sliderInput("tuning.fda", "Tun. Par. (2^i):", min = lambdas[1], max = lambdas[2], value = ifelse(is.null(ind.lam$fda),0,ind.lam$fda), width="200px")
  })
  
  output$run.fda.gcv <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$desc!="gcv") return();
    actionButton('run.fda.gcv', paste('update GCV'))
  })
  
  fda.gcv <- eventReactive(input$run.fda.gcv, {if (memory.size()>700) gc();
    if (input$twoD) lambdas <<- lambdas2 else lambdas <<- lambdas1;
    withProgress(message = 'FDA.GCV: Tuning Parameter', value = 0, { DF <- GCV <- NULL; lams <- 2^(lambdas[1]:lambdas[2])
      for (l in 1:length(lams)) {
        S <- B%*%solve(t(B)%*%B + lams[l]*P)%*%t(B); DF <- c(DF, sum(diag(S)));
        GCV <- c(GCV, sum((Ts-S%*%Ts)^2)/sum(diag(diag(1,ncol(S))-S))^2);
        incProgress(1/length(lambdas[1]:lambdas[2]), detail = paste0("2^", log2(lams[l])));
      }; temp <<- 1; return(list(GCV=GCV, DF=DF, lambdas=lams))
    })  
  })
    
  output$data.desc <- renderPlot({if (memory.size()>700) gc();
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="server") {fname <- names(Xs)[as.numeric(input$s.choice)]} else if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- "Simulation"}
    indx <- as.numeric(input$sts.choice);
    if (length(intersect(input$s.plot,c("bf","bss","pbss")))) b.indx <- input$b.indx[1]:input$b.indx[2]
    if (length(intersect(input$s.plot,c("bf","bss","pbss"))) || input$desc=="gcv") {
      B <<- fda::bsplineS(seq(0,1,length.out = ifelse(input$twoD,sqrt(nrow(Ts)),nrow(Ts))), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1);
      if ("pbss" %in% input$s.plot || input$desc=="gcv") {
        if (input$diffP) {if (input$diff!=0) {L <- diff(diag(input$xdf), differences=input$diff)} else {L <- diag(input$xdf)}; P <<- t(L)%*%L
        } else {P <<- fda::bsplinepen(fda::create.bspline.basis(nbasis=input$xdf,norder=input$xdeg+1))};
      }
      if (input$twoD) {
        B <<- B%x%fda::bsplineS(seq(0,1,length.out = sqrt(nrow(Ts))), breaks=seq(0,1,length.out=input$ydf-input$ydeg+1), norder=input$ydeg+1);
        if ("pbss" %in% input$s.plot || input$desc=="gcv")
          if (input$diffP) {if (input$diff!=0) {L <- diff(diag(input$ydf), differences=input$diff)} else {L <- diag(input$ydf)}; P <<- P%x%diag(ncol(L)) + diag(nrow(P))%x%(t(L)%*%L)
          } else {P2 <- fda::bsplinepen(fda::create.bspline.basis(nbasis=input$ydf,norder=input$ydeg+1)); P <<- P%x%diag(nrow(P2)) + diag(nrow(P))%x%P2};
      }
	    if ("bss" %in% input$s.plot) {
	      cB <- solve(t(B)%*%B)%*%t(B); 
	      if (input$desc=="ts") {S <- B%*%cB; vB <- ifelse(nrow(Ts)>ncol(B), sum((Ts-S%*%Ts)^2)/((nrow(Ts)-ncol(B))*ncol(Ts))*diag(S), 0)}
	    }
      if ("pbss" %in% input$s.plot) {
        cP <- solve(t(B)%*%B+2^(input$tuning.fda)*P)%*%t(B);
        if (input$desc=="ts") {S <- B%*%cP; vP <- ifelse(nrow(Ts)>ncol(B), sum((Ts-S%*%Ts)^2)/((nrow(Ts)-ncol(B))*ncol(Ts))*diag(S%*%t(S)), 0)}
      }
    }
    if ("bf" %in% input$s.plot || input$desc=="basis") {
      Bs <- fda::bsplineS(seq(0,1,length.out=1000), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1)*sd(Ts)
    }
	  #clust <- hclust(dist(t(Ts)),method="ward.D"); clcol <- cutree(clust,k=input$dimn);
	  if (input$desc=="scree") {
      svd.Ts <- svd(Ts); var <- (svd.Ts$d^2/sum(svd.Ts$d^2))*100;
      percentage <- paste0("\n The variation explained in first ",input$dimn," adaptive basis of the initial est. is ", sum(var[1:input$dimn]),"%")
      plot(var, xlab="Component", ylab="Percentage", main=paste("Scree Plot -",fname,percentage), pch=20, frame.plot=F, ylim=c(0,100))
      points(var[1:input$dimn],col=2);
    } else if (input$desc=="gcv") {
      res <- fda.gcv(); ind.m <- which(res$GCV==min(res$GCV)); ind.lam$fda <<- ind.m+lambdas[1]-1; 
      plot(res$lambdas, res$GCV, type="b", log="xy", xlab=bquote(lambda), ylab="GCV (black)", main=paste("Gen. Cross Validation / Deg. of Freedom -", fname), cex.lab=1.5, xaxt="n", pch=20); 
      abline(v=res$lambdas[ind.m], col=1, lty=2); axis(1,res$lambdas,paste0("2^",lambdas[1]:lambdas[2]))
      segments(x0 = res$lambdas[ind.m], y0=res$GCV[ind.m], x1=2^lambdas[1], col=1, lty=2);
      if (temp) {updateSliderInput(session,'tuning.fda',value=log2(res$lambdas[ind.m])); temp <<- 0}
      ind.c <- which(res$lambdas==2^input$tuning.fda); abline(v=res$lambdas[ind.c], col=2); 
      segments(x0 = res$lambdas[ind.c], y0=res$GCV[ind.c], x1=2^lambdas[1], col=2);
      par(new=T); plot(res$lambdas, res$DF, type="b", log="xy", xlab="", ylab="", xaxt="n", yaxt="n", col=4, pch=20); 
      segments(x0 = res$lambdas[ind.m], y0=res$DF[ind.m], x1=2^lambdas[2], col=4, lty=2);
      segments(x0 = res$lambdas[ind.c], y0=res$DF[ind.c], x1=2^lambdas[2], col=2);
      axis(4, round(res$DF),las=2); mtext(side = 4, line=-2, 'DF (blue)', cex=1.5)
    } else if (input$twoD) {
      if ("2f" == input$s.plot) {
        if (input$plot=="contour") contour(matrix(Ts[,indx],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="hmap") filled.contour(matrix(Ts[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="persp") persp(matrix(Ts[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
      } else if ("tf" == input$s.plot) {
        if (input$plot=="contour") contour(matrix(Trs[,indx],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="hmap") filled.contour(matrix(Trs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="persp") persp(matrix(Trs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
      } else if ("bf" == input$s.plot) {
        if (input$plot=="contour") contour(matrix(apply(as.matrix(B[,b.indx]),1,max),nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="hmap") filled.contour(matrix(apply(as.matrix(B[,b.indx]),1,max),nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="persp") persp(matrix(apply(as.matrix(B[,b.indx]),1,max),nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
      } else if ("bss" == input$s.plot) {
        if (input$plot=="contour") contour(matrix(as.matrix(B[,b.indx])%*%(cB%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="hmap") filled.contour(matrix(as.matrix(B[,b.indx])%*%(cB%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="persp") persp(matrix(as.matrix(B[,b.indx])%*%(cB%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
      } else if ("pbss" == input$s.plot) {
        if (input$plot=="contour") contour(matrix(as.matrix(B[,b.indx])%*%(cP%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="hmap") filled.contour(matrix(as.matrix(B[,b.indx])%*%(cP%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
        else if (input$plot=="persp") persp(matrix(as.matrix(B[,b.indx])%*%(cP%*%Ts[,indx])[b.indx,],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
      }; title(main=ifelse("bf"==input$s.plot, "Basis Functions", paste("2D Function -",fname,"- #",indx[1]))) 
    } else if (input$ams.choice=="multiple") {
      if (indx[2]-indx[1]>9) {updateSliderInput(session, 'sts.choice', value=c(indx[1],indx[1]+9)); indx[2] <- indx[1]+9}
      if ("ts" == input$s.plot) {plot.ts(Ts[,seq(indx[1],indx[2])],ann=F,ylim=range(Ts))}
      else if ("tf" == input$s.plot) {plot.ts(Trs[,seq(indx[1],indx[2])],ann=F,ylim=range(Ts))}
      else if ("bss" == input$s.plot) {plot.ts(as.matrix(B[,b.indx])%*%(cB%*%Ts[,seq(indx[1],indx[2])])[b.indx,],ann=F,ylim=range(Ts))}
      else if ("pbss" == input$s.plot) {plot.ts(as.matrix(B[,b.indx])%*%(cP%*%Ts[,seq(indx[1],indx[2])])[b.indx,],ann=F,ylim=range(Ts))}
      if ("bf" == input$s.plot) {
        if (max(b.indx)-b.indx[1]>9) {updateSliderInput(session, 'b.indx', value=c(b.indx[1],b.indx[1]+9)); b.indx <- b.indx[1]:(b.indx[1]+9)}
        plot.ts(as.matrix(B[,b.indx]),ann=F,ylim=range(Ts)); title(main="Basis Functions", xlab="Time")
      } else title(main=paste("Time Series -",fname,"- from",colnames(Ts)[indx[1]],"to",colnames(Ts)[indx[2]]), xlab="Time")
    } else {
      if (input$ams.choice=="all") {indx <- 1:ncol(Ts); m.lab <- fname} else {m.lab <- paste(fname,"-",colnames(Ts)[indx])}; 
      if ("ts" %in% input$s.plot) {clcol <- rep(1,ncol(Ts))} else {clcol <- rep(0,ncol(Ts))}
      ts.plot(Ts[,indx], col=clcol[indx], main=paste("Time Series -",m.lab), ylab="", ylim=range(Ts), gpars=list(xaxt="n")); 
      if ("bss" %in% input$s.plot) {
        for (i in indx) {
          f.est <- as.matrix(B[,b.indx])%*%(cB%*%Ts[,i])[b.indx];
          if (input$ams.choice=="single") if (input$s.CI) 
          {polygon(c(1:nrow(Ts),nrow(Ts):1), c(f.est-2*sqrt(vB),rev(f.est)+2*rev(sqrt(vB))), border=4, lwd=1,col=5); 
           if ("ts" %in% input$s.plot) points(Ts[,i],type="l",col=1)}
          points(f.est, col=4, type="l")
        }
      }
      if ("pbss" %in% input$s.plot) {
        for (i in indx) {
          p.est <- as.matrix(B[,b.indx])%*%(cP%*%Ts[,i])[b.indx];
          if (input$ams.choice=="single") if (input$s.CI) 
            {polygon(c(1:nrow(Ts),nrow(Ts):1), c(p.est-2*sqrt(vP),rev(p.est)+2*rev(sqrt(vP))), border=6, lwd=1,col=7);
             if ("ts" %in% input$s.plot) points(Ts[,i],type="l",col=1); if ("bss" %in% input$s.plot) {points(f.est, col=4, type="l")}}
          points(p.est, col=6, type="l");
        }
      }
      if ("tf" %in% input$s.plot) {for (i in indx) points(Trs[,i],type="l",col=2)}
      if ("bf" %in% input$s.plot) {for (i in b.indx) points(seq(1,nrow(Ts),length.out = 1000), Bs[,i],col=8,type="l",lty=2)}
      axis(1,trunc(summary(1:nrow(Ts))[-4]))
    }
  })
  
  output$runit <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="server") {fname <- names(Xs)[as.numeric(input$s.choice)]} else if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- "Simulation"}
    actionButton('runit', paste0('RUN FPCA (',fname,')'))
  })
  
  output$results <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (!(input$runit)) {textInput("results", "Click on RUN", "to obtain the FPCA")} 
    else {selectInput('results', 'Select Plot: ', choices = list(FPCA=c("FPCA Results"="FPCs", "Tuning Parameter (GCV)"="GCV"), Clustering = c("FPC Scores"="Scores", "Dendrogram"="Dendrogram")), selected="GCV", width="200px")}
  })
  
  runit <- eventReactive(input$runit, {
    if (input$twoD) lambdas <<- lambdas2 else lambdas <<- lambdas1;
    withProgress(message = 'FPCA: Tuning Parameter', value = 0, { GCV <- NULL;
      #if (input$diff!=0) {L <- diff(diag(nrow(Ts)), differences=input$diff)} else {L <- diag(nrow(Ts))}; P <- t(L)%*%L;
      try(ans <- get.FPCk(get.pen(1:nrow(Ts))$EIG.O, t(Ts), input$dimn, 2^(lambdas[1]:lambdas[2]), trunc.SVD=TRUE)); 
      if (nrow(ans$U)==1) updateSelectInput(session, 'results', choices = list(FPCA=c("FPCA Results"="FPCs", "Tuning Parameter (GCV)"="GCV")), selected = "GCV")
    }); 
    temp <<- 1; return(ans)
  })  
  
  output$ams.FPCA <- renderUI({ 
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$results!="FPCs" || input$twoD) return();
    radioButtons('ams.FPCA', 'Plot Choices:', c("All" = "all", "Multiple" = "multiple", "Single" = "single"), selected= "all", inline = TRUE, width="200px")
  })
  
  output$sts.FPCA <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$results!="FPCs" || !input$runit) return();
    if (input$twoD) {
      if (input$sn.plot=="fpcb") {return()}
      else {return(sliderInput("sts.FPCA", "Choose function:", min = 1, max = ncol(Ts), value = 1, step = 1, width="200px"))}
    } else {
      if (input$ams.FPCA=="multiple") 
        if (input$sn.plot=="fpcb") {return()}
        else {return(sliderInput("sts.FPCA", "Choose a range( '< 10' ):", min = 1, max = ncol(Ts), value = c(1,min(10,ncol(Ts))), dragRange=TRUE, step = 1, width="200px"))}
      if (input$ams.FPCA=="single") 
        if (((input$sn.plot=="fpcb" && length(input$sn.plot)==1) || !length(input$sn.plot)) && !length(intersect(input$so.plot,c("bss","pbss")))) {return()}
        else {return(sliderInput("sts.FPCA", "Choose function:", min = 1, max = ncol(Ts), value = 1, step = 1, width="200px"))}
    }
  })

  output$sn.plot <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$results!="FPCs") return();
    choices <- c("Time Series" = "ts", "True Functions" = "tf", "FPCA Basis" = "fpcb", "FPCA results" = "fpcr");
    if (input$f.choice!="sim") {choices <- choices[-2]}
    if (input$twoD) {choices <- c("2D Functions" = "2f", choices[-1]); return(radioButtons('sn.plot', 'Plot:', choices=choices, selected=choices[1], width="250px"))}
    else if (input$ams.FPCA=="multiple") {return(radioButtons('sn.plot', 'Plot:', choices=choices, selected=choices[1], width="250px"))}
    else {checkboxGroupInput('sn.plot', 'Plot:', choices=choices, selected=choices[1], width="250px")}
  })

  output$bn.indx <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || input$results!="FPCs" || !length(intersect(input$sn.plot,c("fpcb","fpcr")))) return();
    sliderInput("bn.indx", "FPC Basis Cont.", min = 1, max = ind.lam$D, value = c(1,ind.lam$D), dragRange=TRUE, step = 1, width="400px")
  })
  
  output$tuning.fpca <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || !input$runit) return();
    if (input$results=="FPCs" && !length(intersect(c("fpcr","fpcb"),input$sn.plot))) return()
    if (input$twoD) lambdas <<- lambdas2 else lambdas <<- lambdas1;
    sliderInput("tuning.fpca", "Tun. Par. (2^i):", min = lambdas[1], max = lambdas[2], value = ifelse(is.null(ind.lam$fpca),0,ind.lam$fpca), width="200px")
  })

  output$so.plot <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$twoD) return();
    if (input$ams.FPCA=="multiple" || input$results!="FPCs") return()
    choices <- c("Basis Func." = "bf", "B-Spline Smoothing" = "bss", "Pen. B-Sp. Smooth." = "pbss")
    if (input$results=="2f") choices <- c("2D Functions" = "2f", choices[-1]); 
    checkboxGroupInput('so.plot', 'FDA Plot:', choices=choices, width="250px")
  })
  
  output$bo.indx <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$twoD) return();
    if (input$ams.FPCA=="multiple" || input$results!="FPCs" || !length(intersect(input$so.plot,c("bf","bss","pbss")))) return()
    sliderInput("bo.indx", "Basis Contr.", min = 1, max = input$xdf, value = c(1,input$xdf), dragRange=TRUE, step = 1, width="200px")
  })
  
  output$tuningo.fda <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$twoD) return();
    if (input$ams.FPCA=="multiple" || input$results!="FPCs" || !"pbss"%in%input$so.plot ) return(); if (input$twoD) lambdas <<- lambdas2 else lambdas <<- lambdas1;
    sliderInput("tuningo.fda", "Tun. Par. (2^i):", min = lambdas[1], max = lambdas[2], value = ifelse(is.null(ind.lam$fda),0,ind.lam$fda), width="200px")
  })
  
  output$diffo <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model)) || input$twoD) return();
    if (input$ams.FPCA=="multiple" || input$results!="FPCs" || !"pbss"%in%input$so.plot || !input$diffP) return(); 
    sliderInput('diffo','Pen.Dif.Ord.', value = 2, min = 0, max=min(4,ifelse(input$twoD,min(input$xdf,input$ydf),input$xdf)-1), step=1, width="210px")
  })
  
  output$FPCA.out = renderPlot({if (memory.size()>700) gc();
    if ((input$f.choice=="upload" && is.null(input$file)) || (input$f.choice=="sim" && !length(input$model))) return();
    if (input$f.choice=="server") {fname <- names(Xs)[as.numeric(input$s.choice)]} else if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- "Simulation"}
    res <- runit(); ind.m <- which(res$GCV==min(res$GCV)); ind.c <- which(res$lambdas==2^input$tuning.fpca); 
    ind.lam$fpca <<- ind.m+lambdas[1]-1; ind.lam$D <<- ncol(res$D); rownames(res$U) <- colnames(Ts)
    res$u <- as.matrix(res$U[,,ind.c]); res$d <- as.matrix(res$D[,,ind.c]); res$v <- as.matrix(res$V[,,ind.c]);
    if (nrow(res$u)>1) {clust <<- hclust(dist(scale(res$u)), method="ward.D"); clcol <<- cutree(clust,k=input$dimn)} else {clcol=1}
    if (input$results=="Scores" && nrow(res$u)>1) {
      #for(i in 2:length(unique(tmp))) for(ii in 1:(i-1)) dis[i,ii] <- sum(as.vector(temp)*dens[,ii]*log(dens[,ii]/dens[,i])) + sum(as.vector(temp)*dens[,i]*log(dens[,i]/dens[,ii]))
      #clust <<- hclust(as.dist(dis),method="ward.D");   
      if (ncol(res$u)==1) {plot(res$u, col=clcol, pch=20, xaxt="n", ylab="FPC Score", xlab=""); axis(1,1:nrow(res$u),clust$labels,las=2)} 
        else {pairs(res$u, col=clcol, main=fname, cex.lab=2, label=paste("FPC Score",1:ncol(res$u)))}
    } else if (input$results=="Dendrogram" && nrow(res$u)>1) { plot(color_labels(as.dendrogram(clust), k = input$dimn), main=fname)
    } else if (input$results=="GCV") {
      plot(res$lambdas, res$GCV, type="b", log="xy", xlab=bquote(lambda), ylab="GCV (black)", main=paste("Gen. Cross Validation / Deg. of Freedom -", fname), cex.lab=1.5, xaxt="n", pch=20);
      axis(1,res$lambdas,paste0("2^",lambdas[1]:lambdas[2]))
      abline(v=res$lambdas[ind.m], col=1, lty=2); 
      segments(x0 = res$lambdas[ind.m], y0=res$GCV[ind.m], x1=2^lambdas[1], col=1, lty=2);
      if (temp) {updateSliderInput(session,'tuning.fpca',value=log2(res$lambdas[ind.m])); ind.c <-ind.m; temp <<- 0}
      abline(v=res$lambdas[ind.c], col=2); 
      segments(x0 = res$lambdas[ind.c], y0=res$GCV[ind.c], x1=2^lambdas[1], col=2);
      par(new=T); plot(res$lambdas, res$DF, type="b", log="xy", xlab="", ylab="", xaxt="n", yaxt="n", col=4, pch=20); 
      segments(x0 = res$lambdas[ind.m], y0=res$DF[ind.m], x1=2^lambdas[2], col=4, lty=2);
      segments(x0 = res$lambdas[ind.c], y0=res$DF[ind.c], x1=2^lambdas[2], col=2);
      axis(4, round(res$DF),las=2); mtext(side = 4, line=-2, 'DF (blue)', cex=1.5)
    } else {
      indx <- as.numeric(input$sts.FPCA); Fs <- t(as.matrix(res$u%*%res$d%*%t(res$v)));
      if (length(intersect(input$sn.plot,c("fpcb","fpcr")))) {
        bn.indx <- input$bn.indx[1]:input$bn.indx[2];
        tFs <- t(as.matrix(res$u[,bn.indx]%*%as.matrix(res$d[bn.indx,bn.indx])%*%t(res$v[,bn.indx])))
        tBs <- t(as.matrix(as.matrix(sqrt(res$d[bn.indx,bn.indx]))%*%t(res$v[,bn.indx])))
      }
      if (length(intersect(input$so.plot,c("bf","bss","pbss")))) bo.indx <- input$bo.indx[1]:input$bo.indx[2];
      if (length(intersect(input$so.plot,c("bss","pbss")))) {
        B <- fda::bsplineS(seq(0,1,length.out = ifelse(input$twoD,sqrt(nrow(Ts)),nrow(Ts))), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1);
        if ("pbss" %in% input$so.plot)
          if (input$diffP) {if (input$diffo!=0) {L <- diff(diag(input$xdf), differences=input$diffo)} else {L <- diag(input$xdf)}; P <- t(L)%*%L
          } else {P <- fda::bsplinepen(fda::create.bspline.basis(nbasis=input$xdf,norder=input$xdeg+1))};
        if (input$twoD) {
          B <- B%x%fda::bsplineS(seq(0,1,length.out = sqrt(nrow(Ts))), breaks=seq(0,1,length.out=input$ydf-input$ydeg+1), norder=input$ydeg+1);
          if ("pbss" %in% input$so.plot)
            if (input$diffP) {if (input$diffo!=0) {L <- diff(diag(input$ydf), differences=input$diffo)} else {L <- diag(input$ydf)}; P <- P%x%diag(ncol(L)) + diag(nrow(P))%x%(t(L)%*%L)
            } else {P2 <- fda::bsplinepen(fda::create.bspline.basis(nbasis=input$ydf,norder=input$ydeg+1)); P <- P%x%diag(nrow(P2)) + diag(nrow(P))%x%P2};
        }
        if ("bss" %in% input$so.plot) {cB <- solve(t(B)%*%B)%*%t(B)}; if ("pbss" %in% input$so.plot) {cP <- solve(t(B)%*%B+2^(input$tuningo.fda)*P)%*%t(B)}
      }
      if ("bf" %in% input$so.plot) {
        Bs <- fda::bsplineS(seq(0,1,length.out=1000), breaks=seq(0,1,length.out=input$xdf-input$xdeg+1), norder=input$xdeg+1)*sd(Ts)
      }
      if (input$twoD) {
        if ("2f" == input$sn.plot) {
          if (input$plot=="contour") contour(matrix(Ts[,indx],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="hmap") filled.contour(matrix(Ts[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="persp") persp(matrix(Ts[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
        } else if ("tf" == input$sn.plot) {
          if (input$plot=="contour") contour(matrix(Trs[,indx],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="hmap") filled.contour(matrix(Trs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="persp") persp(matrix(Trs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
        } else if ("fpcb" == input$sn.plot) {
          if (input$plot=="contour") contour(matrix(apply(tBs,1,max),nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="hmap") filled.contour(matrix(apply(tBs,1,max),nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="persp") persp(matrix(apply(tBs,1,max),nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
        } else if ("fpcr" == input$sn.plot) {
          if (input$plot=="contour") contour(matrix(tFs[,indx],nr=sqrt(nrow(Ts))),lwd=2, xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="hmap") filled.contour(matrix(tFs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", cex.lab=2)
          else if (input$plot=="persp") persp(matrix(tFs[,indx],nr=sqrt(nrow(Ts))), xlab="x", ylab="y", theta=input$p.ang1, phi=input$p.ang2, cex.lab=2)
        }; title(main=ifelse("fpcb"==input$sn.plot, "Basis Functions", paste("2D Function -",fname,"- #",indx[1])))
      } else if (input$ams.FPCA=="multiple") {
        if ("ts" == input$sn.plot) {plot.ts(Ts[,seq(indx[1],indx[2])],ann=F,ylim=range(Ts))}
        else if ("tf" == input$sn.plot) {plot.ts(Trs[,seq(indx[1],indx[2])],ann=F,ylim=range(Trs))}
        else if ("fpcr" == input$sn.plot) {plot.ts(tFs[,seq(indx[1],indx[2])],ann=F,axes=FALSE,frame=TRUE,ylim=range(Fs))}
        if ("fpcb" == input$sn.plot) {plot.ts(tBs,ann=F,axes=FALSE,frame=TRUE,ylim=range(tBs)); title(main="Basis Functions")}
        else {title(main=paste(input$results,"-",fname,"- from",colnames(Ts)[indx[1]],"to",colnames(Ts)[indx[2]]), xlab="Time")}
      } else {
        if (input$ams.FPCA=="all") {indx <- 1:ncol(Ts); m.lab <- fname} else {m.lab <- paste(fname,"-",colnames(Ts)[indx])};
        ts.plot(Fs[,indx], col=clcol[indx], main=paste(input$results,"-",m.lab), xlab="Time", ylab="", ylim=range(Fs), gpars=list(xaxt="n"));
        if ("ts" %in% input$sn.plot) {clcol <- rep(1,ncol(Ts))} else {clcol <- rep(0,ncol(Ts))}
        ts.plot(Ts[,indx], col=clcol[indx], main=paste("Time Series -",m.lab), ylab="", ylim=range(Ts), gpars=list(xaxt="n")); 
        if ("tf" %in% input$sn.plot) {for (i in indx) points(Trs[,i],type="l",col=2)}
        if ("fpcr" %in% input$sn.plot) {for (i in indx) points(tFs[,i],type="l",col=3, lwd=2)}
        if ("fpcb" %in% input$sn.plot) {for (i in bn.indx) points(sqrt(res$d[i,i])*res$v[,i],col=8,type="l",lty=3, lwd=2)}
        if ("bf" %in% input$so.plot) {for (i in bo.indx) points(seq(1,nrow(Ts),length.out = 1000), Bs[,i],col=8,type="l",lty=2)}
        if ("bss" %in% input$so.plot) {for (i in indx) points(as.matrix(B[,bo.indx])%*%(cB%*%Ts[,i])[bo.indx], col=4, type="l")}
        if ("pbss" %in% input$so.plot) {for (i in indx) points(as.matrix(B[,bo.indx])%*%(cP%*%Ts[,i])[bo.indx], col=6, type="l")}
        axis(1,trunc(summary(1:nrow(Ts))[-4]))
      }
    }
  })
}

shinyApp(ui=ui,server=server)