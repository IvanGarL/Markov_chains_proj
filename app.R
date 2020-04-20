#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(markovchain)
library(d3heatmap)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Análisis de "),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h2("Fase 1: "),
            radioButtons("ini_calidad",
                         "Estado inicial de calidad:",
                         c("Buena" = 1,
                           "Regular" = 2,
                           "Mala" = 3,
                           "Crítica" = 4)),
            sliderInput("semanas",
                        "Periodo de tiempo de observacion (semanas):",
                        min = 4,
                        max = 40,
                        value = 8,
                        step = 4),
            h2("Fase 2:"),
            numericInput("lim_pac",
                         "Limite de compromiso con pacientes:",
                         value = 20,
                         min = 10,
                         max = 30),
            numericInput("lim_ped",
                         "Mínimo posible de unidades de inventario físico:",
                         value = 15,
                         min = 10,
                         max = 30),
            numericInput("max_inv",
                         "Máximo posible de unidades de inventario físico:",
                         value = 40,
                         min = 35,
                         max = 60),
            sliderInput("pbim",
                        "Probabilidad de daño de medicamentos:",
                        min = 0.0,
                        max = 0.78,
                        value = 0.1,
                        step = 0.02),
            
            h2("Fase 3:"),
            sliderInput("edad",
                        "Edad inicial de paciente:",
                        min = 30,
                        max = 79,
                        value = 1),
            radioButtons("est_ini",
                         "Estado inicial de paciente:",
                         c("Sana" = "Sana",
                           "VPH" = "VPH",
                           "Estadío I" = "Estadio1",
                           "Estadío II" = "Estadio2", 
                           "Estadio III" = "Estadio3", 
                           "Estadio IV" = "Estadio4"))
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                id = 'dataset',
                tabPanel("Fase 1", plotOutput("sesionesplot")),
                tabPanel("Fase 2", plotOutput("distPlot"), plotOutput("costsPlot")),
                tabPanel("Fase 3", 
                         textOutput("e_optimo"),
                         h4("Mapa de calor de Politica optima"),
                         d3heatmapOutput("heatmap", width = "100%", height="600px"),
                         h4("Tabla de decisiones optimas"),
                         DT::dataTableOutput("tabla_opt")
                         )
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    #Funciones ------------EMPIEZA FASE 1:
    
    #calcula el valor esperado basado en las probabilidades en estado estable y un vector x.
    valorEsperado <- function (staStates, vector){
        
        ve <- 0
        
        for(i in 1:length(staStates)){
            ve = ve + (staStates[i]*vector[i]);
        }
        
        return(ve)
    }
    
    #Retorna el alpha respectivo al estado inicial
    getAlpha <- function(estado){
        if(estado == 1){
            return(c(1,0,0,0))
        }
        else if(estado == 2){
            return(c(0,1,0,0))
        }
        else if(estado == 3){
            return(c(0,0,1,0))
        }
        else if(estado == 4){
            return(c(0,0,0,1))
        }
    }
    
    #Calcula el valor esperado de sesiones para una region en un periodo de semanas
    #   Recuerdese que las sesiones se hacen cada dos semanas, por eso el ciclo  salta de 2 en 2
    valorEsperadoSesiones <- function(alpha, cmtd_region, periodo, sesiones){
        
        ve_total = 0
        for(i in seq(2, periodo, by = 2)){
            
            #Se empiezan las sesiones a partir de la 2da semana.
            P_i <- alpha * (cmtd_region^i)
            ve_region_i <- valorEsperado(sesiones, P_i)
            
            ve_total <- ve_total + ve_region_i
        }
        
        return(ve_total)
    }
    
    
    #---server output - sesionesplot:
    output$sesionesplot <- renderPlot({
        
        #INPUTS:
        #alpha & periodo de observacion
        alpha_in = getAlpha(input$ini_calidad)
        periodo_in = input$semanas
        
        #Carga de matrices de probabilidad
        atlantica <- read.csv(file = "atlantica.csv", sep = ",")
        bogota <- read.csv(file = "bogota.csv", sep = ",")
        central <- read.csv(file = "central.csv", sep = ",")
        oriental <- read.csv(file = "oriental.csv", sep = ",")
        pacifica <- read.csv(file = "pacifica.csv", sep = ",")
        
        #Conversion de Data frames a matrices para inicializar las markovchain
        matAtl = as.matrix(atlantica[1:4,1:4], rownames.force = TRUE)
        matBog = as.matrix(bogota[1:4,1:4], rownames.force = TRUE)
        matCen = as.matrix(central[1:4,1:4], rownames.force = TRUE)
        matOri = as.matrix(oriental[1:4,1:4], rownames.force = TRUE)
        matPac = as.matrix(pacifica[1:4,1:4], rownames.force = TRUE)
        
        #Notacion de los estados :
        # 1 : buena
        # 2 : regular
        # 3 : mala
        # 4 : critica
        
        #Guardo los estados en un vector para la clase "markovchain"
        estados = c(1:4)
        estadosChar = as.character(estados)
        
        #Inicializacion de las markovchain
        cmtdAtl <- new(Class="markovchain", states = estadosChar, transitionMatrix = matAtl)
        cmtdBog <- new(Class="markovchain", states = estadosChar, transitionMatrix = matBog) 
        cmtdCen  <- new(Class="markovchain", states = estadosChar, transitionMatrix = matCen) 
        cmtdOri <- new(Class="markovchain", states = estadosChar, transitionMatrix = matOri) 
        cmtdPac <- new(Class="markovchain", states = estadosChar, transitionMatrix = matPac) 
        
        #Vector de sesiones segun el enunciado:
        ses <- c(0,1,3,5)
        
        #Calculo el valor esperado total de sesiones para cada sesion
        ve_sesiones_Atl <- valorEsperadoSesiones(alpha = alpha_in, cmtd_region = cmtdAtl, periodo = periodo_in, sesiones = ses)
        ve_sesiones_Bog <- valorEsperadoSesiones(alpha = alpha_in, cmtd_region = cmtdBog, periodo = periodo_in, sesiones = ses)
        ve_sesiones_Cen <- valorEsperadoSesiones(alpha = alpha_in, cmtd_region = cmtdCen, periodo = periodo_in, sesiones = ses)
        ve_sesiones_Ori <- valorEsperadoSesiones(alpha = alpha_in, cmtd_region = cmtdOri, periodo = periodo_in, sesiones = ses)
        ve_sesiones_Pac <- valorEsperadoSesiones(alpha = alpha_in, cmtd_region = cmtdPac, periodo = periodo_in, sesiones = ses)
        
        #Se dibuja el plot con los resultados obtenidos:
        data_sesiones <- c(ve_sesiones_Atl, ve_sesiones_Bog, ve_sesiones_Cen, ve_sesiones_Ori, ve_sesiones_Pac)
        
        max = max(data_sesiones) + 5
        titulo = paste("Valor esperado de sesiones por region en un periodo de ", periodo_in, " semanas")
        
        bar = barplot(data_sesiones, 
                main= titulo,
                ylab="VE sesiones", 
                names.arg = c("Atlántica", "Bogotá", "Central", "Oriental", "Pacífica"),
                col = "gray",
                ylim = c(0,max)
        )
        
        text( bar, data_sesiones, labels = round(data_sesiones, 3), pos = 3, adj = 1)
        
    })
    #Funciones ----------------FIN FASE 1:
    
    #Funciones ------------EMPIEZA FASE 2:
    
    #Retorna la cadena de markov para la formulación correspondiente a los parámetros.
    getMarkovChain <- function(lda, u, m, k, pbim){
        estados = c(-u:k)
        Mat = matrix(0, nrow=length(estados), ncol=length(estados))
        
        for(i in estados){
            for(j in estados){
                #a
                if(j>=0 & i<m){
                    for(l in (0:k-j)){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+(dpois(k-j-l, lambda = lda)*dbinom(x = l,size = k, prob = pbim))
                    }
                }
                #b
                if(j>(-u) & j<0 & i<m){
                    for(l in 0:k){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+(dpois(k-j-l, lambda = lda)*dbinom(x = l,size = k, prob = pbim))
                    }
                }
                #c
                if(j>=0 & j<=i & i>=m){
                    for(l in 0:(i-j)){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+(dpois(i-j-l, lambda = lda)*dbinom(x = l,size = i, prob = pbim))
                    }
                }
                #d
                if(j>(-u) & j<0 & i>=m & j<=i){
                    for(l in 0:i){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+(dpois(i-j-l, lambda = lda)*dbinom(x = l,size = i, prob = pbim))
                    }
                }
                #e
                if(j==(-u) & i<m){
                    for(l in 0:k){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+((1-ppois(k-j-l-1, lambda = lda, lower.tail = TRUE))*dbinom(x = l,size = k, prob = pbim))
                    }
                }
                
                #f
                if(j==(-u) & i>=m){
                    for(l in 0:i){
                        Mat[i+(u+1),j+(u+1)]= Mat[i+(u+1),j+(u+1)]+((1-ppois(i-j-l-1, lambda = lda, lower.tail = TRUE))*dbinom(x = l,size = i, prob = pbim))
                    }
                }
            }
        }
        estados_char<-as.character(estados)
        cmtdActual<-new("markovchain", states=as.character(estados_char), transitionMatrix=Mat)
        
        return(cmtdActual)
    }
    
    #Estima el costo de la politica segun los parametros 
    estimateCosts <- function(stSts, cost_prov, cost_domc, cost_invt, u, k, m){
        
        estados = c(-u:k)
        lim <- length(estados)
        costoEstimado<-0
        
        for(e in 1:(lim-1))
        {
            if(e<(u+1)) #Si tengo faltantes y pido
            {
                costoEstimado<- (costoEstimado + ((((((lim+1)-e)*cost_prov)+ (((u+1)-e)*cost_domc)))*stSts[e]))
            }
            if(e>=(u+1) & e<(k-m+1)) #Si pido pero no tengo faltantes
            {
                costoEstimado<- (costoEstimado + ((((((lim+1)-e)*cost_prov)+ ((e-(u+1))*cost_invt)))*stSts[e]))
            }
            
            if(e>=(k-m+1))  #Si no pido y no tengo faltantes
            {
                costoEstimado<- (costoEstimado + ((((e-(u+1))*cost_invt))*stSts[e]))
            }
        }
        
        return(costoEstimado)
    }
    
    #---server output - distplot:
    output$distPlot <- renderPlot({
        
        #INPUTS:
        #Limites & probabilidad
        lda_f = 24.99
        u_in = input$lim_pac # Limite de compromiso con pacientes
        m_in = input$lim_ped #Min de inv. fisico posible
        k_in =  input$max_inv #Max de inv. fisico posible
        pbim_in = input$pbim #Probabilidad de tener que botar las medicinas
        
        cmtd_all <-getMarkovChain(lda = lda_f, u = u_in, m = m_in, k = k_in, pbim = pbim_in )
        
        estables_plot <- steadyStates(cmtd_all)
        
        barplot(estables_plot, 
                main="Distribucion de probabilidades en estado estable",
                ylab="Probabilidad", 
                names.arg = c(-u_in:k_in),
                col = "gray"
        )
        
    })
    
    #---server output - costsplot:
    output$costsPlot <- renderPlot({
        #Valores fijos de parametros: limites y costos
        #Limites, media y pbinomial
        lda_f = 24.99 #media de la demanda -> Poisson
        u_f = 20 #limite de compromiso con clientes
        m_f = 15 #Minimo de inv. para reorden
        k_f = 40 #Max de capacidad de inv. posible
        pbim_f = 0.1 #Probabilidad de da;o en la binomial
        
        #Costos
        prov_f = 25000 #Proveedores
        domc_f = 5000 #Domicilio
        invt_f = 7500 # Inventario
        
        #INPUTS:
        #Limites & media
        u_in = input$lim_pac # Limite de compromiso con pacientes
        m_in = input$lim_ped #Min de inv. fisico posible
        k_in =  input$max_inv #Max de inv. fisico posible
        pbim_in = input$pbim #Probabilidad de tener que botar las medicinas
        
        
        #Calculamos cada CMTD variando unicamente una categoria para cada una respectivamente
        cmtd_inic <-getMarkovChain(lda = lda_f, u = u_f, m = m_f, k = k_f, pbim = pbim_f)
        cmtd_limPac <-getMarkovChain(lda = lda_f, u = u_in, m = m_f, k = k_f, pbim = pbim_f)
        cmtd_limPed <-getMarkovChain(lda = lda_f, u = u_f, m = m_in, k = k_f, pbim = pbim_f)
        cmtd_limMax <-getMarkovChain(lda = lda_f, u = u_f, m = m_f, k = k_in, pbim = pbim_f)
        cmtd_pBim <-getMarkovChain(lda = lda_f, u = u_f, m = m_f, k = k_f, pbim = pbim_in)
        
        #Calculamos sus respectivos estados estables para poder calcular los costos
        estables_inic <- steadyStates(cmtd_inic)
        estables_limPac <- steadyStates(cmtd_limPac)
        estables_limPed <- steadyStates(cmtd_limPed)
        estables_limMax <- steadyStates(cmtd_limMax)
        estables_pBim <- steadyStates(cmtd_pBim)
        
        #Calculamos los costos de cada cmtd
        costos_inic <- estimateCosts(stSts = estables_inic, cost_prov = prov_f, cost_domc = domc_f, cost_invt = invt_f, u = u_f, k = k_f, m = m_f)
        costos_limPac <- estimateCosts(stSts = estables_limPac, cost_prov = prov_f, cost_domc = domc_f, cost_invt = invt_f, u = u_in, k = k_f, m = m_f)
        costos_limPed <- estimateCosts(stSts = estables_limPed, cost_prov = prov_f, cost_domc = domc_f, cost_invt = invt_f, u = u_f, k = k_f, m = m_in)
        costos_limMax<- estimateCosts(stSts = estables_limMax, cost_prov = prov_f, cost_domc = domc_f, cost_invt = invt_f, u = u_f, k = k_in, m = m_f)
        costos_pBim<- estimateCosts(stSts = estables_pBim, cost_prov = prov_f, cost_domc = domc_f, cost_invt = invt_f, u = u_f, k = k_f, m = m_f)
        
        data <- c(costos_inic, costos_limPac, costos_limPed, costos_limMax, costos_pBim)
        max = max(data) + 100000
        
        
        bar = barplot(data, 
                main="Variacion de costos con respecto a la variación de sus parámetros",
                ylab="Valor del costo", 
                names.arg = c("Ningún parámetro variado", "Límite de inv. negativo", "Mínimo inv. físico", "Máximo inv. físico", "Probabilidad daño mercancia"),
                col = "darkblue",
                ylim = c(0,max)
              )
        
        text( bar, data, labels = round(data, 3), pos = 3, adj = 1)
        
    })
    #Funciones ----------------FIN FASE 2:
    
    #Funciones ------------EMPIEZA FASE 3:
    setUpSDP <- function(){
        #Espacio de estados
        
        estados<-c(1,2,3,4,5,6,7)
        names(estados)<-c("Sana","VPH", "Estadio1",  "Estadio2", "Estadio3", "Estadio4", "Fallecida")
        
        #INPUTS:
        edad_inicial = (input$edad - 29)
        estado_inic = estados[input$est_ini]
        
        #Conjunto de epocas
        epocas<-seq(1,50, by=1)
        
        #Conjunto de decisiones
        A<-c("Esperar", "Medicamento", "Radioterapia", "Quimioterapia")
        
        #Carga de matrices de probabilidad
        esperar <- read.csv(file = "m_esperar.csv", sep = ",")
        medicamento <- read.csv(file = "m_medicamento.csv", sep = ",")
        quimio <- read.csv(file = "m_quimio.csv", sep = ",")
        radio <- read.csv(file = "m_radio.csv", sep = ",")
        
        #Convertimos cada matriz a Dataframe
        matEsp = as.matrix(esperar, rownames.force = FALSE)
        matMed = as.matrix(medicamento, rownames.force = TRUE)
        matQuim = as.matrix(quimio, rownames.force = TRUE)
        matRad = as.matrix(radio, rownames.force = TRUE)
        
        #Rotulamos cada columna
        rownames(matEsp)<-names(estados)
        colnames(matEsp)<-names(estados)
        rownames(matMed)<-names(estados)
        colnames(matMed)<-names(estados)
        rownames(matQuim)<-names(estados)
        colnames(matQuim)<-names(estados)
        rownames(matRad)<-names(estados)
        colnames(matRad)<-names(estados)
        
        #Definicion de vectores de retorno inmediato para cada decision en cada estado.
        rEsperar<-c(1, 0.98, 0.65, 0.67, 0.5, 0.2, 0)
        names(rEsperar)<-estados
        
        rMedicamento<-(0.95)
        
        rRadioterapia<-c(0.8, 0.7)
        names(rRadioterapia)<-c("Estado1", "Estado2")
        
        rQuimioterapia<-c(0.75, 0.65, 0.5)
        names(rQuimioterapia)<-c("Estado2", "Estado3", "Estado4")
        
        
        #Crear tabla de estados y retornos innmediatos
        tablaRetornos<-matrix(0,ncol=5, nrow=7, byrow = TRUE)
        colnames(tablaRetornos)<-c("Estado", "RetEsperar", "RetMedicamento", "RetRadioterapia", "RetQuimioterapia")
        for (i in (1:7)) {
            tablaRetornos[i,1]=names(estados[i])
        }
        for (i in (1:7)) {
            tablaRetornos[i,2]=rEsperar[i]
        }
        tablaRetornos[2,3]=0.95
        for (i in (1:2)) {
            tablaRetornos[i+2,4]=rRadioterapia[i]
        }
        for (i in (1:3)) {
            tablaRetornos[i+3,5]=rQuimioterapia[i]
        }
        
        Tmax <- 50
        #creamos la  matriz de retornos para cada una de las epocas
        Ft <- matrix(0,nrow=7,ncol=Tmax)
        
        #Matriz para guardar las decisiones óptimas para cada una de las épocas
        Mat_Dec<-matrix(0,nrow=7,ncol=Tmax)
        
        #Matriz Q de retornos posibles para la ultima epoca 
        Q<- matrix(0,nrow=7,ncol=length(A))
        colnames(Q)<-A
        Q<- tablaRetornos[,2:5]
        
        #Se guardan la solucion optima de la ultima etapa para cada estado
        Ft[,Tmax] <- apply(Q, 1, max)
        
        #Decisiones optimas para cada estado para cada etapa/epoca
        dec_opt_por_etapa<-c()
        
        #Recorrido para determinar cuales son las decisiones optimas para cada estado
        for(j in 1:7){
            dec_opt_por_etapa<- c(dec_opt_por_etapa,A[which.max(Q[j,])])  
        }
        
        #Guardamos el vector de decisiones optimas en la última época
        Mat_Dec[,Tmax]<-dec_opt_por_etapa
        
        #Se define el factor de descuento
        factordesc <- 1
        
        #Iteraciones desde la época 49 hasta la época 1
        for (t in (Tmax-1):1) {
            
            #Inicializamos la matriz de retornos de la presente época
            Q <- array(0, dim=c(7,length(A)))
            
            #Recorro sobre decisiones y estados
            for (i in 1:length(A)) {
                for (j in 1:7) {
                    #Si la decision es esperar
                    if(A[i]=="Esperar"){
                        MatProb<-matEsp
                        Q[j,i] <- as.numeric(tablaRetornos[j,2]) + factordesc*
                            (MatProb[j,] %*% 
                                 as.numeric(Ft[,t+1]))     
                    }
                    #Si la decision es tomar el medicamento
                    else if(A[i]=="Medicamento"){
                        MatProb<-matMed
                        Q[j,i] <- as.numeric(tablaRetornos[j,3]) + factordesc*
                            (MatProb[j,] %*% 
                                 as.numeric(Ft[,t+1]))     
                    }
                    #Si la decision es ir a Radioterapia
                    else if(A[i]=="Radioterapia"){
                        MatProb<-matRad
                        Q[j,i] <- as.numeric(tablaRetornos[j,4]) + factordesc*
                            (MatProb[j,] %*% 
                                 as.numeric(Ft[,t+1]))      
                    }
                    #Si la decision es ir a Quimioterapia
                    else if(A[i]=="Quimioterapia"){
                        MatProb<-matQuim
                        Q[j,i] <- as.numeric(tablaRetornos[j,5]) + factordesc*
                            (MatProb[j,] %*% 
                                 as.numeric(Ft[,t+1]))
                    }
                    #Si no es posible tomar una decisión se penaliza el retorno inmediato(-999999999).
                    else{
                        Q[j,i]<- -999999999
                    }
                }  
            }
            
            #Encontramos el vector de valores asociados a las decisiones optimas
            Ft[,t] <- apply(Q, 1, max)
            #Identificar la decision optima para cada estado en la presente etapa
            dec_opt_por_etapa<-c()
            for(j in 1:7){
                dec_opt_por_etapa<- c(dec_opt_por_etapa,A[which.max(Q[j,])])  
            }
            #Guardamos el vector de decisiones optimas de la presente época
            Mat_Dec[,t]<-dec_opt_por_etapa
        }
        
        matrices <- list(Ft, Mat_Dec, estado_inic, A, edad_inicial)
        names(matrices) <- c("Ft", "Mat_Dec", "estado_inic", "A", "edad_inicial")
        
        return(matrices)
    } 
    
    drawHeatMap <- function (){
        mats <- setUpSDP()
        Mat_Dec <- mats[["Mat_Dec"]]
        A <- unlist(mats["A"][1], use.names = FALSE)
        
        
        # PRIMER OUTPUT
        # VALOR ESPERADO DE VIDA  
        #
        Mat_Dec_map = matrix(0, nrow = 7, ncol = 50)
        for(state in 1:7){
            for(time in 1:50){
                dec = which(A==Mat_Dec[state,time])
                Mat_Dec_map[state,time] = dec
            }
        }
        
        #Se crea la grafica para exponer la politica optima.
        return(d3heatmap(Mat_Dec_map, scale = "column",
                         colors = "Spectral",
                         dendrogram = "none",
                         Rowv = FALSE,
                         Colv = FALSE)
        )
    }
    
    #---server output - e_optimo:
    output$e_optimo <- renderText({
        mats <- setUpSDP()
        Ft <- mats[["Ft"]]
        estado_inic <- mats[["estado_inic"]]
        edad_inicial <- mats[["edad_inicial"]]
        
        resp = paste("El valor esperado de la esperanza de vida segun los parametros es: ", Ft[estado_inic, edad_inicial])
        return(resp) #Valor esperado de la espectativa de vida para los input de la app
    })
    
    #---server output - heatmap:
    output$heatmap <- renderD3heatmap({
        drawHeatMap()
    })
    
    
    #---server output - tabla_opt:
    output$tabla_opt <- DT::renderDataTable({
        mats <- setUpSDP()
        Mat_Dec <- mats[["Mat_Dec"]]
        estado_inic <- unlist(mats["estado_inic"][1], use.names = FALSE)
        
        #
        #PROCESO DE HALLAR LA POLITICA OPTIMA
        #
        #Creamos la matriz para guardar la política óptima
        Pol_Opt = matrix(0,nrow = 1,ncol=3)
        colnames(Pol_Opt)=c("Epoca", "Estado","Decision")
        
        #Muestras aleatorias.
        getRandomState <- function(est_act, decs){
            
            if(est_act == 1 & decs == "Esperar"){
                samp = sample(1:2, replace = F)
                return(samp[1])
            }
            if(est_act == 2 & (decs == "Esperar" || decs == "Medicamento")){
                samp = sample(1:3, replace = F)
                return(samp[1])
            }
            if(est_act == 3 & (decs == "Esperar" || decs == "Radioterapia")){
                samp = sample(2:4, replace = F)
                return(samp[1])
            }
            if(est_act == 4 & (decs == "Esperar" || decs == "Radioterapia" || decs == "Quimioterapia")){
                samp = sample(2:4, replace = F)
                return(samp[1])
            }
            if(est_act == 5 & (decs == "Esperar" || decs == "Quimioterapia")){
                samp = sample(5:6, replace = F)
                return(samp[1])
            }
            if(est_act == 6 & (decs == "Esperar" || decs == "Quimioterapia")){
                samp = sample(6:7, replace = F)
                return(samp[1])
            }
        }
        
        #creamos vectores auxiliares para guardar temporalmente 
        #la decisión y estado de llegada para cada época
        decisiones=c()
        estadollegada= -1;
        estado_epoca_actual = estado_inic
        Tmax <- 50
        
        for(epoca in 1:Tmax){
            decision <- Mat_Dec[estado_epoca_actual, epoca]
            decisiones[epoca] = decision
            
            estadollegada = getRandomState(estado_epoca_actual, decision)
            
            Pol_Opt = rbind(Pol_Opt, c(epoca,estado_epoca_actual, decision))
            estado_epoca_actual <- estadollegada
        }
        
        return(DT::datatable(Pol_Opt))
    })
    #Funciones ----------------FIN FASE 3:
}

# Run the application 
shinyApp(ui = ui, server = server)
