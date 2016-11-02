#####################################
#####################################
## Functions for CellAuto




##########
## Init landscape

##' .. Generate initial landscape for cellular automaton simulation ..
##'
##' .. Take a data.frame of species parameters and generate initial condition with random assignment of specie sto cells with a given size of the landscape (defined as square landscape of size NN repeated Nlandscape times to cover the abiotic stress gradient). The abiotic stress gradient is linear form min_ss to max_ss and run from along NN*Nlandscape. ..
##' @title InitLandscape
##' @param df_sp data.frame with species parameter
##' @param NN size of square landscape
##' @param Nlandscape number of time the landscape is repeated
##' @param min_ss minimun abiotic stress value
##' @param max_ss maximum abiotic stress value
##' @param per_occup percentage of cell occuped by a species
##' @export
##' @return list with first element a matrix (NN, NN*Nlandscape) with species number (-100 no species), second element matrix off successional stage of the cells, third element a vector of teh abiotic gradient stress (size NN*Nlandscape), fourst element species early successional competitive ability, fiveth element species late successional competitive ability, sixt element species abiotic stress tolerance.
##' @author Georges Kunstler
InitLandscape <- function(df_sp, NN = 100, Nlandscape = 1, min_ss = 0, max_ss = 0, per_occup = 0.5){
   matrix_sp <- matrix(0, nrow=NN,ncol=Nlandscape*NN)
   matrix_suc <- matrix(1, nrow=NN,ncol=Nlandscape*NN)
   climate_grad <-  seq(from= min_ss,to= max_ss,length=Nlandscape*NN) # not used yet
   ## Select cells to init landscape with species
   init.temp <- sample(1:(NN*NN*Nlandscape),size=round(NN*NN*Nlandscape*per_occup))
   ## fill with species
   nsp <- nrow(df_sp)
   matrix_sp[init.temp] <- sample(1:nsp, length(init.temp), replace = TRUE)

   ## species param
   sp_e <- c(-100, df_sp$E)
   sp_l <- c(-100, df_sp$L)
   sp_s <- c(-100, df_sp$S)

   return(list(mat_sp = matrix_sp, mat_suc = matrix_suc, ss = climate_grad,
               c_e = sp_e, c_l = sp_l, c_s = sp_s))
}

## generate species
##' .. Create the data.frame describing a set of nsp species with a triangular constrain on their three parameters (early successional competitive ability, late successional competitive ability, abiotic stress tolerance). ..
##'
##' .. Generate nsp equally spaced species in traiangular space of three parameters (early successional competitive ability, late successional competitive ability, abiotic stress tolerance) representing teh species strategies. These species are randomly drawn from Nval value of the parameters. ..
##' @title GenerateRandSp
##' @param nsp number of species
##' @param Nval number of value to drwan from (must be greater than nsp)
##' @return data.frame of the three species parameters (E, L, S)
##' @export
##' @author Georges Kunstler
GenerateRandSp <- function(nsp, Nval = 1000){
   if(nsp > 400) stop("error nsp must be smaller than 400")
   E <- matrix(rep(0:Nval, Nval + 1), Nval +1, Nval+1)/Nval
   L <- matrix(rep(0:Nval, Nval + 1), Nval + 1, Nval +1, byrow= TRUE)/Nval
   S <- (E+L)
   E.r <-  E[rev(0:(Nval + 1)), ]
   L.r <-  L[rev(0:(Nval + 1)), ]
   S.r <-  S[rev(0:(Nval + 1)), ]
   E.r[upper.tri(E.r,diag = FALSE)] <-  NA
   L.r[upper.tri(L.r,diag = FALSE)] <-  NA
   S.r[upper.tri(S.r,diag = FALSE)] <-  NA
   E.b <- E.r[rev(0:(Nval + 1)), ]
   L.b <- L.r[rev(0:(Nval + 1)), ]
   S.b <- S.r[rev(0:(Nval + 1)), ]

   df <- na.omit(data.frame(E = as.vector(E.b),
                            L = as.vector(L.b),
                            S = as.vector(S.b)))
  df.sp <- df[sample(1:dim(df)[1], nsp+50), ]
  df.sp <- df.sp[!duplicated(df.sp$E) & !duplicated(df.sp$L), ]
if(nrow(df.sp) > nsp) df.sp <- df.sp[sample(1:nrow(df.sp), nsp), ]
   return(df.sp)
}

##' .. Create a plot of species three parameters early successional competitive ability, late successional competitive ability, abiotic stress tolerance..
##'
##' .. content for \details{} ..
##' @title TernPlotSp
##' @param Nval number of species used to generate the plot.
##' @return plot
##' @export
##' @author Georges KUnstler
TernPlotSp <- function(Nval = 1000){
   require(raster)
   E <- matrix(rep(0:Nval, Nval + 1), Nval +1, Nval+1)/Nval
   L <- matrix(rep(0:Nval, Nval + 1), Nval + 1, Nval +1, byrow= TRUE)/Nval
   S <- (E+L)
   E.r <-  E[rev(0:(Nval + 1)), ]
   L.r <-  L[rev(0:(Nval + 1)), ]
   S.r <-  S[rev(0:(Nval + 1)), ]
   E.r[upper.tri(E.r,diag = FALSE)] <-  NA
   L.r[upper.tri(L.r,diag = FALSE)] <-  NA
   S.r[upper.tri(S.r,diag = FALSE)] <-  NA
   E.b <- E.r[rev(0:(Nval + 1)), ]
   L.b <- L.r[rev(0:(Nval + 1)), ]
   S.b <- S.r[rev(0:(Nval + 1)), ]

 trade.array <- stack(flip(raster(E.b), 2),
                      flip(raster(L.b), 2),
                      flip(raster(1 - S.b), 2))
plotRGB(trade.array,scale=1,axes=FALSE)
## text(0.15,0.1,labels="Stress tolerant", cex = 1.5)
## text(0.85,0.1,labels="Comp late succ", cex = 1.5)
## text(0.15,0.9,labels="Comp early succ", cex = 1.5)
}

##' .. Fill the matrix with species parameters ..
##'
##' .. content for \details{} ..
##' @title return_mat_fill_c
##' @param m matrix with species number
##' @param cc species parameters value
##' @return matrix
##' @author Georges Kunstler
return_mat_fill_c <- function(m, cc){
  d <-  dim(m)
  cc[1] <- NA
  mat <- matrix(cc[m+1], d[1], d[2])
  mat
}
##' .. Fill the matrix with species parameters adn return it as a raster object ..
##'
##' .. content for \details{} ..
##' @title return_mat_fill_c_raster
##' @param m matrix with species number
##' @param cc species parameters value
##' @return raster
##' @author Georges Kunstler
return_mat_fill_c_raster <- function(m, cc){
  require(raster)
  mat <- return_mat_fill_c(m, cc)
  return(raster(mat,
                xmn=1, xmx=dim(m)[2],
                ymn=1, ymx=dim(m)[1]))
}
##' .. Convert matrix to raster for plot ..
##'
##' .. content for \details{} ..
##' @title return_mat_raster
##' @param m matrix
##' @return raster
##' @author Georges Kunstler
return_mat_raster <- function(m){
  require(raster)
  d <-  dim(m)
  return(raster(m,
                xmn=1, xmx=d[2],
                ymn=1, ymx=d[1]))
}

### plot image of landscape Succ
##' .. plot matrix as a raster..
##'
##' .. content for \details{} ..
##' @title image_return_mat
##' @param m matrix
##' @param ...
##' @return plot
##' @export
##' @author Georges Kunstler
image_return_mat<-  function(m, ...){
    rast <-  return_mat_raster(m)
    image(rast,axes=FALSE ,asp=1, ...)
 }

##' .. Plot a landscape list with one species parameters. ..
##'
##' .. content for \details{} ..
##' @title image_landscape_e
##' @param list_res landscape list
##' @param c_ species parameters
##' @return plot
##' @export
##' @author Georges Kunstler
image_landscape_e<-  function(list_res, c_){
    rast_e <-  return_mat_fill_c_raster(list_res[[1]], c_)
    image(rast_e,axes=FALSE ,asp=1)
 }

##' .. Plot RGB of species three parameters over the landscape (one color for each parameter). ..
##'
##' .. content for \details{} ..
##' @title image_landscape
##' @param list_res landcsape list
##' @param c_e vector of species parameters c
##' @param c_l vector of species parameters l
##' @param c_s  vector of species parameters s
##' @param ...
##' @export
##' @return plot
##' @author Georges Kunstler
image_landscape <-  function(list_res, c_e, c_l, c_s, ...){
    rast_e <-  return_mat_fill_c_raster(list_res[[1]], c_e)
    rast_l <-  return_mat_fill_c_raster(list_res[[1]], c_l)
    rast_s <-  return_mat_fill_c_raster(list_res[[1]], 1 - c_s)
    rast.temp <- stack(rast_e,
                       rast_l,
                       rast_s)
    plotRGB(rast.temp,scale=1,axes=FALSE ,asp=1, ...)
 }


## test convergence
##' .. Run sequence of Nrun simulation Niter times to explore convergence. ..
##'
##' .. content for \details{} ..
##' @title eval_converg
##' @param list_init Initila list of the landscape
##' @param Niter number of iteration
##' @param Nrun number of cellular automaton update run per iteration
##' @param p_d mortality rate paramter
##' @param p_s succession rate paramter
##' @param K strength of the competitive asymmetry parameter
##' @export
##' @return list of landscape list
##' @author Georges Kunstler
eval_converg<- function(list_init, Niter, Nrun, p_d, p_s, K){

  list_res <- vector('list', Niter+1)
  list_res[[1]] <- list_init
  print(1)
  list_res[[2]] <- UpdateIterR(list_init$mat_sp, list_init$mat_suc,
                       list_init$c_e, list_init$c_l, list_init$c_s,
                       list_init$ss,
                       p_d, p_s, K , n = Nrun)
  for (i in 2:Niter){
     print(i)
     list_res[[i+1]] <- UpdateIterR(list_res[[i]]$sp, list_res[[i]]$suc,
                       list_init$c_e, list_init$c_l, list_init$c_s,
                       list_init$ss,
                       p_d, p_s, K , n = Nrun)
  }
  return(list_res)
}


## coexistence
##' .. Run simulation of cellur automaton for all combination mortality rate and succession rate to evalute their effect on species coexistence ..
##'
##' .. content for \details{} ..
##' @title eval_coex
##' @param p_d vector of mortality rate to explore
##' @param p_s vector of succession rate to explore
##' @param K strength of the competitive asymmetry parameter
##' @param list_int list of initial landscape
##' @param Nrun Number of run of cellular automaton
##' @return list of landscape results
##' @author Georges Kunstler
eval_coex <- function(list_int_coex, Nrun = 1000,
                       vec_p_d = seq(0.001, 0.5, length = 20),
                       vec_p_s  = seq(0.001, 0.5, length = 20),
                       vec_K = c(1, 5, 10, 100)){

list_coex <- list()
for ( K in 1:length(vec_K)){
 list_coex_K<- list()
   for (d in 1:length(vec_p_d)){
       list_coex_d<- list()
       for (s in 1:length(vec_p_s)){
           res <- UpdateIterR(list_init_coex$mat_sp, list_init_coex$mat_suc,
                            list_init_coex$c_e, list_init_coex$c_l,
                            list_init_coex$c_s,
                            list_init_coex$ss,
                            vec_p_d[d], vec_p_s[s], vec_K[K], n = Nrun)
           list_coex_d[[s]] <- res
           print("done s")
       }
       list_coex_K[[d]] <- list_coex_d
       print("done d")
   }
 list_coex[[K]] <- list_coex_K
 print("done K")
}
return(list_coex)
}


## function to process output

##' .. Compute abundance of each species. ..
##'
##' .. content for \details{} ..
##' @title table_level
##' @param x vector or matrix of species number
##' @param nsp number of species
##' @return
##' @author Georges Kunstler
table_level <- function(x, nsp){
  v <- rep(0, nsp+1)
  names(v) <- 0:nsp
  r <- table(x)
  v[names(r)] <-  r
  return(v)
}

##' .. Compute species abundance per column ..
##'
##' .. content for \details{} ..
##' @title table_level_row
##' @param mm matrix of species
##' @param nsp number of species
##' @return matrix
##' @author Georges Kunstler
table_level_row <- function(mm, nsp){
t(apply(mm, MARGIN = 1, table_level, nsp))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title mean_landscape_l
##' @param list_res list of landscape
##' @param c_l species parameters
##' @return vector
##' @author Georges Kunstler
mean_landscape_l <-  function(list_res, c_l){
    mat_l <-  return_mat_fill_c(list_res[[1]], c_l)
    apply(mat_l, MARGIN = 2, FUN = mean, na.rm = TRUE)
}
