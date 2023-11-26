rm(list = ls(all = TRUE))

######################################################
############### ESTIMAÇÃO CLÁSSICA ###################
######################################################

#--------------------------------
# INFORMAÇÕES SOBRE OS PARÂMETROS
#--------------------------------

alpha <- 2
lambda <- 1
n <- 50
r <- 50 # falhas
D <- r
c <- (n - r) # censuras

amostra <- 500 # qtd de simulaçoes

#--------------------------------
# Descrições dos vetores
#--------------------------------

X <- matrix(0, nrow = amostra, ncol = r) # Em cada LINHA a uma amostra
resul <- matrix(0, nrow = amostra, ncol = 14)

#--------------------------------
# LOG VEROSSIMILHANÇA (ECPHA)
#--------------------------------

lPH <- function(theta, x) {
    D <- length(x)
    alpha <- theta[1]
    lambda <- theta[2]

    if (theta[1] <= 0) {
        return(-Inf)
    }
    if (theta[2] <= 0) {
        return(-Inf)
    }

    lprogH <- (D * log(alpha) + D * log(lambda) + D * Rd * log(1 - exp(-lambda / T1)) + sum((alpha * R + alpha - 1) * log(1 - exp(-lambda / x))) - sum(log(x^2)) - lambda * sum(1 / x))

    if (is.na(lprogH) == TRUE) {
        return(-Inf)
    } else {
        return(lprogH)
    }
}

#--------------------
#       Simulação
## -------------------

for (j in 1:amostra) {
    R <- c(rep(1, c), rep(0, (r - c)))
    Rd <- n - D - sum(R[1:r])

    T1 <- exp(0.5)
    U <- vector()
    V <- vector()

    W <- runif(r)

    for (i in 1:r)
    {
        V[i] <- W[i]^(1 / (i + sum(R[(1 + r - i):r])))
    }

    for (i in 1:r)
    {
        U[i] <- 1 - (prod(V[(r - i + 1):r]))
    }

    Tinv <- function(u) (-lambda / log(-exp(log((1 - u)) / alpha) + 1))
    x <- sort(sapply(U, Tinv))
    Tfinal <- 1



    mlPH <- optim(x = x, c(alpha, lambda), fn = lPH, method = "BFGS", control = list(fnscale = -1), hessian = T)
    fit <- mlPH$par
    sd <- sqrt(diag(solve(-mlPH$hessian)))

    while (sd[1] == "NaN" || sd[2] == "NaN" || sd[1] == "Inf" || sd[2] == "Inf") {
        U <- vector()
        V <- vector()

        W <- runif(r)

        for (i in 1:r)
        {
            V[i] <- W[i]^(1 / (i + sum(R[(1 + r - i):r])))
        }

        for (i in 1:r)
        {
            U[i] <- 1 - (prod(V[(r - i + 1):r]))
        }

        Tinv <- function(u) (-lambda / (log(1 - u^(1 / alpha))))
        x <- sort(sapply(U, Tinv))

        mlPH <- optim(x = x, c(alpha, lambda), fn = lPH, method = "BFGS", control = list(fnscale = -1), hessian = T)
        fit <- mlPH$par
        sd <- sqrt(diag(solve(-mlPH$hessian)))
    }
    X[j, ] <- rbind(x)

    resul[j, 1] <- fit[1] ## alpha.estimado
    resul[j, 2] <- fit[2] ## lamb.estimado

    resul[j, 3] <- sd[1] ## desvio.alpha.estimado
    resul[j, 4] <- sd[2] ## desvio.lamb.estimado

    resul[j, 5] <- resul[j, 1] - 1.96 * resul[j, 3] ## IC.inf.alpha
    resul[j, 6] <- resul[j, 1] + 1.96 * resul[j, 3] ## jC.sup.alpha

    resul[j, 7] <- resul[j, 2] - 1.96 * resul[j, 4] ## jC.jnf.lamb
    resul[j, 8] <- resul[j, 2] + 1.96 * resul[j, 4] ## jC.sup.lamb


    #---------------------------------------
    # Calculando a prop. de cobertura e viés
    #---------------------------------------

    if (alpha < resul[j, 5] | alpha > resul[j, 6]) resul[j, 9] <- 0 else resul[j, 9] <- 1
    if (lambda < resul[j, 7] | lambda > resul[j, 8]) resul[j, 10] <- 0 else resul[j, 10] <- 1

    resul[j, 11] <- (resul[j, 1] - alpha)
    resul[j, 12] <- (resul[j, 2] - lambda)

    resul[j, 13] <- (alpha - resul[j, 1])^2 + (resul[j, 11])^2
    resul[j, 14] <- (lambda - resul[j, 2])^2 + +(resul[j, 12])^2
}

resultadosECPHA <- matrix(apply(resul, 2, mean), ncol(resul), 1)
rownames(resultadosECPHA) <- c(
    "Alpha", "Lamb", "sd.Alpha", "sd.lamb", "IC.inf.Alpha", "IC.Sup.Alpha",
    "IC.inf.lamb", "IC.sup.lamb", "Prob. de Cobertura de Alpha", "Prob.Corbertura de Lamb",
    "Vies Alpha", "Vies Lambda", "E.Q.M alpha", "E.Q.M lamb"
)
colnames(resultadosECPHA) <- c("Resultados")
resultadosECPHA

######################################################
####################### BAYESIANA ####################
######################################################

## install.packages("MCMCpack")
require(MCMCpack)

# Defina os valores iniciais e parâmetros MCMC
theta <- c(alpha, lambda) # Valores iniciais
alpha_prior <- 0.01
lambda_prior <- 0.01


#-----------------------------------------------------
####################### PRIORIS ######################
#-----------------------------------------------------


log_prior_alpha <- function(alpha, alpha_prior) {
    if (alpha <= 0) {
        return(-Inf) # Probabilidade zero para alpha negativo ou zero
    }

    log_prior_val <- dgamma(alpha, shape = alpha_prior, rate = alpha_prior, log = TRUE)
    return(log_prior_val)
}

log_prior_lambda <- function(lambda, lambda_prior) {
    if (lambda <= 0) {
        return(-Inf) # Probabilidade zero para lambda negativo ou zero
    }

    log_prior_val <- dgamma(lambda, shape = lambda_prior, rate = lambda_prior, log = TRUE)
    return(log_prior_val)
}


#-----------------------------------------------------
####################### POSTERIORIS ######################
#-----------------------------------------------------
# Valores iniciais

# Função log posteriori
log_posteriori <- function(theta, x, alpha_prior, lambda_prior) {
    alpha <- theta[1]
    lambda <- theta[2]

    # Log-priori para alpha
    log_prior_alpha_val <- log_prior_alpha(alpha, alpha_prior)

    # Log-priori para lambda
    log_prior_lambda_val <- log_prior_lambda(lambda, lambda_prior)

    # Log-verossimilhança
    log_likelihood <- lPH(theta, x)

    # Log-posteriori é a soma das contribuições
    log_posteriori_val <- log_prior_alpha_val + log_prior_lambda_val + log_likelihood
    return(log_posteriori_val)
}


posteriori <- function(theta, x, alpha_prior, lambda_prior) {
    log_posteriori_val <- log_posteriori(theta, x, alpha_prior, lambda_prior)

    # Exponenciação para obter a posteriori
    posteriori_val <- exp(log_posteriori_val)

    return(posteriori_val)
}


#--------------------------------
###### AJUSTE DO MODELO #########
#--------------------------------
## Vetores

a.bayes <- double(amostra)
b.bayes <- double(amostra)
sd.a.bayes <- double(amostra)
sd.b.bayes <- double(amostra)
ICred.ai <- double(amostra)
ICred.as <- double(amostra)
ICred.bi <- double(amostra)
ICred.bs <- double(amostra)
probc_a <- rep(1, times = amostra)
probc_b <- rep(1, times = amostra)
Vicio_a <- double(amostra)
Vicio_b <- double(amostra)
EQM_a <- double(amostra)
EQM_b <- double(amostra)


#--------------------------------
####### Início do Loop ##########
#--------------------------------

inicio <- Sys.time() # Inicia o contador de tempo

for (v in 1:amostra) {
    x <- X[v, ]

    # Valores iniciais e parâmetros MCMC

    theta_init <- c(alpha, lambda)
    alpha_prior <- 0.01
    lambda_prior <- 0.01

    # Número de iterações

    num_iter <- 10000

    results <- matrix(0, nrow = num_iter, ncol = 4)

    # Resultados do MCMC

    results <- MCMCmetrop1R(
        logfun = FALSE, fun = posteriori, theta.init = theta_init, x = X[v, ], alpha_prior = alpha_prior, lambda_prior = lambda_prior,
        mcmc = num_iter, thin = 1, burnin = 2000
    )
    # ---------------------------------------
    ############### Estatísticas ############
    #----------------------------------------

    a.b <- mean(results[, 1])
    b.b <- mean(results[, 2])
    a.b
    b.b

    sd.a <- sd(results[, 1])
    sd.b <- sd(results[, 2])
    sd.a
    sd.b

    a.bayes[v] <- a.b
    b.bayes[v] <- b.b
    sd.a.bayes[v] <- sd.a
    sd.b.bayes[v] <- sd.b


    ICred.ai[v] <- quantile(results[, 1], 0.025)
    ICred.as[v] <- quantile(results[, 1], 0.975)
    ICred.bi[v] <- quantile(results[, 2], 0.025)
    ICred.bs[v] <- quantile(results[, 2], 0.975)

    probc_a[v] <- mean(alpha > ICred.ai[v] & alpha < ICred.as[v])
    probc_b[v] <- mean(lambda > ICred.bi[v] & lambda < ICred.bs[v])

    # Imprime o contador do loop e o tempo decorrido
    tempo_decorrido <- as.numeric(Sys.time() - inicio, units = "secs")
    minutos <- floor(tempo_decorrido / 60)
    segundos <- round(tempo_decorrido %% 60)
    cat("Tempo decorrido: ", minutos, " minuto(s) e ", segundos, " segundo(s)\n")
    cat("Executando a amostra número: ", v, "\n")
}

# ---------------------------------------
# Calculando a prop. de cobertura e viés
#----------------------------------------

Vicio_a <- mean(a.bayes - alpha)
Vicio_b <- mean(b.bayes - lambda)
EQM_a <- mean((a.bayes - alpha)^2)
EQM_b <- mean((b.bayes - lambda)^2)

# ---------------------------------------
####### Imprimindo os resultados ########
#----------------------------------------

Bayes <- rbind(
    mean(a.bayes), mean(b.bayes), mean(sd.a.bayes), mean(sd.b.bayes), mean(ICred.ai), mean(ICred.as), mean(ICred.bi), mean(ICred.bs), mean(probc_a), mean(probc_b),
    mean(Vicio_a), mean(Vicio_b), mean(EQM_a), mean(EQM_b)
)
colnames(resultadosECPHA) <- c("Clássico")
colnames(Bayes) <- c("Bayesiano")
cbind(Bayes, resultadosECPHA)
