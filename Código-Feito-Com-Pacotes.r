
rm(list=ls(all=TRUE))
# Carregar o pacote
library(MCMCpack)
#--------------------------------
#INFORMAÇÕES SOBRE OS PARÂMETROS
#--------------------------------


# Definir os valores verdadeiros de Alpha e Lambda
alpha<-2 ; lambda<-1 # valor verdadeiro

n=1000
r=1000 #falhas
D=r
c= (n-r) #censuras
amostra<- 2000 #qtd de simulaçoes
espaçamento <- 2
#--------------------------------
# Descrições dos vetores
#--------------------------------

X<-matrix(0,nrow=amostra,ncol=r)#Em cada LINHA a uma amostra
resul<- matrix(0,nrow=amostra,ncol=14)

#--------------------------------

#----------------------------
#LOG VEROSSIMILHANÇA (ECPHA)
#----------------------------

lPH=function(theta,x)
    
{
    
    
    D=length(x)
    alpha=theta[1]
    lambda=theta[2]
    
    if(theta[1]<=0)return(-Inf)
    if(theta[2]<=0)return(-Inf)
    
    lprogH= (D*log(alpha) + D*log(lambda) + D*Rd*log(1-exp(-lambda/T1))+ 
                 sum((alpha*R + alpha -1)*log(1-exp(-lambda/x)))- sum(log(x^2)) -
                 lambda*sum(1/x))
    
    if(is.na(lprogH)==TRUE)
        
    {
        return(-Inf)
    }
    
    else
    {   
        return(lprogH) 
    }
    
}

#Simulação
##-------------------

for(j in 1:amostra){
    
    R=c(rep(1,c),rep(0,(r-c)))
    Rd= n - D - sum(R[1:r]); 
    
    T1<- exp(0.5); 
    U<- vector() ; V<-vector()
    
    W = runif(r);
    
    for(i in 1:r)
    { 
        V[i] = W[i]^(1/(i+sum(R[(1+r-i):r])))
    } ;
    
    for(i in 1:r)
    {U[i]= 1 - (prod(V[(r-i+1):r]))};
    
    Tinv= function(u)(-lambda/(log(1-u^(1/alpha))))
    x=sort(sapply(U,Tinv))}

# Definir a função de log-posterior
# Definir x como uma variável global

# Definir a função de log-posterior
logpost <- function(theta) {
  lPH(theta, x)
}

# Ajustar o modelo usando MCMC
start <- c(alpha = 2, lambda = 2)  # valores iniciais para alpha e lambda
#fit <- MCMCmetrop1R(burnin = 500,fun = logpost, theta.init = start, thin = 10, mcmc = 1000)
fit <- MCMCmetrop1R(fun = logpost, theta.init = start, mcmc = amostra, burnin = 1000,thin = espaçamento)

# Calcular as estatísticas
Alpha <- fit[,1]
Lamb <- fit[,2]
sd.Alpha <- sd(Alpha)
sd.lamb <- sd(Lamb)
IC.inf.Alpha <- quantile(Alpha, 0.025)
IC.Sup.Alpha <- quantile(Alpha, 0.975)
IC.inf.lamb <- quantile(Lamb, 0.025)
IC.sup.lamb <- quantile(Lamb, 0.975)
Prob.Cobertura.Alpha <- mean(Alpha > IC.inf.Alpha & Alpha < IC.Sup.Alpha)
Prob.Cobertura.Lamb <- mean(Lamb > IC.inf.lamb & Lamb < IC.sup.lamb)

# Calcular o viés e o EQM

Vies.Alpha <- mean(Alpha) - alpha  # você precisa fornecer o valor verdadeiro de Alpha
Vies.Lamb <- mean(Lamb) - lambda  # você precisa fornecer o valor verdadeiro de Lamb
EQM.Alpha <- mean((Alpha - alpha)^2)  # você precisa fornecer o valor verdadeiro de Alpha
EQM.Lamb <- mean((Lamb - lambda )^2)  # você precisa fornecer o valor verdadeiro de Lamb

# Exibir os resultados
cat("Alpha: ", mean(Alpha), "\n",
    "Lambda: ", mean(Lamb), "\n",
    "Desvio padrão de Alpha: ", sd(Alpha), "\n",
    "Desvio padrão de Lambda: ", sd(Lamb), "\n",
    "IC inferior de Alpha: ", quantile(Alpha, 0.025), "\n",
    "IC superior de Alpha: ", quantile(Alpha, 0.975), "\n",
    "IC inferior de Lambda: ", quantile(Lamb, 0.025), "\n",
    "IC superior de Lambda: ", quantile(Lamb, 0.975), "\n",
    "Probabilidade de cobertura de Alpha: ", mean(Alpha > IC.inf.Alpha & Alpha < IC.Sup.Alpha), "\n",
    "Probabilidade de cobertura de Lambda: ", mean(Lamb > IC.inf.lamb & Lamb < IC.sup.lamb), "\n",
    "Viés de Alpha: ", Vies.Alpha, "\n",
    "Viés de Lambda: ", Vies.Lamb, "\n",
    "EQM de Alpha: ", EQM.Alpha, "\n",
    "EQM de Lambda: ", EQM.Lamb, "\n")

length(Alpha)

# Traceplots
plot(fit)
# Autocorrelation plots
acf(fit)

