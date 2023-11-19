# Análise Bayesiana de Sobrevivência

Este é um script R (primeira versão) para realizar uma análise bayesiana de sobrevivência usando o pacote MCMCpack. O objetivo é estimar os parâmetros de um modelo de sobrevivência com base em dados censurados.

Ainda tem bastante coisa pra corrigir xd

## Parâmetros do Modelo

- **Valores Verdadeiros de Alpha e Lambda:** Os valores reais dos parâmetros alpha e lambda no modelo de sobrevivência.

- **Número de Observações (n):** O tamanho total da amostra.

- **Número de Falhas (r):** O número de falhas na amostra.

- **Número de Censuras (c):** O número de censuras na amostra.

- **Número de Simulações (amostra):** O número de simulações a serem executadas para estimar os parâmetros.

- **Espaçamento (espaçamento):** O espaçamento entre as amostras para reduzir a autocorrelação.

## Função de Verossimilhança

A função `lPH` implementa a log-verossimilhança para o modelo de sobrevivência. Esta função é maximizada para obter estimativas dos parâmetros alpha e lambda.

## Simulação e Ajuste do Modelo

O loop `for` simula dados de sobrevivência usando um método específico e ajusta o modelo usando MCMCmetrop1R. Se a cadeia convergir (sem NaN ou Inf), os resultados são armazenados em uma matriz `resul`.

## Estatísticas e Resultados

O script calcula estatísticas como média, desvio padrão, intervalos de confiança e probabilidade de cobertura para os parâmetros alpha e lambda. Além disso, são apresentados viés e erro quadrático médio.

## Visualizações

O script inclui gráficos de traceplot e autocorrelação para avaliar a convergência da cadeia MCMC.

Lembre-se de ajustar os valores verdadeiros de alpha e lambda e outros parâmetros conforme necessário para suas necessidades específicas.
