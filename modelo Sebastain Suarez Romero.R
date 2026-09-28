# =====================================================================
# Econometría II - Examen Parcial I (2026)
# Modelo de corte transversal: determinantes de la FECUNDIDAD
# Datos: 'swiss' (viene incluido en R, no hay que descargar nada)
#        47 provincias francófonas de Suiza, c. 1888
# =====================================================================

# ---- Paquetes (se instalan solos si faltan) --------------------------
for (p in c("lmtest", "car")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# =====================================================================
# (a) ESPECIFICACIÓN DEL MODELO
# =====================================================================
# Teoría: modelo de fecundidad de Becker (cantidad vs. calidad de hijos)
# y transición demográfica.
#
# Fertility        : índice de fecundidad general (Ig), variable dependiente
# Education        : % de reclutas con educación más allá de la primaria (-)
#                    Más capital humano -> mayor costo de oportunidad del
#                    tiempo -> menos hijos, más "calidad".
# Examination      : % de reclutas con nota más alta en el examen militar (-)
#                    Proxy de nivel/calidad educativa y modernización.
# Agriculture      : % de hombres ocupados en agricultura (+)
#                    Los hijos son mano de obra en economías agrarias.
#                    (OJO: en estos datos suele salir negativo; si pasa,
#                    coméntelo como hallazgo contrario a la teoría.)
# Catholic         : % de población católica (+)
#                    Normas culturales pro-natalistas.
# Infant.Mortality : nacidos vivos que mueren antes del año, % (+)
#                    Efecto "reemplazo/seguro": más mortalidad infantil ->
#                    más nacimientos para compensar.
#
# Modelo poblacional:
# Fertility_i = b0 + b1*Education_i + b2*Examination_i + b3*Agriculture_i
#               + b4*Catholic_i + b5*Infant.Mortality_i + u_i
#
# Signos esperados: b1<0, b2<0, b3>0, b4>0, b5>0
# Hipótesis (por coeficiente):  H0: bj = 0   vs   H1: bj != 0

data(swiss)
datos <- swiss
cat("\nObservaciones:", nrow(datos), "\n")
str(datos)

# =====================================================================
# (b) ANÁLISIS DESCRIPTIVO CON TESTS ESTADÍSTICOS
# =====================================================================
vars <- c("Fertility", "Education", "Examination",
          "Agriculture", "Catholic", "Infant.Mortality")
X <- datos[, vars]

asimetria <- function(x) { m <- mean(x); mean((x - m)^3) / (mean((x - m)^2))^1.5 }
curtosis  <- function(x) { m <- mean(x); mean((x - m)^4) / (mean((x - m)^2))^2 }

desc <- data.frame(
  Media    = sapply(X, mean),
  Mediana  = sapply(X, median),
  DesvEst  = sapply(X, sd),
  CV       = sapply(X, function(x) sd(x) / mean(x)),
  Min      = sapply(X, min),
  Max      = sapply(X, max),
  Asimetria = sapply(X, asimetria),
  Curtosis  = sapply(X, curtosis),
  ShapiroW = sapply(X, function(x) shapiro.test(x)$statistic),
  ShapiroP = sapply(X, function(x) shapiro.test(x)$p.value)
)
cat("\n--- Estadísticos descriptivos y test de normalidad (Shapiro-Wilk) ---\n")
print(round(desc, 4))
# H0 Shapiro: la variable es normal. p < 0.05 -> se rechaza normalidad.

cat("\n--- Matriz de correlaciones ---\n")
print(round(cor(X), 3))

cat("\n--- Test de correlación de Pearson de cada X con Fertility ---\n")
for (v in vars[-1]) {
  ct <- cor.test(X$Fertility, X[[v]])
  cat(sprintf("%-17s r = %6.3f  t = %6.3f  p-valor = %.4f\n",
              v, ct$estimate, ct$statistic, ct$p.value))
}

# Gráficos descriptivos
par(mfrow = c(2, 3))
for (v in vars) hist(X[[v]], main = v, xlab = "", col = "steelblue", border = "white")
par(mfrow = c(1, 1))
pairs(X, main = "Diagramas de dispersión", pch = 19, col = "steelblue")

# =====================================================================
# (c) ESTIMACIÓN DE LOS PARÁMETROS (MCO)
# =====================================================================
modelo <- lm(Fertility ~ Education + Examination + Agriculture +
               Catholic + Infant.Mortality, data = datos)
res <- summary(modelo)
cat("\n--- Resultados de la estimación MCO ---\n")
print(res)
cat("\nIntervalos de confianza al 95%:\n")
print(round(confint(modelo), 4))

# =====================================================================
# (d) TABLA ANOVA
# =====================================================================
y    <- datos$Fertility
n    <- nobs(modelo)
k    <- length(coef(modelo)) - 1          # regresores (sin constante)
gl_r <- n - k - 1
SCT  <- sum((y - mean(y))^2)
SCE  <- sum((fitted(modelo) - mean(y))^2) # suma de cuadrados explicada
SCR  <- sum(resid(modelo)^2)              # suma de cuadrados residual
CME  <- SCE / k
CMR  <- SCR / gl_r
Fcal <- CME / CMR

anova_tab <- data.frame(
  Fuente = c("Regresión", "Residual", "Total"),
  SC     = c(SCE, SCR, SCT),
  gl     = c(k, gl_r, n - 1),
  CM     = c(CME, CMR, NA),
  F      = c(Fcal, NA, NA),
  Pvalor = c(pf(Fcal, k, gl_r, lower.tail = FALSE), NA, NA)
)
cat("\n--- Tabla ANOVA ---\n")
print(anova_tab, digits = 5, row.names = FALSE, na.print = "")
cat(sprintf("\nR2 = %.4f | R2 ajustado = %.4f | Error estándar residual = %.4f\n",
            res$r.squared, res$adj.r.squared, sqrt(CMR)))

# =====================================================================
# (e) SIGNIFICANCIA INDIVIDUAL (test t)
# =====================================================================
# H0: bj = 0  vs  H1: bj != 0   | Rechazo si |t| > t crítico
alfa <- 0.05
tcrit <- qt(1 - alfa / 2, gl_r)
tab_t <- data.frame(
  Coef      = coef(modelo),
  ErrorEst  = res$coefficients[, 2],
  t_calc    = res$coefficients[, 3],
  t_crit    = tcrit,
  Pvalor    = res$coefficients[, 4],
  Decision  = ifelse(abs(res$coefficients[, 3]) > tcrit,
                     "Rechaza H0 (significativo)", "No rechaza H0")
)
signo_esp <- c(NA, "-", "-", "+", "+", "+")
tab_t$SignoEsperado <- signo_esp
tab_t$SignoObtenido <- ifelse(coef(modelo) > 0, "+", "-")
cat("\n--- Significancia individual (alfa = 5%) ---\n")
print(round_df <- within(tab_t, {
  Coef <- round(Coef, 4); ErrorEst <- round(ErrorEst, 4)
  t_calc <- round(t_calc, 3); t_crit <- round(t_crit, 3)
  Pvalor <- round(Pvalor, 4)
}))

# =====================================================================
# (f) SIGNIFICANCIA GLOBAL (test F)
# =====================================================================
# H0: b1 = b2 = b3 = b4 = b5 = 0   vs   H1: al menos un bj != 0
Fcrit <- qf(1 - alfa, k, gl_r)
cat("\n--- Significancia global ---\n")
cat(sprintf("F calculado = %.3f | F crítico (%d, %d; 5%%) = %.3f | p-valor = %.6f\n",
            Fcal, k, gl_r, Fcrit, pf(Fcal, k, gl_r, lower.tail = FALSE)))
cat(ifelse(Fcal > Fcrit,
           "Decisión: se rechaza H0 -> el modelo es globalmente significativo.\n",
           "Decisión: no se rechaza H0 -> el modelo no es globalmente significativo.\n"))

# =====================================================================
# (g) VERIFICACIÓN DE SUPUESTOS
# =====================================================================
e <- resid(modelo)

# 1. Normalidad de los residuos (H0: normales)
cat("\n=== 1. Normalidad de los errores ===\n")
print(shapiro.test(e))
S <- asimetria(e); K <- curtosis(e)
JB <- n / 6 * (S^2 + (K - 3)^2 / 4)
cat(sprintf("Jarque-Bera = %.4f | p-valor = %.4f\n", JB, 1 - pchisq(JB, 2)))

# 2. Homocedasticidad (H0: varianza constante)
cat("\n=== 2. Homocedasticidad ===\n")
print(bptest(modelo))                                   # Breusch-Pagan

# Fix for '0 (non-NA) cases' error in White's test with bptest
# Create a data frame for White's test variables explicitly
white_test_data <- data.frame(
  res_sq = resid(modelo)^2,
  fitted_vals = fitted(modelo),
  fitted_vals_sq = fitted(modelo)^2
)
# Run bptest with the new data frame for White's test (auxiliary regression on fitted values and their squares)
print(bptest(res_sq ~ fitted_vals + fitted_vals_sq, data = white_test_data))  # tipo White

print(gqtest(modelo, order.by = ~ fitted(modelo)))      # Goldfeld-Quandt

# 3. Multicolinealidad (VIF > 10 problemático; > 5 alerta)
cat("\n=== 3. Multicolinealidad (VIF) ===\n")
print(round(vif(modelo), 3))

# 4. Especificación funcional (H0: forma funcional correcta)
cat("\n=== 4. Especificación (RESET de Ramsey) ===\n")
print(resettest(modelo, power = 2:3, type = "fitted"))

# 5. Media cero de los errores
cat("\n=== 5. Media de los residuos ===\n")
cat("Media =", format(mean(e), scientific = TRUE), " | t-test p-valor =",
    round(t.test(e, mu = 0)$p.value, 4), "\n")

# 6. Observaciones influyentes
cat("\n=== 6. Observaciones influyentes ===\n")
cook <- cooks.distance(modelo)
cat("Cook > 4/n:\n"); print(round(cook[cook > 4 / n], 4))

# Gráficos de diagnóstico
par(mfrow = c(2, 2)); plot(modelo); par(mfrow = c(1, 1))
hist(e, breaks = 12, freq = FALSE, col = "grey85", main = "Histograma de residuos")
curve(dnorm(x, mean(e), sd(e)), add = TRUE, col = "red", lwd = 2)

# 7. Si hay heterocedasticidad: errores estándar robustos (HC3)
cat("\n=== Errores estándar robustos (HC3), por si falla homocedasticidad ===\n")
print(coeftest(modelo, vcov. = hccm(modelo, type = "hc3")))

# Nota: autocorrelación (Durbin-Watson) no aplica en corte transversal
# porque el orden de las observaciones no tiene significado.
