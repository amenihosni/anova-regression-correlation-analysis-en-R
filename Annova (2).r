# Gas Turbine CO and NOx Emission Dataset Analysis

# Installation des librairies nécessaires:
install.packages(c("tidyverse", "ggplot2", "dplyr", "caret", "MASS", "reshape2", 
                 "GGally","robustbase", "ggcorrplot", "car", "lmtest", "nortest", "moments", "boot"))

# Charger les bibliothèques:
library(tidyverse)
library(robustbase)
library(ggplot2)
library(dplyr)
library(caret)
library(MASS)
library(reshape2)
library(GGally)
library(ggcorrplot)
library(car)
library(lmtest)
library(nortest)
library(moments)
library(boot)
# Charger les jeux de données:
gt_2011 <- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")
gt_2012 <- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")
gt_2013 <- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")
gt_2014 <- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")
gt_2015 <- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")
ai4i2020<- read.table(file = file.choose(), header = TRUE, sep = ",", dec = ".", na.strings = "")


# Combiner les jeux de données en un seul:
gt_combined <- rbind(gt_2011, gt_2012, gt_2013, gt_2014, gt_2015)

# Creation d'un fichier pour stocker les graphes:
plot_path <- "C:\\Users\\ameni\\OneDrive\\Desktop\\stat\\projet stat\\plots"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# Exploration initiale des données:
summary(gt_combined)
str(gt_combined)

summary(ai4i2020)
str(ai4i2020)

#------------------------Préparation des données-------------------------------------

# 1. Vérifier les valeurs manquantes dans chaque colonne:
sapply(gt_combined, function(x) sum(is.na(x))) # Pas de valeurs manquantes
sapply(ai4i2020, function(x) sum(is.na(x)))

# 2. Traitement des valeurs aberrantes (Capping des outliers):
# Boxplot pour chaque colonne numérique avant le traitement des valeurs aberrantes
# Vérification de la structure du jeu de données
str(gt_combined)
str(gt_2011)
str(gt_2012)
names(gt_2011)
names(gt_2012)
gt_combined <- rbind(gt_2011, gt_2012, gt_2013, gt_2014, gt_2015)
str(gt_combined)

# Identification des colonnes numériques
numeric_columns <- names(gt_combined)[vapply(gt_combined, is.numeric, logical(1))]


# Afficher les colonnes numériques pour validation
print(numeric_columns)

# Boxplot pour chaque colonne numérique avant le traitement des valeurs aberrantes
for (col in numeric_columns) {
  png(filename = paste0(plot_path, "/Boxplot_Before_", col, ".png"))
  boxplot(
    gt_combined[[col]], 
    main = paste("Boxplot for", col, "Before Handling Outliers (Capping)"), 
    ylab = col, 
    col = "lightblue", 
    outcol = "red"
  )
  dev.off()
}

# Définir une fonction pour capper les valeurs aberrantes basées sur l'IQR
cap_outliers <- function(data, column) {
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  data[[column]] <- ifelse(data[[column]] < lower_bound, lower_bound,
                           ifelse(data[[column]] > upper_bound, upper_bound, data[[column]]))
  return(data)
}

# Appliquer le capping à toutes les colonnes numériques
for (col in numeric_columns) {
  gt_combined <- cap_outliers(gt_combined, col)
}

# Boxplot pour chaque colonne numérique après le capping
for (col in numeric_columns) {
  png(filename = paste0(plot_path, "Boxplot_After_", col, ".png"))
  boxplot(
    gt_combined[[col]], 
    main = paste("Boxplot for", col, "After Handling Outliers (Capping)"), 
    ylab = col, 
    col = "lightblue", 
    outcol = "red"
  )
  dev.off()
}

# 3. Suppression des caractéristiques fortement corrélées

# Charger les bibliothèques nécessaires
library(ggcorrplot)
library(gridExtra)

# Calculer la matrice de corrélation avant la suppression
cor_matrix <- cor(gt_combined[, numeric_columns], use = "complete.obs")

# Tracer la matrice de corrélation avant suppression
plot_before <- ggcorrplot(cor_matrix, hc.order = TRUE, type = "upper", lab = TRUE, 
                          title = "Correlation Matrix Before")

# Identifier les caractéristiques fortement corrélées
library(caret)
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.8)

# Supprimer les colonnes fortement corrélées dans gt_combined
gt_combined <- gt_combined[, -highly_correlated]

# Calculer la matrice de corrélation après suppression
cor_matrix_reduced <- cor(gt_combined[, sapply(gt_combined, is.numeric)], use = "complete.obs")

# Tracer la matrice de corrélation après suppression
plot_after <- ggcorrplot(cor_matrix_reduced, hc.order = TRUE, type = "upper", lab = TRUE, 
                         title = "Correlation Matrix After")

# Combiner les deux graphiques dans une seule figure
png(filename = paste0(plot_path, "Correlation_Matrices_Before_After.png"), width = 1200, height = 600)
grid.arrange(plot_before, plot_after, ncol = 2)
dev.off()





# 4. Transformation de la variable cible (CO)
# Appliquer une transformation logarithmique pour résoudre les problèmes de non-normalité
if ("CO" %in% names(gt_combined)) {
  gt_combined$CO <- log1p(gt_combined$CO)  # Transformation log1p pour éviter des valeurs infinies pour 0
}

# 5. Test de normalité avec Kolmogorov-Smirnov
ks_test <- ks.test(gt_combined$CO, "pnorm", mean = mean(gt_combined$CO), sd = sd(gt_combined$CO))
print(ks_test)

# Si p-value < 0.05, la distribution n'est pas normale et nécessite une transformation.
if (ks_test$p.value < 0.05) {
  # Appliquer une transformation Box-Cox pour rendre les données plus normales
  # S'assurer que CO est positif avant la transformation Box-Cox
  if (any(gt_combined$CO <= 0)) {
    gt_combined$CO <- gt_combined$CO + abs(min(gt_combined$CO)) + 1
  }
  
  # Appliquez la transformation Box-Cox
  boxcox_trans <- boxcox(lm(CO ~ ., data = gt_combined))
  lambda <- boxcox_trans$x[which.max(boxcox_trans$y)]  # Trouver la valeur lambda optimale
  
  # Appliquer la transformation Box-Cox en fonction du lambda calculé
  gt_combined$CO <- if (lambda == 0) {
    log(gt_combined$CO)
  } else {
    (gt_combined$CO^lambda - 1) / lambda
  }
}



# 6. Ajouter des termes quadratiques ou interactions
gt_combined$TIT_squared <- gt_combined$TIT^2
gt_combined$AT_TIT <- gt_combined$AT * gt_combined$TIT

# 7. Division des données en ensembles d'entraînement et de test pour la régression linéaire
set.seed(123)
train_indices <- sample(1:nrow(gt_combined), 0.8 * nrow(gt_combined))
train_data <- gt_combined[train_indices, ]
test_data <- gt_combined[-train_indices, ]

# 8. Standardiser le jeu de données en utilisant les statistiques des données d'entraînement
numeric_columns <- names(gt_combined)[sapply(gt_combined, is.numeric)]

# Calculer la moyenne et l'écart type à partir des données d'entraînement
means <- sapply(train_data[, numeric_columns], mean)
sds <- sapply(train_data[, numeric_columns], sd)

# Appliquer la standardisation aux ensembles d'entraînement et de test
train_data[, numeric_columns] <- scale(train_data[, numeric_columns], center = means, scale = sds)
test_data[, numeric_columns] <- scale(test_data[, numeric_columns], center = means, scale = sds)

# Vérifier les données standardisées
summary(train_data)

# ----------------- Analyse statistique exploratoire et tests -----------------------------------
# 1. Analyse de corrélation:
# 1.1 Corrélation de Pearson:
correlation_matrix <- cor(train_data[, numeric_columns], use = "complete.obs")
print("Matrice de corrélation de Pearson :")
print(correlation_matrix)
png(filename = paste0(plot_path, "pearson_correlation.png"))
heatmap(correlation_matrix, 
        main = "Matrice de corrélation de Pearson", 
        col = colorRampPalette(c("white", "lightblue"))(100),
        scale = "none",  # Pas de mise à l'échelle des valeurs
        cexRow = 0.8, cexCol = 0.8)
dev.off()

# 1.2 Corrélation de Spearman:
spearman_correlation <- cor(train_data[, numeric_columns], method = "spearman")
print("Matrice de corrélation de Spearman :")
print(spearman_correlation)
png(filename = paste0(plot_path, "spearman_correlation.png"))
heatmap(spearman_correlation, 
        main = "Matrice de corrélation de Spearman", 
        col = colorRampPalette(c("white", "lightblue"))(100),
        scale = "none",  # Pas de mise à l'échelle des valeurs
        cexRow = 0.8, cexCol = 0.8)
dev.off()

# 1.3 Nuages de points pour CO par rapport aux autres variables :
for (var in numeric_columns) {
  if (var != "CO") {
    # Créez le graphique pour chaque variable
    plot <- ggplot(gt_combined, aes(x = .data[[var]], y = CO)) +
      geom_point() +
      geom_smooth(method = "lm", col = "red") +
      ggtitle(paste("Scatterplot of CO vs", var)) +
      xlab(var) +
      ylab("CO") +
      theme_minimal()
    
    # Enregistrer le graphique dans un fichier PNG avec un fond blanc
    ggsave(filename = paste0("scatterplot_CO_vs_", var, ".png"),
           plot = plot,
           path = plot_path,
           width = 8, height = 6,
           bg = "white")  # Définit un fond blanc
  }
}
# ----------------- Tests statistiques -----------------------------------
# 2.1 T-Tests
# Test t indépendant pour CO entre deux groupes
gt_combined$FailureCategory <- ifelse(gt_combined$CO > mean(gt_combined$CO), "High", "Low")
t_test_result <- t.test(CO ~ FailureCategory, data = gt_combined)
print(t_test_result)

# Graphique pour le test t (Boxplot pour CO par catégorie de défaillance)
ggplot(gt_combined, aes(x = FailureCategory, y = CO)) + 
  geom_boxplot(fill = "lightblue", color = "black") +
  ggtitle("Boxplot de CO par catégorie de défaillance") +
  xlab("Catégorie de défaillance") +
  ylab("CO")

# 2.2 Test du Chi-Carré
contingency_table <- table(gt_combined$FailureCategory, gt_combined$NOX > median(gt_combined$NOX))
chi_sq_test <- chisq.test(contingency_table)
print(chi_sq_test)

# Graphique pour le test du Chi-Carré (Barplot pour Catégorie de défaillance et NOX > médiane)
ggplot(as.data.frame(contingency_table), aes(x = Var1, fill = Var2)) +
  geom_bar(position = "dodge") +
  ggtitle("Test du Chi-Carré : Catégorie de défaillance vs NOX > médiane") +
  xlab("Catégorie de défaillance") +
  ylab("Fréquence") +
  scale_fill_manual(values = c("lightblue", "lightgreen"))

# 2.3 Test de Kruskal-Wallis
kruskal_test <- kruskal.test(CO ~ FailureCategory, data = gt_combined)
print(kruskal_test)

# Graphique pour le test de Kruskal-Wallis (Boxplot pour CO par catégorie de défaillance)
ggplot(gt_combined, aes(x = FailureCategory, y = CO)) + 
  geom_boxplot(fill = "lightblue", color = "black") +
  ggtitle("Boxplot de CO par catégorie de défaillance (Test de Kruskal-Wallis)") +
  xlab("Catégorie de défaillance") +
  ylab("CO")


--------------------------Regression Linéaire---------------------------------
# Modèle de régression par étapes.
stepwise_model <- step(
  lm(CO ~ GTEP + AT + AP + AH + AFDP + TAT + TEY + NOX, data = train_data), 
  direction = "both",
  trace = FALSE  # Suppress detailed output
)
summary(stepwise_model)

# Diagnostics des résidus pour le modèle par étapes
Residuals_stepwise <- resid(stepwise_model)

# Tracer les résidus par rapport aux valeurs ajustées avec couleurs
png(filename = paste0(plot_path, "Stepwise_Fitted_vs_Residuals.png"))
plot(stepwise_model$fitted.values, Residuals_stepwise, 
     main = "Residuals vs Fitted Values for Stepwise Model", 
     xlab = "Fitted Values", ylab = "Residuals",
     col = "lightblue", pch = 16)  # Points colorés en bleu clair
abline(h = 0, col = "red")  # Ligne de référence en rouge
dev.off()

# Tracer le graphique QQ pour vérifier la normalité des résidus avec couleurs
png(filename = paste0(plot_path, "Stepwise_QQ_Plot.png"))
qqnorm(Residuals_stepwise, col = "lightblue")  # Points colorés en bleu clair
qqline(Residuals_stepwise, col = "red")  # Ligne QQ en rouge
dev.off()

# Test de Breusch-Pagan pour l'hétéroscédasticité
#bptest_stepwise <- bptest(stepwise_model)
#print(bptest_stepwise)

#---------------------------------- ANOVA Analysis--------------------------------------

# 1. Analyse ANOVA pour le modèle à sélection pas à pas (Stepwise Model)
anova_results <- anova(stepwise_model)
print(anova_results)  # Les résultats montrent que toutes les variables sélectionnées ont un effet significatif sur CO.

# Graphique de l'ANOVA pour le modèle pas à pas
png(filename = paste0(plot_path, "Stepwise_ANOVA.png"))
boxplot(CO ~ fitted(stepwise_model), data = train_data, 
        main = "ANOVA: Modèle à sélection pas à pas", 
        xlab = "Valeurs ajustées", ylab = "CO", 
        col = "lightblue")
dev.off()

# 4. Comparer plusieurs modèles à l'aide de l'ANOVA
anova_results_multiple <- anova(anova_model_1, anova_model_2, anova_model_3, anova_model_4, anova_model_5)
print(anova_results_multiple)  # Display ANOVA results for multiple models.

# Graphique de l'ANOVA pour chaque modèle comparé
png(filename = paste0(plot_path, "Multiple_Model_ANOVA.png"))
par(mfrow = c(1, 5))  # Plusieurs graphiques sur une ligne
models <- list(anova_model_1, anova_model_2, anova_model_3, anova_model_4, anova_model_5)
for (i in 1:5) {
  boxplot(CO ~ fitted(models[[i]]), data = train_data,
          main = paste("Modèle", i), xlab = "Valeurs ajustées", ylab = "CO",
          col = "lightblue")
}
dev.off()

# 5. Analyse ANOVA pour chaque variable prédictive
variables <- c("AP", "AH", "TAT", "AFDP")
for (var in variables) {
  # Check if variable can be used in Bartlett's test
  if (length(unique(gt_combined[[var]])) > 10) {
    print(paste("Variable", var, "has too many unique values. Binning it for Bartlett's test."))
    gt_combined[[var]] <- cut(gt_combined[[var]], breaks = 10)  # Bin the variable
  }
  
  # Check for sufficient group sizes
  if (all(table(gt_combined[[var]]) >= 2)) {
    # Bartlett test
    bartlett_test <- bartlett.test(as.formula(paste("CO ~ as.factor(", var, ")")), data = gt_combined)
    print(paste("Bartlett test for", var, ":"))
    print(bartlett_test)
    
    # Graphique pour chaque test de Bartlett
    png(filename = paste0(plot_path, var, "_Bartlett_Test.png"))
    boxplot(CO ~ as.factor(gt_combined[[var]]), data = gt_combined, 
            main = paste("Bartlett Test for", var), 
            xlab = var, ylab = "CO", col = "lightblue")
    dev.off()
  } else {
    print(paste("Skipped Bartlett test for", var, "due to insufficient group sizes."))
  }
}

