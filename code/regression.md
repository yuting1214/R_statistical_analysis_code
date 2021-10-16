# Formula 

## Encoding variables
```
The models fit by, e.g., the lm and glm functions are specified in a compact symbolic form. The ~ operator is basic in the formation of such models. An expression of the form y ~ model is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model. Such a model consists of a series of terms separated by + operators. The terms themselves consist of variable and factor names separated by : operators. Such a term is interpreted as the interaction of all the variables and factors appearing in the term.

In addition to + and :, a number of other operators are useful in model formulae. The * operator denotes factor crossing: a*b interpreted as a+b+a:b. The ^ operator indicates crossing to the specified degree. For example (a+b+c)^2 is identical to (a+b+c)*(a+b+c) which in turn expands to a formula containing the main effects for a, b and c together with their second-order interactions. The %in% operator indicates that the terms on its left are nested within those on the right. For example a + b %in% a expands to the formula a + a:b. The - operator removes the specified terms, so that (a+b+c)^2 - a:b is identical to a + b + c + b:c + a:c. It can also used to remove the intercept term: when fitting a linear model y ~ x - 1 specifies a line through the origin. A model with no intercept can be also specified as y ~ x + 0 or y ~ 0 + x.
```
# Model Diagnosis
http://www.sthda.com/english/articles/39-regression-model-diagnostics/161-linear-regression-assumptions-and-diagnostics-in-r-essentials/

https://daviddalpiaz.github.io/appliedstats/model-diagnostics.html

# Model Selection

## 1. Best Subset Selection
http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/155-best-subsets-regression-essentials-in-r/
http://www.science.smith.edu/~jcrouser/SDS293/labs/lab8-r.html
http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/
```

## 2. PRESS for model selection
https://stevencarlislewalker.wordpress.com/2013/06/18/calculating-the-press-statistic-in-r/
https://statisticaloddsandends.wordpress.com/2018/07/30/the-press-statistic-for-linear-regression/


library(leaps) #regsubsets()
```
## 2. Stepwise Regression
```
library(leaps) #regsubsets()
library(MASS) #stepAIC()
library(stats) #step(), base function
library(olsrr)

```

## 3. Type I, II, III SS
https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
http://www.utstat.utoronto.ca/reid/sta442f/2009/typeSS.pdf
```
library(car)
Anova(lm(Y ~ . - batch, data = data_final), type=3)
```

# Multicollinearity Diagnostics

## Varaince Inflation Factor (VIF)



 
