# This function takes as input a vector of log probabilities, and returns it in a normalized way 
# so that the sum(yi in cluster k) = 1, where k =1,..,K_star

normalize_weights <- function(prob){
  const = mean(prob[prob != -Inf])
  weight = prob
  weight[prob != -Inf] = exp( prob[prob != -Inf] + const )/ sum(exp(prob[prob != -Inf] + const))
  weight[prob == -Inf] = 0
  return(weight)
}

# valutare se serve davvero escludere NA (ossia se non possono esserci)

#norm_Probs = t(sapply(1:n, function(i) normalize_weights(Probs, i)))
#rowSums(norm_Probs)