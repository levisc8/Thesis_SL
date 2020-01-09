
picks <- 1:7

fact <- function(n, r) {
  factorial(n) / (factorial(r) * factorial(n - r))
}

out <- 0

for(i in picks) {
  
  out <- out + fact(7, i)
  
}

out

picks <- 1:4

out <- 0
for(i in picks) {
  out  <- out + fact(5, i)
}

out