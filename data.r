coronary.dat <- read.table('coronary.txt', header=T)
rhc.dat <- read.csv('rhc-small.txt', header=T)

observed <- table(rhc.dat)

graph <- matrix(c(0, 1, 1, 0, 0, 0,
                  1, 0, 1, 0, 0, 0,
                  1, 1, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 1,
                  0, 0, 0, 1, 1, 0), 6)

