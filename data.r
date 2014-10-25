coronary.dat <- read.table('coronary.txt', header=T)

observed <- table(coronary.dat)

graph <- matrix(c(0, 1, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0), 6)

