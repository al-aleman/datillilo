# Load the required packages
library(pcadapt)
library(qvalue)
library(LEA)

# Read the data from a BED file
filename <- read.pcadapt("valida.oneSNP.bed", type = "bed")

# Run PCAdapt with 20 principal components
x <- pcadapt(input = filename, K = 20)

# Plot the screeplot
plot(x, option = "screeplot")

# Define population names
poplist.names <- c(rep("North", 54),rep("Central", 72),rep("South", 21))

# Plot PCAdapt scores for different populations
plot(x, option = "scores", pop = poplist.names)

# Run PCAdapt again with 3 principal components
x <- pcadapt(filename, K = 3)

# Plot the Manhattan plot
plot(x, option = "manhattan")

# Identify outliers based on adjusted p-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
outliers
# Expected result should be close to:
# [1]  285  458  542  646  716  758  913 1085 1501 1743 1987 2087 2235 2497 2909 3044 3920 3976 3983 4193

# Convert a VCF file to LEA format
vcf2lfmm("valida.oneSNP.vcf", output.file = "valida.oneSNP.lfmm", force = TRUE)

# Run LFMM analysis with 3 latent factors and other parameters
project <- lfmm("valida.oneSNP.lfmm", "valida.oneSNP.env", K = 3, repetitions = 10, CPU = 24, iterations = 10000, burnin = 5000, project = "new")

# Compute z-scores for each latent factor
z.one <- z.scores(project, K = 3, d = 1)
z.one <- apply(z.one, 1, median)
lambda.one <- median(z.one^2) / qchisq(0.5, df = 1)
p.one.adj <- pchisq(z.one^2 / lambda.one, df = 1, lower = FALSE)

# Repeat the above steps for other latent factors
# (from d = 2 to d = 20)
z.two = z.scores(project, K = 3, d = 2)
z.two <- apply(z.two, 1, median)
lambda.two = median(z.two^2)/qchisq(0.5, df = 1)
p.two.adj = pchisq(z.two^2/lambda.two, df = 1, lower = FALSE)
z.three = z.scores(project, K = 3, d = 3)
z.three <- apply(z.three, 1, median)
lambda.three = median(z.three^2)/qchisq(0.5, df = 1)
p.three.adj = pchisq(z.three^2/lambda.three, df = 1, lower = FALSE)
z.four = z.scores(project, K = 3, d = 4)
z.four <- apply(z.four, 1, median)
lambda.four = median(z.four^2)/qchisq(0.5, df = 1)
p.four.adj = pchisq(z.four^2/lambda.four, df = 1, lower = FALSE)
z.five = z.scores(project, K = 3, d = 5)
z.five <- apply(z.five, 1, median)
lambda.five = median(z.five^2)/qchisq(0.5, df = 1)
p.five.adj = pchisq(z.five^2/lambda.five, df = 1, lower = FALSE)
z.six = z.scores(project, K = 3, d = 6)
z.six <- apply(z.six, 1, median)
lambda.six = median(z.six^2)/qchisq(0.5, df = 1)
p.six.adj = pchisq(z.six^2/lambda.six, df = 1, lower = FALSE)
z.seven = z.scores(project, K = 3, d = 7)
z.seven <- apply(z.seven, 1, median)
lambda.seven = median(z.seven^2)/qchisq(0.5, df = 1)
p.seven.adj = pchisq(z.seven^2/lambda.seven, df = 1, lower = FALSE)
z.eight = z.scores(project, K = 3, d = 8)
z.eight <- apply(z.eight, 1, median)
lambda.eight = median(z.eight^2)/qchisq(0.5, df = 1)
p.eight.adj = pchisq(z.eight^2/lambda.eight, df = 1, lower = FALSE)
z.nine = z.scores(project, K = 3, d = 9)
z.nine <- apply(z.nine, 1, median)
lambda.nine = median(z.nine^2)/qchisq(0.5, df = 1)
p.nine.adj = pchisq(z.nine^2/lambda.nine, df = 1, lower = FALSE)
z.ten = z.scores(project, K = 3, d = 10)
z.ten <- apply(z.ten, 1, median)
lambda.ten = median(z.ten^2)/qchisq(0.5, df = 1)
p.ten.adj = pchisq(z.ten^2/lambda.ten, df = 1, lower = FALSE)
z.eleven = z.scores(project, K = 3, d = 11)
z.eleven <- apply(z.eleven, 1, median)
lambda.eleven = median(z.eleven^2)/qchisq(0.5, df = 1)
p.eleven.adj = pchisq(z.eleven^2/lambda.eleven, df = 1, lower = FALSE)
z.twelve = z.scores(project, K = 3, d = 12)
z.twelve <- apply(z.twelve, 1, median)
lambda.twelve = median(z.twelve^2)/qchisq(0.5, df = 1)
p.twelve.adj = pchisq(z.twelve^2/lambda.twelve, df = 1, lower = FALSE)
z.thirteen = z.scores(project, K = 3, d = 13)
z.thirteen <- apply(z.thirteen, 1, median)
lambda.thirteen = median(z.thirteen^2)/qchisq(0.5, df = 1)
p.thirteen.adj = pchisq(z.thirteen^2/lambda.thirteen, df = 1, lower = FALSE)
z.fourteen = z.scores(project, K = 3, d = 14)
z.fourteen <- apply(z.fourteen, 1, median)
lambda.fourteen = median(z.fourteen^2)/qchisq(0.5, df = 1)
p.fourteen.adj = pchisq(z.fourteen^2/lambda.fourteen, df = 1, lower = FALSE)
z.fifteen = z.scores(project, K = 3, d = 15)
z.fifteen <- apply(z.fifteen, 1, median)
lambda.fifteen = median(z.fifteen^2)/qchisq(0.5, df = 1)
p.fifteen.adj = pchisq(z.fifteen^2/lambda.fifteen, df = 1, lower = FALSE)
z.sixteen = z.scores(project, K = 3, d = 16)
z.sixteen <- apply(z.sixteen, 1, median)
lambda.sixteen = median(z.sixteen^2)/qchisq(0.5, df = 1)
p.sixteen.adj = pchisq(z.sixteen^2/lambda.sixteen, df = 1, lower = FALSE)
z.seventeen = z.scores(project, K = 3, d = 17)
z.seventeen <- apply(z.seventeen, 1, median)
lambda.seventeen = median(z.seventeen^2)/qchisq(0.5, df = 1)
p.seventeen.adj = pchisq(z.seventeen^2/lambda.seventeen, df = 1, lower = FALSE)
z.eighteen = z.scores(project, K = 3, d = 18)
z.eighteen <- apply(z.eighteen, 1, median)
lambda.eighteen = median(z.eighteen^2)/qchisq(0.5, df = 1)
p.eighteen.adj = pchisq(z.eighteen^2/lambda.eighteen, df = 1, lower = FALSE)
z.nineteen = z.scores(project, K = 3, d = 19)
z.nineteen <- apply(z.nineteen, 1, median)
lambda.nineteen = median(z.nineteen^2)/qchisq(0.5, df = 1)
p.nineteen.adj = pchisq(z.nineteen^2/lambda.nineteen, df = 1, lower = FALSE)
z.twenty = z.scores(project, K = 3, d = 20)
z.twenty <- apply(z.twenty, 1, median)
lambda.twenty = median(z.twenty^2)/qchisq(0.5, df = 1)
p.twenty.adj = pchisq(z.twenty^2/lambda.twenty, df = 1, lower = FALSE)

# Perform q-value adjustment for each p-value
q.one <- qvalue(p.one.adj)$qvalues
q.two<-qvalue(p.two.adj)$qvalues
q.three<-qvalue(p.three.adj)$qvalues
q.four<-qvalue(p.four.adj)$qvalues
q.five<-qvalue(p.five.adj)$qvalues
q.six<-qvalue(p.six.adj)$qvalues
q.seven<-qvalue(p.seven.adj)$qvalues
q.eight<-qvalue(p.eight.adj)$qvalues
q.nine<-qvalue(p.nine.adj)$qvalues
q.ten<-qvalue(p.ten.adj)$qvalues
q.eleven<-qvalue(p.eleven.adj)$qvalues
q.twelve<-qvalue(p.twelve.adj)$qvalues
q.thirteen<-qvalue(p.thirteen.adj)$qvalues
q.fourteen<-qvalue(p.fourteen.adj)$qvalues
q.fifteen<-qvalue(p.fifteen.adj)$qvalues
q.sixteen<-qvalue(p.sixteen.adj)$qvalues
q.seventeen<-qvalue(p.seventeen.adj)$qvalues
q.eighteen<-qvalue(p.eighteen.adj)$qvalues
q.nineteen<-qvalue(p.nineteen.adj)$qvalues
q.twenty<-qvalue(p.twenty.adj)$qvalues

# Identify significant markers based on q-value and z-score thresholds
which(q.one < 0.01 & abs(z.one) > 2)
which(q.one<0.01 & abs(z.one)>2)
which(q.two<0.01 & abs(z.two)>2)
which(q.three<0.01 & abs(z.three)>2)
which(q.four<0.01 & abs(z.four)>2)
which(q.five<0.01 & abs(z.five)>2)
which(q.six<0.01 & abs(z.six)>2)
which(q.seven<0.01 & abs(z.seven)>2)
which(q.eight<0.01 & abs(z.eight)>2)
which(q.nine<0.01 & abs(z.nine)>2)
which(q.ten<0.01 & abs(z.ten)>2)
which(q.eleven<0.01 & abs(z.eleven)>2)
which(q.twelve<0.01 & abs(z.twelve)>2)
which(q.thirteen<0.01 & abs(z.thirteen)>2)
which(q.fourteen<0.01 & abs(z.fourteen)>2)
which(q.fifteen<0.01 & abs(z.fifteen)>2)
which(q.sixteen<0.01 & abs(z.sixteen)>2)
which(q.seventeen<0.01 & abs(z.seventeen)>2)
which(q.eighteen<0.01 & abs(z.eighteen)>2)
which(q.nineteen<0.01 & abs(z.nineteen)>2)
which(q.twenty<0.01 & abs(z.twenty)>2)

# Collect the indices of significant markers for each latent factor
BIOL1<-which(q.one<0.01 & abs(z.one)>2)
BIOL2<-which(q.two<0.01 & abs(z.two)>2)
BIOL3<-which(q.three<0.01 & abs(z.three)>2)
BIOL4<-which(q.four<0.01 & abs(z.four)>2)
BIOL5<-which(q.five<0.01 & abs(z.five)>2)
BIOL6<-which(q.six<0.01 & abs(z.six)>2)
BIOL7<-which(q.seven<0.01 & abs(z.seven)>2)
BIOL8<-which(q.eight<0.01 & abs(z.eight)>2)
BIOL9<-which(q.nine<0.01 & abs(z.nine)>2)
BIOL10<-which(q.ten<0.01 & abs(z.ten)>2)
BIOL11<-which(q.eleven<0.01 & abs(z.eleven)>2)
BIOL12<-which(q.twelve<0.01 & abs(z.twelve)>2)
BIOL13<-which(q.thirteen<0.01 & abs(z.thirteen)>2)
BIOL14<-which(q.fourteen<0.01 & abs(z.fourteen)>2)
BIOL15<-which(q.fifteen<0.01 & abs(z.fifteen)>2)
BIOL16<-which(q.sixteen<0.01 & abs(z.sixteen)>2)
BIOL17<-which(q.seventeen<0.01 & abs(z.seventeen)>2)
BIOL18<-which(q.eighteen<0.01 & abs(z.eighteen)>2)
BIOL19<-which(q.nineteen<0.01 & abs(z.nineteen)>2)
BIOL20<-which(q.twenty<0.01 & abs(z.twenty)>2)

# Combine the indices of significant markers for all latent factors
joint_list <- c(BIOL1,BIOL2,BIOL3,BIOL4,BIOL5,BIOL6,BIOL7,BIOL8,BIOL9,BIOL10,BIOL11,BIOL12,BIOL13,BIOL14,BIOL15,BIOL16,BIOL17,BIOL18,BIOL19,BIOL20)
unique(joint_list)
# Expected result should be close to:
# [1]  613  899 1079 1229 1366 1394 3406 3606 3735 3789 3822 3898 3108 4365 1390 1801 1993 2430  267 2087
#[21] 2968 3301 3414 3932  557  746  870  961 1393 1828 2039  318 1563 2021  488 1179 1282 1289 1372 1391
#[41] 3066 4221  371  701 2678 3006 3284 4052 4056   51 4360  345 2030 3656 3933  734 2318 2083 3508 3882
#[61] 4218 3944 4257

lfmm.results <- cbind(z.one, q.one, z.two, q.two, z.three, q.three, z.four, q.four, z.five, q.five, z.six, q.six, z.seven, q.seven, z.eight, q.eight, z.nine, q.nine, z.ten, q.ten, z.eleven, q.eleven, z.twelve, q.twelve, z.thirteen, q.thirteen, z.fourteen, q.fourteen, z.fifteen, q.fifteen, z.sixteen, q.sixteen, z.seventeen, q.seventeen, z.eighteen, q.eighteen, z.nineteen, q.nineteen, z.twenty, q.twenty)
