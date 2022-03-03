# Half-life (CMI paper)
# 2 microliter
posFinger <- c(18,6,5,7,5,2) # H3N2
# posFinger <- c(18,15,8,5,5,2)
time <- c(1,3,5,10,15,30)

data <- data.frame(cbind(x=time,y=posFinger))

# Produce starting values for nls
exp.model <- lm(log(y)~x, data=data)
summary(exp.model)
alpha.0 <- exp(coef(exp.model)[1])
beta.0 <- coef(exp.model)[2]
-log(2)/coef(exp.model)[2]
start <- list(alpha=18, beta=coef(model)[2])

# Fitting
model <- nls(y~alpha*exp(beta*x),data=data, start=start)
coef(model)[1]
coef(model)[2]
# half-life
-log(2)/coef(model)[2]

# Plotting
timevalues <- seq(1,30)
df.model <- data.frame(x=timevalues,y=predict(model, list(x = timevalues)))


axis.title.size <- 25
axis.text.size <- 20
ggplot(data, aes(x=x,y=y)) +
  geom_bar(stat="identity",alpha=1) +
  geom_line(data=df.model, aes(x=x,y=y),size=1.25,col="red") + 
  xlab("Time (min)") + 
  ylab("Fingers with positive virus culture (n)") +
  scale_x_continuous(breaks=time) +
  theme_bw() +
  theme(axis.title.x = element_text(size=axis.title.size),
        axis.title.y = element_text(size=axis.title.size),
        axis.text.x = element_text(size=axis.text.size),
        axis.text.y = element_text(size=axis.text.size))


# Plot with starting values
counts.exp <- exp(predict(exp.model,list(x=timevalues))) 
plot(time,posFinger,pch=16)
lines(timevalues,counts.exp, lwd=2,col="red",xlab="Time (min)", ylab="Number of positive fingers")



# 30 microliter
posFinger <- c(12,9)
time <- c(0,15)

data <- data.frame(cbind(x=time,y=posFinger))

exp.model <- lm(log(y)~x, data=data)
summary(exp.model)
(alpha.0 <- exp(coef(exp.model)[1]))
beta.0 <- coef(exp.model)[2]
-log(2)/coef(exp.model)[2]

# Plot fitted curvetimevalues <- seq(1,30)
df.model <- data.frame(x=time,y=predict(exp.model, list(x = time)))


axis.title.size <- 25
axis.text.size <- 20
# Plot fitted curve
ggplot(data, aes(x=x,y=y)) +
  geom_bar(stat="identity",alpha=1) +
  geom_line(data=df.model, aes(x=x,y=y),size=1.25,col="red") + 
  xlab("Time (min)") + 
  ylab("Fingers with positive virus culture (n)") +
  scale_x_continuous(breaks=time) +
  theme_bw() +
  theme(axis.title.x = element_text(size=axis.title.size),
        axis.title.y = element_text(size=axis.title.size),
        axis.text.x = element_text(size=axis.text.size),
        axis.text.y = element_text(size=axis.text.size))
