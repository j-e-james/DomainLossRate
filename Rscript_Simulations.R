library(ggplot2)


Simulations <- data.frame(
  c(10,0.35191658230411016,0.0030707601310254874,0.21580693982529456,0.04657263527675831),
  c(50,0.3459554896142433,0.0030800734258826204,0.21156315651136406,0.044758969193051924),
  c(100,0.3395920787648244,0.003094813623586921,0.2068900846561273,0.04280350712901952),
  c(250,0.32297673803859367,0.0031590278640300086,0.19429947504805872,0.037752286003951196),
  c(500,0.3010102006120367,0.0032731578154292554,0.17750658041638687,0.03150858609111922),
  c(750,0.28544075829383886,0.003428897172322091,0.16519314887249384,0.02728877643440991),
  c(1000,0.2722644409114997,0.00356943592751935,0.15505491624046172,0.0240420270503366),
  c(1500,0.25186938775510204,0.0038689274764594646,0.13541246167608126,0.01833653477717618),
  c(2000,0.23911344019728728,0.004250782666582779,0.121054051187394,0.014654083308880206)
)
  
Simulations <- data.frame(t(Simulations))
colnames(Simulations) <- c('Time', 'MeanISD', 'StdErrorISD', 'StdDev', 'VarISD')

BYTime <- Simulations$Time/1000
summary(lm(Simulations$MeanISD ~ Simulations$Time))
plot(Simulations$Time, Simulations$MeanISD, xlab = 'Age (MY)', ylab = 'ISD')

ISDplot <- ggplot(data = Simulations, aes(Time, MeanISD)) +
  geom_point(color = "darkslategray4", size=2.5)+
  geom_errorbar(aes(ymax=MeanISD+StdErrorISD, ymin =MeanISD-StdErrorISD))+
  geom_hline(yintercept=0.18, linetype = 'dashed', col = 'midnightblue')

plot <- ISDplot + labs(x = "Time, MY") +
  labs(y = "ISD") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
plot  

ISDplot <- ggplot(data = Simulations, aes(Time, StdDev)) +
  geom_point(color = "darkslategray4", size=2.5)

plot <- ISDplot + labs(x = "Time, MY") +
  labs(y = "Standard Deviation, ISD") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
plot  

Simulations <- data.frame(
  c(10,1.0467775303643725,0.006505332833897853,0.4572281871037261,0.20905761508215998),
  c(50,1.0291548965222959,0.006218587060877838,0.4257348434566201,0.1812501569330328),
  c(100,1.0135679458239277,0.005962513957918534,0.39685446583576944,0.1574934670537939),
  c(250,0.9840389920424404,0.005397598903545933,0.33141433084859534,0.10983545869182222),
  c(500,0.9559524779361847,0.004683473552473112,0.25420520770873817,0.06462028762624272),
  c(750,0.946476393024245,0.00456488009318715,0.221337850928009,0.048990444253429535),
  c(1000,0.9427757638529259,0.004648489355915341,0.2042692421951098,0.04172592330696443),
  c(1500,0.9390272373540857,0.005097646349275933,0.18273480188086183,0.033392007818437826),
  c(2000,0.9341580135440181,0.005531797404448656,0.16465811038677775,0.027112293316144287)
)

Simulations <- data.frame(t(Simulations))
colnames(Simulations) <- c('Time', 'MeanClustering', 'StdErrorClustering', 'StdDev', 'VarClustering')

BYTime <- Simulations$Time/1000
summary(lm(Simulations$MeanClustering ~ Simulations$Time))
plot(Simulations$Time, Simulations$MeanClustering, xlab = 'Age (MY)', ylab = 'Clustering')

ISDplot <- ggplot(data = Simulations, aes(Time, MeanClustering)) +
  geom_point(color = "darkslategray4", size=2.5)+
  geom_errorbar(aes(ymax=MeanClustering+StdErrorClustering, ymin =MeanClustering-StdErrorClustering))+
  geom_hline(yintercept=0.853, linetype = 'dashed', col = 'midnightblue')

plot <- ISDplot + labs(x = "Time, MY") +
  labs(y = "Clustering") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
plot  

ISDplot <- ggplot(data = Simulations, aes(Time, StdDev)) +
  geom_point(color = "darkslategray4", size=2.5)

plot <- ISDplot + labs(x = "Time, MY") +
  labs(y = "Standard Deviation, Clustering") + theme_bw() + theme_linedraw() + theme(text = element_text(size=30))
plot  
