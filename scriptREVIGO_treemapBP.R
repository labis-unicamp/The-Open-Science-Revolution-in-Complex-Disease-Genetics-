# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0000278","mitotic cell cycle",1.0859809163429945,0.9871819545762662,0.01778156,"mitotic cell cycle"),
c("GO:0001678","intracellular glucose homeostasis",0.05154073064080684,0.9884614944550029,-0,"intracellular glucose homeostasis"),
c("GO:0005975","carbohydrate metabolic process",5.177832504046095,0.9726998201694254,0.09347189,"carbohydrate metabolic process"),
c("GO:0006520","amino acid metabolic process",5.5080894579711925,0.972484408206426,0.12609755,"carbohydrate metabolic process"),
c("GO:0006629","lipid metabolic process",5.980051423570241,0.9721925168465625,0.12927537,"carbohydrate metabolic process"),
c("GO:0006383","transcription by RNA polymerase III",0.17492343100195132,0.9341340719816013,0.09476263,"transcription by RNA polymerase III"),
c("GO:0006457","protein folding",1.1925033490398829,0.9287984365662165,0.32896472,"transcription by RNA polymerase III"),
c("GO:0006662","glycerol ether metabolic process",0.009180846479952699,0.9716936004001495,0.02665081,"glycerol ether metabolic process"),
c("GO:0006760","folic acid-containing compound metabolic process",0.43430572691358005,0.95344666726222,0.03441951,"folic acid-containing compound metabolic process"),
c("GO:0006915","apoptotic process",0.3770497347375319,0.982919860382217,0.01472057,"apoptotic process"),
c("GO:0007049","cell cycle",2.051845347681118,0.9863699772166102,0.01766711,"cell cycle"),
c("GO:0007155","cell adhesion",1.0284049482842672,0.9872467520127078,-0,"cell adhesion"),
c("GO:0007160","cell-matrix adhesion",0.09374554956620869,0.9806988646638687,0.01294664,"cell-matrix adhesion"),
c("GO:0007338","single fertilization",0.065043712889874,0.992050557448283,-0,"single fertilization"),
c("GO:0008152","metabolic process",58.08643542977757,1,-0,"metabolic process"),
c("GO:0008202","steroid metabolic process",0.41840538689736173,0.9483966462945939,0.04171815,"steroid metabolic process"),
c("GO:0036151","phosphatidylcholine acyl-chain remodeling",0.0018386306489342266,0.9356816253299357,0.32335008,"steroid metabolic process"),
c("GO:0015811","L-cystine transport",0.015937260310373652,0.9592368795271007,-0,"L-cystine transport"),
c("GO:0006621","protein retention in ER lumen",0.03276060768047464,0.9399029693471893,0.39844138,"L-cystine transport"),
c("GO:0006897","endocytosis",0.5813395677106135,0.9268405547409807,0.31233359,"L-cystine transport"),
c("GO:0008645","hexose transmembrane transport",0.09563340727401666,0.9462297530894285,0.21047889,"L-cystine transport"),
c("GO:0015914","phospholipid transport",0.22272290516835383,0.9500377150946,0.18508379,"L-cystine transport"),
c("GO:0016192","vesicle-mediated transport",2.591260690701909,0.9345637261518414,0.28163132,"L-cystine transport"),
c("GO:0016049","cell growth",0.04912368205010616,0.9847496849729558,0.01226061,"cell growth"),
c("GO:0016477","cell migration",0.48278207313949656,0.9787423822106434,0.01508773,"cell migration"),
c("GO:0019048","modulation by virus of host process",0.006549660183151241,0.9890763233383977,-0,"modulation by virus of host process"),
c("GO:0030198","extracellular matrix organization",0.23021526352893723,0.9499871461378534,0.0140387,"extracellular matrix organization"),
c("GO:0000722","telomere maintenance via recombination",0.014509675603035163,0.9316786697194065,0.36745757,"extracellular matrix organization"),
c("GO:0007005","mitochondrion organization",0.6594013763196482,0.9468967934969572,0.3313676,"extracellular matrix organization"),
c("GO:0048790","maintenance of presynaptic active zone structure",0.011085933658366477,0.948739200594785,0.24168573,"extracellular matrix organization"),
c("GO:0051259","protein complex oligomerization",0.23557116753440024,0.9489635091646754,0.30260809,"extracellular matrix organization"),
c("GO:0061024","membrane organization",1.185564795091629,0.9467716010200933,0.38799265,"extracellular matrix organization"),
c("GO:0030574","collagen catabolic process",0.04133842268922402,0.9731354719885076,-0,"collagen catabolic process"),
c("GO:0006511","ubiquitin-dependent protein catabolic process",1.2078769595261523,0.9401113447793931,0.35789438,"collagen catabolic process"),
c("GO:0019557","L-histidine catabolic process to glutamate and formate",0.039701622981672124,0.9347706779347026,0.27965609,"collagen catabolic process"),
c("GO:0032781","positive regulation of ATP-dependent activity",0.05862204305979448,0.9112339323460823,-0,"positive regulation of ATP-dependent activity"),
c("GO:0002230","positive regulation of defense response to virus by host",0.0059047857119052345,0.8485241115131007,0.10684169,"positive regulation of ATP-dependent activity"),
c("GO:0008285","negative regulation of cell population proliferation",0.12504657495107693,0.8692013488323905,0.39868766,"positive regulation of ATP-dependent activity"),
c("GO:0010543","regulation of platelet activation",0.009264532480038059,0.865379762433177,0.11746209,"positive regulation of ATP-dependent activity"),
c("GO:0030155","regulation of cell adhesion",0.18103497036112626,0.8963597632605241,0.14784165,"positive regulation of ATP-dependent activity"),
c("GO:0030334","regulation of cell migration",0.2739092009852698,0.8933340427821649,0.17159261,"positive regulation of ATP-dependent activity"),
c("GO:0031630","regulation of synaptic vesicle fusion to presynaptic active zone membrane",0.0025868819438150904,0.9051005857008574,0.10194411,"positive regulation of ATP-dependent activity"),
c("GO:0032224","positive regulation of synaptic transmission, cholinergic",0.0036772612978684532,0.8540990141381414,0.10397527,"positive regulation of ATP-dependent activity"),
c("GO:0042177","negative regulation of protein catabolic process",0.034412175505688654,0.8830595453169198,0.1883348,"positive regulation of ATP-dependent activity"),
c("GO:0042391","regulation of membrane potential",0.3881799727488847,0.8853075741683958,0.15468386,"positive regulation of ATP-dependent activity"),
c("GO:0048742","regulation of skeletal muscle fiber development",0.008171691773041008,0.8857975980379102,0.11664158,"positive regulation of ATP-dependent activity"),
c("GO:2000045","regulation of G1/S transition of mitotic cell cycle",0.03849556003926547,0.9009606714519047,0.12765205,"positive regulation of ATP-dependent activity"),
c("GO:2000377","regulation of reactive oxygen species metabolic process",0.019230550549026927,0.9020614646495709,0.11472706,"positive regulation of ATP-dependent activity"),
c("GO:0033058","directional locomotion",0.00033720535328512596,0.9966640357510198,0,"directional locomotion"),
c("GO:0042590","antigen processing and presentation of exogenous peptide antigen via MHC class I",0.002148761119838795,0.925557238209291,-0,"antigen processing and presentation of exogenous peptide antigen via MHC class I"),
c("GO:0044281","small molecule metabolic process",14.792283225323425,0.9780663948484096,0.06195454,"small molecule metabolic process"),
c("GO:0044237","cellular metabolic process",39.640441131550894,0.9607218374798676,0.16054973,"small molecule metabolic process"),
c("GO:0050896","response to stimulus",16.94347616187058,1,-0,"response to stimulus"),
c("GO:0051156","glucose 6-phosphate metabolic process",0.2974889621857917,0.935181529032261,0.099233,"glucose 6-phosphate metabolic process"),
c("GO:0006468","protein phosphorylation",3.756826993126081,0.9107092662646804,0.39688397,"glucose 6-phosphate metabolic process"),
c("GO:0046491","L-methylmalonyl-CoA metabolic process",0.03570438580112436,0.9249819049301151,0.36491916,"glucose 6-phosphate metabolic process"),
c("GO:0051602","response to electrical stimulus",0.0013192851778162592,0.942417876173021,-0,"response to electrical stimulus"),
c("GO:0007187","G protein-coupled receptor signaling pathway, coupled to cyclic nucleotide second messenger",0.046172519870625386,0.838383141247805,0.25849792,"response to electrical stimulus"),
c("GO:0007267","cell-cell signaling",0.6278049885815375,0.8371071845164769,0.3400966,"response to electrical stimulus"),
c("GO:0009611","response to wounding",0.15748228404298487,0.9164570826057615,0.18567426,"response to electrical stimulus"),
c("GO:0010045","response to nickel cation",0.008036317361138221,0.9274653235178735,0.14270573,"response to electrical stimulus"),
c("GO:0038061","non-canonical NF-kappaB signal transduction",0.00812738742005464,0.8429625836578812,0.18580609,"response to electrical stimulus"),
c("GO:0048011","neurotrophin TRK receptor signaling pathway",0.004853788004950864,0.8399731631879044,0.22840628,"response to electrical stimulus"),
c("GO:0071230","cellular response to amino acid stimulus",0.04028742498226964,0.9089415214468597,0.34913084,"response to electrical stimulus"),
c("GO:0061384","heart trabecula morphogenesis",0.004130150239506871,0.8906173398189734,-0,"heart trabecula morphogenesis"),
c("GO:0008544","epidermis development",0.07531247737093637,0.8709912413869323,0.38291575,"heart trabecula morphogenesis"),
c("GO:0014009","glial cell proliferation",0.0023505920612211335,0.8130917455846588,0.38223189,"heart trabecula morphogenesis"),
c("GO:0016188","synaptic vesicle maturation",0.011329607599791493,0.859111055395595,0.38547788,"heart trabecula morphogenesis"),
c("GO:0048058","compound eye corneal lens development",3.692029415530576E-05,0.8568773299784243,0.26575519,"heart trabecula morphogenesis"),
c("GO:0048755","branching morphogenesis of a nerve",0.00046765705930053964,0.8399347988068976,0.37188298,"heart trabecula morphogenesis"),
c("GO:0072537","fibroblast activation",0.0003962778239336152,0.8732106214112562,0.00878759,"fibroblast activation"),
c("GO:0001816","cytokine production",0.001796787648891547,0.8435907774693517,0.28836,"fibroblast activation"),
c("GO:0003016","respiratory system process",0.002513041355504479,0.8664920716873523,0.39849793,"fibroblast activation"),
c("GO:0007625","grooming behavior",0.0011962175306319065,0.8697711611714057,0.2700237,"fibroblast activation"),
c("GO:0032610","interleukin-1 alpha production",0.00011371450599834175,0.8628245668912423,0.24570163,"fibroblast activation"),
c("GO:0050975","sensory perception of touch",0.0008294759420225361,0.8693193475847263,0.27868324,"fibroblast activation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "uniqueness",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()
