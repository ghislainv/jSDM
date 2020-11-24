\name{birds}
\alias{birds}
\docType{data}
\title{birds dataset}

\description{
The Swiss breeding bird survey ("Monitoring Häufige Brutvögel" MHB) has monitored the populations of 158 common species since 1999. The MHB sample consists of 267 1-km squares that are laid out as a grid across Switzerland. Fieldwork is conducted by about 200 skilled birdwatchers, most of them volunteers. Avian populations are monitored using a simplified territory mapping protocol, where each square is surveyed up to three times during the breeding season (only twice above the tree line). Surveys are conducted along a transect that does not change over the years. 

The \code{birds} dataset has the data for 2014, except one quadrat not surveyed in 2014. 
}


\usage{data("birds")}
\format{ A data frame with 266 observations on the following 166 variables.

158 bird species named in latin and whose occurrences are indicated as the number of visits to each site during which the species was observed, including 13 species not recorded in the year 2014 and 5 covariates collected on the 266 1x1 km quadrat as well as their identifiers and coordinates : 
    
    \tabular{ll}{
    \bold{\code{siteID}} \tab an alphanumeric site identifier \cr
    \tab \cr
    \bold{\code{coordx}} \tab a numeric vector indicating the x coordinate of the centre of the quadrat; the coordinate reference system intentionally not specified. \cr
    \tab \cr
    \bold{\code{coordy}} \tab a numeric vector indicating the y coordinate of the centre of the quadrat.\cr
        \tab  \cr
    \bold{\code{elev}} \tab a numeric vector indicating the mean elevation of the quadrat (m).\cr
        \tab  \cr
    \bold{\code{rlength}} \tab  the length of the route walked in the quadrat (km).\cr
    \tab  \cr
    \bold{\code{nsurvey}} \tab  a numeric vector indicating the number of replicate surveys planned in the quadrat; above the tree-line 2, otherwise 3. \cr 
            \tab \cr
    \bold{\code{forest}} \tab a numeric vector indicating the percentage of forest cover in the quadrat.\cr 
    \tab \cr 
    \bold{\code{obs14}} \tab a categorical vector indicating the identifying number of the observer.\cr 
}
}

\source{Swiss Ornithological Institute}

\details{
\tabular{ll}{
english name \tab latin name \cr 
Little Grebe \tab Tachybaptus ruficollis \cr                 
Great Crested Grebe \tab Podiceps cristatus \cr             
Grey Heron \tab Ardea cinerea \cr                           
Little Bittern \tab Ixobrychus minutus \cr                  
White Stork \tab Ciconia ciconia \cr                        
Mute Swan \tab Cygnus olor \cr                              
Greylag Goose \tab Anser anser \cr                          
Mallard \tab Anas platyrhynchos \cr                         
Common Merganser \tab Mergus merganser \cr                  
European Honey Buzzard \tab Pernis apivorus \cr             
Red Kite \tab Milvus milvus \cr                             
Black Kite \tab Milvus migrans \cr                          
Northern Goshawk \tab Accipiter gentilis \cr                
Eurasian Sparrowhawk \tab Accipiter nisus \cr               
Common Buzzard \tab Buteo buteo \cr                         
Golden Eagle \tab Aquila chrysaetos \cr                     
Peregrine Falcon \tab Falco peregrinus \cr                  
Eurasian Hobby \tab Falco subbuteo \cr                      
Common Kestrel \tab Falco tinnunculus \cr                   
Western Capercaillie \tab Tetrao urogallus \cr              
Black Grouse \tab Tetrao tetrix \cr                         
Rock Ptarmigan \tab Lagopus muta \cr                        
Hazel Grouse \tab Bonasa bonasia \cr                        
Rock Partridge \tab Alectoris graeca \cr                    
Common Quail \tab Coturnix coturnix \cr                     
Common Pheasant \tab Phasianus colchicus \cr                
Water Rail \tab Rallus aquaticus \cr                        
Common Moorhen \tab Gallinula chloropus \cr                 
Eurasian Coot \tab Fulica atra \cr                          
Northern Lapwing \tab Vanellus vanellus \cr                 
Little Ringed Plover \tab Charadrius dubius \cr             
Common Sandpiper \tab Actitis hypoleucos \cr                
Eurasian Woodcock \tab Scolopax rusticola \cr               
Yellow-legged Gull \tab Larus michahellis \cr               
Feral Pigeon \tab Columba livia domestica \cr               
Stock Dove \tab Columba oenas \cr                           
Common Wood Pigeon \tab Columba palumbus \cr                
European Turtle Dove \tab Streptopelia turtur \cr           
Eurasian Collared Dove \tab Streptopelia decaocto \cr       
Common Cuckoo \tab Cuculus canorus \cr                      
Barn Owl \tab Tyto alba \cr                                 
Scops Owl \tab Otus scops \cr                               
Eurasian Eagle-Owl \tab Bubo bubo \cr                       
Eurasian Pygmy Owl \tab Glaucidium passerinum \cr           
Little Owl \tab Athene noctua \cr                           
Tawny Owl \tab Strix aluco \cr                              
Long-eared Owl \tab Asio otus \cr                           
Boreal Owl \tab Aegolius funereus \cr                       
Alpine Swift \tab Apus melba \cr                            
Common Swift \tab Apus apus \cr                             
Common Kingfisher \tab Alcedo atthis \cr                    
Eurasian Hoopoe \tab Upupa epops \cr                        
Eurasian Wryneck \tab Jynx torquilla \cr                    
European Green Woodpecker \tab Picus viridis \cr            
Grey-headed Woodpecker \tab Picus canus \cr                 
Black Woodpecker \tab Dryocopus martius \cr                 
Great Spotted Woodpecker \tab Dendrocopos major \cr         
Middle Spotted Woodpecker \tab Dendrocopos medius \cr       
White-backed Woodpecker \tab Dendrocopos leucotos \cr       
Lesser Spotted Woodpecker \tab Dendrocopos minor \cr        
Eurasian Three-toed Woodpecker \tab Picoides tridactylus \cr
Woodlark \tab Lullula arborea \cr                         
Eurasian Skylark \tab Alauda arvensis \cr                   
Barn Swallow \tab Hirundo rustica \cr                       
Eurasian Crag Martin \tab Ptyonoprogne rupestris \cr        
Common House Martin \tab Delichon urbicum \cr               
Eurasian Golden Oriole \tab Oriolus oriolus \cr             
Northern Raven \tab Corvus corax \cr                        
Carrion/Hooded Crow \tab Corvus corone \cr                  
Rook \tab Corvus frugilegus \cr                             
Western Jackdaw \tab Corvus monedula \cr                    
Eurasian Magpie \tab Pica pica \cr                          
Spotted Nutcracker \tab Nucifraga caryocatactes \cr         
Eurasian Jay \tab Garrulus glandarius \cr                   
Red-billed Chough \tab Pyrrhocorax pyrrhocorax \cr          
Alpine Chough \tab Pyrrhocorax graculus \cr                 
Great Tit \tab Parus major \cr                              
Blue Tit \tab Parus caeruleus \cr                           
Coal Tit \tab Parus ater \cr                                
European Crested Tit \tab Parus cristatus \cr               
Marsh Tit \tab Parus palustris \cr                          
Willow Tit \tab Parus montanus \cr                          
Long-tailed Bushtit \tab Aegithalos caudatus \cr            
Eurasian Nuthatch \tab Sitta europaea \cr                   
Eurasian Treecreeper \tab Certhia familiaris \cr            
Short-toed Treecreeper \tab Certhia brachydactyla \cr       
Wallcreeper \tab Tichodroma muraria \cr                     
White-throated Dipper \tab Cinclus cinclus \cr              
Winter Wren \tab Troglodytes troglodytes \cr                
European Robin \tab Erithacus rubecula \cr                  
Common Nightingale \tab Luscinia megarhynchos \cr           
Bluethroat \tab Luscinia svecica \cr                        
Black Redstart \tab Phoenicurus ochruros \cr                
Common Redstart \tab Phoenicurus phoenicurus \cr            
Whinchat \tab Saxicola rubetra \cr                          
Eurasian Stonechat \tab Saxicola torquatus \cr              
Northern Wheatear \tab Oenanthe oenanthe \cr                
Rufous-tailed Rock Thrush \tab Monticola saxatilis \cr      
Ring Ouzel \tab Turdus torquatus \cr                        
Common Blackbird \tab Turdus merula \cr                     
Fieldfare \tab Turdus pilaris \cr                           
Song Thrush \tab Turdus philomelos \cr                      
Mistle Thrush \tab Turdus viscivorus \cr                    
Common Grasshopper Warbler \tab Locustella naevia \cr       
Savi's Warbler \tab Locustella luscinioides \cr             
Great Reed Warbler \tab Acrocephalus arundinaceus \cr       
Eurasian Reed Warbler \tab Acrocephalus scirpaceus \cr      
Marsh Warbler \tab Acrocephalus palustris \cr               
Melodious Warbler \tab Hippolais polyglotta \cr             
Icterine Warbler \tab Hippolais icterina \cr                
Eurasian Blackcap \tab Sylvia atricapilla \cr               
Garden Warbler \tab Sylvia borin \cr                      
Common Whitethroat \tab Sylvia communis \cr                 
Lesser Whitethroat \tab Sylvia curruca \cr                  
Willow Warbler \tab Phylloscopus trochilus \cr              
Common Chiffchaff \tab Phylloscopus collybita \cr           
Western Bonelli's Warbler \tab Phylloscopus bonelli \cr     
Wood Warbler \tab Phylloscopus sibilatrix \cr               
Goldcrest \tab Regulus regulus \cr                          
Firecrest \tab Regulus ignicapilla \cr                      
Spotted Flycatcher \tab Muscicapa striata \cr               
European Pied Flycatcher \tab Ficedula hypoleuca \cr        
Collared Flycatcher \tab Ficedula albicollis \cr            
Red-breasted Flycatcher \tab Ficedula parva \cr             
Dunnock \tab Prunella modularis \cr                         
Alpine Accentor \tab Prunella collaris \cr                  
Meadow Pipit \tab Anthus pratensis \cr                      
Tawny Pipit \tab Anthus campestris \cr                      
Tree Pipit \tab Anthus trivialis \cr                        
Water Pipit \tab Anthus spinoletta \cr                      
White Wagtail \tab Motacilla alba \cr                       
Grey Wagtail \tab Motacilla cinerea \cr                     
Western Yellow Wagtail \tab Motacilla flava \cr             
Woodchat Shrike \tab Lanius senator \cr                     
Red-backed Shrike \tab Lanius collurio \cr                  
Common Starling \tab Sturnus vulgaris \cr                   
House Sparrow \tab Passer domesticus \cr                    
Italian Sparrow \tab Passer hispaniolensis italiae \cr      
Eurasian Tree Sparrow \tab Passer montanus \cr              
White-winged Snowfinch \tab Montifringilla nivalis \cr      
Hawfinch \tab Coccothraustes coccothraustes \cr             
European Greenfinch \tab Carduelis chloris \cr              
European Goldfinch \tab Carduelis carduelis \cr             
Eurasian Siskin \tab Carduelis spinus \cr                   
Common Linnet \tab Carduelis cannabina \cr                  
Lesser Redpoll \tab Carduelis cabaret \cr                   
Citril Finch \tab Serinus citrinella \cr                    
European Serin \tab Serinus serinus \cr                     
Eurasian Bullfinch \tab Pyrrhula pyrrhula \cr               
Common Rosefinch \tab Carpodacus erythrinus \cr             
Red Crossbill \tab Loxia curvirostra \cr                    
Common Chaffinch \tab Fringilla coelebs \cr                 
Corn Bunting \tab Emberiza calandra \cr                     
Yellowhammer \tab Emberiza citrinella \cr                   
Cirl Bunting \tab Emberiza cirlus \cr                       
Ortolan Bunting \tab Emberiza hortulana \cr                 
Rock Bunting \tab Emberiza cia \cr                          
Common Reed Bunting \tab Emberiza schoeniclus \cr   
}
}
\references{
\tabular{l}{
Kéry and Royle (2016) Applied Hierarachical Modeling in Ecology Section 11.3 \cr
}}


\examples{
data(birds, package="jSDM")
head(birds)
# find species not recorded in 2014
which(colSums(birds[,1:158])==0)
}

\keyword{datasets}